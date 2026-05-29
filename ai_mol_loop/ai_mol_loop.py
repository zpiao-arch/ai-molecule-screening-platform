#!/usr/bin/env python3
"""A lightweight CLI for an AI molecule design closed-loop workflow.

Stages 1-3 define the file contracts between generation, filtering, scoring,
ranking, and feedback. Stage 4 uses real local chemistry libraries when they
are installed, with RDKit as the first supported backend for descriptors,
fingerprints, scaffolds, and ligand SDF export.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
import site
import sysconfig


DEFAULT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DESIGN_REPOS = DEFAULT_ROOT / "ai_molecule_design_repos"
DEFAULT_EVAL_TOOLS = DEFAULT_ROOT / "ai_drug_eval_tools"
DEFAULT_TARGET_CATALOG = Path(__file__).resolve().parent / "targets" / "influenza" / "target_catalog.json"
DEFAULT_EVIDENCE_DIR = Path(__file__).resolve().parent / "targets" / "influenza" / "evidence"
DEFAULT_EVIDENCE_SOURCES = Path(__file__).resolve().parent / "targets" / "influenza" / "evidence_sources.json"
DEFAULT_KNOWN_DRUGS = Path(__file__).resolve().parent / "targets" / "influenza" / "known_drugs.csv"
DEFAULT_PDB_STRUCTURES = Path(__file__).resolve().parent / "targets" / "influenza" / "pdb_structures.csv"

DEFAULT_WEIGHTS = {
    "validity": 0.10,
    "qed": 0.18,
    "sa": 0.14,
    "lipinski": 0.10,
    "docking": 0.32,
    "pose": 0.10,
    "novelty": 0.06,
}

DEFAULT_SEEDS = [
    ("seed_001", "CCO", "ethanol-like small polar fragment"),
    ("seed_002", "c1ccccc1", "benzene aromatic core"),
    ("seed_003", "CC(=O)Oc1ccccc1C(=O)O", "aspirin-like reference"),
    ("seed_004", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "caffeine-like reference"),
    ("seed_005", "CCN(CC)CC", "tertiary amine fragment"),
]

SUBSTITUENTS = [
    "C",
    "CC",
    "CCC",
    "O",
    "N",
    "F",
    "Cl",
    "C#N",
    "C(=O)O",
    "C(=O)N",
    "S(=O)(=O)N",
    "CO",
    "OC",
    "N1CCOCC1",
]

TEMPLATE_MOLECULES = [
    "c1ccccc1F",
    "c1ccccc1Cl",
    "c1ccccc1C#N",
    "c1ccccc1C(=O)O",
    "c1ccccc1C(=O)N",
    "CCOC(=O)c1ccccc1",
    "CCN(CC)C(=O)c1ccccc1",
    "COc1ccccc1O",
    "CC(C)NC(=O)c1ccccc1",
    "O=C(NCCO)c1ccccc1",
    "CCS(=O)(=O)Nc1ccccc1",
    "CC(C)Oc1ccccc1C#N",
]


def clamp(value: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, value))


def stable_hash(text: str) -> int:
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return int(digest[:16], 16)


def ensure_project_dirs(project: Path) -> None:
    for name in [
        "candidates",
        "scores",
        "ranked",
        "feedback",
        "reports",
        "seeds",
        "external",
        "runs",
        "briefs",
        "prompts",
        "targets",
        "evidence",
        "filtered",
        "stage3",
        "stage4",
        "stage4_5",
        "stage4_6",
        "stage5",
        "stage6",
        "stage7",
    ]:
        (project / name).mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: Sequence[Dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields))
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_json(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, data: Dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def project_path(raw: str) -> Path:
    return Path(raw).expanduser().resolve()


def default_config() -> Dict[str, object]:
    return {
        "schema_version": "0.1",
        "target": {
            "name": "demo_target",
            "protein_pdb": "",
            "reference_ligand_sdf": "",
            "pocket": {
                "center": [0.0, 0.0, 0.0],
                "size": [20.0, 20.0, 20.0],
                "source": "manual_or_p2rank",
            },
        },
        "weights": dict(DEFAULT_WEIGHTS),
        "thresholds": {
            "advance_total": 0.64,
            "pose_min": 0.50,
            "max_candidates_per_round": 10000,
        },
        "tool_paths": {
            "REINVENT4": str(DEFAULT_DESIGN_REPOS / "REINVENT4"),
            "DrugEx": str(DEFAULT_DESIGN_REPOS / "DrugEx"),
            "GraphINVENT": str(DEFAULT_DESIGN_REPOS / "GraphINVENT"),
            "MolBART": str(DEFAULT_DESIGN_REPOS / "MolBART-master"),
            "guacamol": str(DEFAULT_DESIGN_REPOS / "guacamol-master"),
            "DockStream": str(DEFAULT_DESIGN_REPOS / "DockStream"),
            "AutoDock-Vina": str(DEFAULT_EVAL_TOOLS / "AutoDock-Vina"),
            "gnina": str(DEFAULT_EVAL_TOOLS / "gnina"),
            "posebusters": str(DEFAULT_EVAL_TOOLS / "posebusters"),
            "openfe": str(DEFAULT_EVAL_TOOLS / "openfe"),
            "p2rank": str(DEFAULT_EVAL_TOOLS / "p2rank"),
        },
        "notes": [
            "Proxy scoring is only for wiring the closed loop.",
            "Replace proxy docking with Vina/GNINA/DockStream CSV outputs before making scientific claims.",
        ],
    }


def explain() -> None:
    print(
        """
AI分子设计闭环可以拆成 6 个可替换模块：

1. 定义任务
   输入靶点结构、口袋位置、参考配体、期望性质和约束。

2. 生成候选分子
   REINVENT4 / DrugEx / GraphINVENT / MolBART 负责生成 SMILES 或结构库。
   本 CLI 先用轻量模板生成器占位，真实工具输出只要整理成 candidates CSV 即可接入。

3. 快速过滤
   去重、基础 SMILES 合法性、类药性、合成可及性、简单理化性质。

4. 结构验证与评分
   RDKit 生成真实描述符、指纹、骨架、3D SDF；Vina/GNINA 可选对接打分。
   PoseBusters 检查 docked pose 是否物理合理。
   高价值 top hits 可进入 OpenFE 自由能计算。

5. 决策排序
   把 docking、pose、类药性、合成可及性、新颖性合成一个多目标分数。

6. 反馈下一轮
   top 分子作为下一轮 seed、强化学习奖励或微调数据，继续生成-验证-筛选。

CLI 文件流：
  candidates/round_N_candidates.csv
  scores/round_N_scores.csv
  ranked/round_N_ranked.csv
  feedback/round_N_feedback.json
  seeds/round_N_seeds.csv

重要边界：
  proxy score 只用于打通工程闭环，不等于真实药效证明。
  Stage 4 的 RDKit 输出是真实化学库计算结果，但仍不是药效证明。
  真实项目里，Vina/GNINA/PoseBusters/OpenFE 的结果要逐级替换 proxy 字段。
""".strip()
    )


def init_project(args: argparse.Namespace) -> None:
    project = project_path(args.project)
    ensure_project_dirs(project)
    config_path = project / "config.json"
    if config_path.exists() and not args.force:
        raise SystemExit(f"Config already exists: {config_path}. Use --force to overwrite.")

    write_json(config_path, default_config())
    seed_rows = [
        {"id": seed_id, "smiles": smiles, "note": note}
        for seed_id, smiles, note in DEFAULT_SEEDS
    ]
    write_csv(project / "seeds" / "round_0_seeds.csv", seed_rows, ["id", "smiles", "note"])
    write_text(
        project / "README.md",
        f"""# AI Molecule Closed Loop Project

This directory was initialized by `ai_mol_loop.py`.

Run a demo:

```bash
python3 {Path(__file__).resolve()} run-demo {project} --rounds 2 --n 24 --top 6
```

Main files:

- `config.json`: target, weights, local tool paths.
- `seeds/`: molecules used to start each round.
- `candidates/`: generated candidate SMILES.
- `scores/`: proxy or imported docking/scoring results.
- `ranked/`: selected and ranked candidates.
- `feedback/`: next-round feedback package.
""",
    )
    print(f"Initialized project: {project}")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def load_target_catalog(path: Optional[str] = None) -> Dict[str, object]:
    catalog_path = Path(path).expanduser().resolve() if path else DEFAULT_TARGET_CATALOG
    if not catalog_path.exists():
        raise SystemExit(f"Missing target catalog: {catalog_path}")
    return read_json(catalog_path)


def target_total_score(target: Dict[str, object], weights: Dict[str, object]) -> float:
    scores = target.get("scores", {})
    if not isinstance(scores, dict):
        return 0.0
    total = 0.0
    total_weight = 0.0
    for name, raw_weight in weights.items():
        try:
            weight = float(raw_weight)
            score = float(scores.get(name, 0.0))
        except (TypeError, ValueError):
            continue
        total += weight * score
        total_weight += weight
    if total_weight <= 0:
        return 0.0
    return round(total / total_weight, 4)


def normalize_query(text: str) -> str:
    mapping = {
        "甲流": "influenza_a",
        "乙流": "influenza_b",
        "流感": "influenza",
        "a型流感": "influenza_a",
        "b型流感": "influenza_b",
        "h1n1": "h1n1",
        "h3n2": "h3n2",
    }
    lowered = text.strip().lower().replace(" ", "_")
    return mapping.get(lowered, lowered)


def target_matches_query(target: Dict[str, object], query: str) -> bool:
    normalized = normalize_query(query)
    if not normalized or normalized in {"all", "influenza"}:
        return True
    scope = target.get("pathogen_scope", [])
    fields = [
        str(target.get("id", "")),
        str(target.get("display_name", "")),
        str(target.get("short_name", "")),
        str(target.get("target_family", "")),
        str(target.get("recommendation", "")),
    ]
    if isinstance(scope, list):
        fields.extend(str(item) for item in scope)
    haystack = " ".join(fields).lower().replace(" ", "_")
    return normalized in haystack


def sorted_targets(catalog: Dict[str, object], query: str) -> List[Dict[str, object]]:
    weights = catalog.get("scoring_weights", {})
    if not isinstance(weights, dict):
        weights = {}
    targets = catalog.get("targets", [])
    if not isinstance(targets, list):
        return []
    selected = [target for target in targets if isinstance(target, dict) and target_matches_query(target, query)]
    for target in selected:
        target["_target_score"] = target_total_score(target, weights)
    selected.sort(key=lambda target: float(target.get("_target_score", 0.0)), reverse=True)
    return selected


def render_target_selection_report(
    query: str,
    targets: Sequence[Dict[str, object]],
    catalog: Dict[str, object],
    evidence_map: Optional[Dict[str, Dict[str, str]]] = None,
) -> str:
    evidence_map = evidence_map or {}
    lines = [
        "# Target Selection Report",
        "",
        f"- Query: `{query}`",
        f"- Catalog domain: `{catalog.get('domain', 'unknown')}`",
        f"- Catalog updated: `{catalog.get('updated', 'unknown')}`",
        "",
        "## Ranked Targets",
        "",
        "| rank | target id | target | score | evidence score | readiness | recommendation | PDB evidence | PubMed evidence | known drugs |",
        "|---:|---|---|---:|---:|---|---|---:|---:|---|",
    ]
    for idx, target in enumerate(targets, start=1):
        drugs = target.get("known_drugs", [])
        if isinstance(drugs, list):
            drug_text = ", ".join(str(item) for item in drugs)
        else:
            drug_text = str(drugs)
        evidence = evidence_map.get(str(target.get("id", "")), {})
        lines.append(
            f"| {idx} | `{target.get('id', '')}` | {target.get('display_name', '')} | "
            f"{float(target.get('_target_score', 0.0)):.3f} | {evidence.get('evidence_score', '')} | "
            f"{evidence.get('readiness', '')} | {target.get('recommendation', '')} | "
            f"{evidence.get('pdb_entries', '')} | {evidence.get('pubmed_articles', '')} | {drug_text} |"
        )

    if targets:
        top = targets[0]
        lines.extend(
            [
                "",
                "## Recommended First Target",
                "",
                f"- Target: `{top.get('id')}` ({top.get('display_name')})",
                f"- Score: {float(top.get('_target_score', 0.0)):.3f}",
                f"- Reason: {top.get('recommendation_reason', '')}",
                "",
                "## Why This Matters",
                "",
                "This module chooses a computationally practical and clinically grounded target before molecule generation. "
                "It prevents the generator from producing molecules for an undefined or poorly validated pocket.",
                "",
                "## Caveats",
                "",
                "- Target selection is a prioritization step, not a therapeutic claim.",
                "- Docking and proxy scores are not measured potency.",
                "- Evidence counts are metadata support counts, not proof of efficacy.",
                "- M2 is retained as a historical control because adamantane resistance makes it unsuitable as the primary MVP target.",
                "- PA endonuclease targets need metal-aware receptor and ligand preparation.",
            ]
        )
    return "\n".join(lines) + "\n"


def evidence_summary_by_target(path: Optional[Path] = None) -> Dict[str, Dict[str, str]]:
    summary_path = path or (DEFAULT_EVIDENCE_DIR / "evidence_summary.csv")
    rows = read_csv(summary_path)
    return {row.get("target_id", ""): row for row in rows if row.get("target_id")}


def target_select(args: argparse.Namespace) -> None:
    project = project_path(args.project) if args.project else None
    if project:
        ensure_project_dirs(project)
    catalog = load_target_catalog(args.catalog)
    targets = sorted_targets(catalog, args.disease)
    if args.target:
        wanted = args.target.lower()
        targets = [
            target
            for target in targets
            if wanted in str(target.get("id", "")).lower()
            or wanted in str(target.get("display_name", "")).lower()
            or wanted in str(target.get("short_name", "")).lower()
        ]
    if not targets:
        raise SystemExit(f"No targets matched disease/query: {args.disease}")

    fields = [
        "rank",
        "target_id",
        "display_name",
        "score",
        "recommendation",
        "pdb_evidence_count",
        "pubmed_evidence_count",
        "evidence_score",
        "readiness",
        "evidence_report",
        "known_drugs",
        "recommended_pdb",
        "recommendation_reason",
    ]
    rows = []
    evidence_map = evidence_summary_by_target()
    for idx, target in enumerate(targets[: args.top], start=1):
        structures = target.get("representative_structures", [])
        recommended_pdb = ""
        if isinstance(structures, list) and structures:
            recommended_pdb = str(structures[0].get("pdb_id", ""))
        drugs = target.get("known_drugs", [])
        evidence = evidence_map.get(str(target.get("id", "")), {})
        evidence_dir = evidence.get("evidence_dir", "")
        rows.append(
            {
                "rank": idx,
                "target_id": target.get("id", ""),
                "display_name": target.get("display_name", ""),
                "score": float(target.get("_target_score", 0.0)),
                "recommendation": target.get("recommendation", ""),
                "pdb_evidence_count": evidence.get("pdb_entries", ""),
                "pubmed_evidence_count": evidence.get("pubmed_articles", ""),
                "evidence_score": evidence.get("evidence_score", ""),
                "readiness": evidence.get("readiness", ""),
                "evidence_report": str(Path(evidence_dir) / "evidence_report.md") if evidence_dir else "",
                "known_drugs": "; ".join(str(item) for item in drugs) if isinstance(drugs, list) else str(drugs),
                "recommended_pdb": evidence.get("primary_pdb", "") or recommended_pdb,
                "recommendation_reason": target.get("recommendation_reason", ""),
            }
        )

    if project:
        output_csv = project / "targets" / "target_selection.csv"
        output_report = project / "reports" / "target_selection.md"
        write_csv(output_csv, rows, fields)
        write_text(output_report, render_target_selection_report(args.disease, targets[: args.top], catalog, evidence_map))
        print(f"Wrote target ranking: {output_csv}")
        print(f"Wrote target report: {output_report}")

    print("Ranked targets:")
    for row in rows:
        print(
            f"{row['rank']}. {row['target_id']} | score={row['score']:.3f} | "
            f"{row['recommendation']} | drugs={row['known_drugs']}"
        )


def get_target_by_id(catalog: Dict[str, object], target_id: str) -> Dict[str, object]:
    targets = catalog.get("targets", [])
    if not isinstance(targets, list):
        raise SystemExit("Target catalog has no targets list.")
    for target in targets:
        if isinstance(target, dict) and target.get("id") == target_id:
            return target
    raise SystemExit(f"Target not found: {target_id}")


def target_to_brief(
    target: Dict[str, object],
    disease: str,
    free_text: str,
    max_heavy_atoms: int,
    max_molecular_weight: float,
) -> Dict[str, object]:
    structures = target.get("representative_structures", [])
    first_structure = structures[0] if isinstance(structures, list) and structures else {}
    if not isinstance(first_structure, dict):
        first_structure = {}
    site = target.get("binding_site", {})
    if not isinstance(site, dict):
        site = {}
    guidance = target.get("generation_guidance", {})
    if not isinstance(guidance, dict):
        guidance = {}
    validation = target.get("validation_plan", {})
    if not isinstance(validation, dict):
        validation = {}
    known_drugs = target.get("known_drugs", [])

    return {
        "schema_version": "0.1",
        "target_catalog_id": target.get("id", ""),
        "target_name": target.get("display_name", ""),
        "disease_context": disease,
        "target_rationale": target.get("recommendation_reason", ""),
        "free_text_requirement": free_text
        or f"Generate diverse virtual-screening candidates for {target.get('display_name', '')}. Use known drugs as controls, not as claims of new efficacy.",
        "protein": {
            "name": target.get("display_name", ""),
            "gene": "",
            "pdb_id": first_structure.get("pdb_id", ""),
            "structure_file": "",
        },
        "binding_site": {
            "description": site.get("description", ""),
            "reference_ligand": site.get("reference_ligand", first_structure.get("ligand", "")),
            "key_residues": site.get("key_residues", []),
            "center": [],
            "size": [],
            "source": site.get("source", ""),
            "box_strategy": site.get("box_strategy", ""),
        },
        "design_intent": {
            "style": guidance.get("style", "target_guided_virtual_screening"),
            "primary_goal": guidance.get("primary_goal", ""),
            "desired_activity": "in silico enrichment and ranking, not measured potency.",
            "must_have": guidance.get("must_have", []),
            "avoid": guidance.get("avoid", []),
            "desired_properties": guidance.get("desired_properties", []),
            "selectivity_notes": guidance.get("selectivity_notes", ""),
        },
        "generation_constraints": {
            "max_heavy_atoms": max_heavy_atoms,
            "max_molecular_weight": max_molecular_weight,
            "prefer_synthesizable": True,
            "preserve_reference_scaffold": False,
            "diversity_required": True,
        },
        "validation_plan": {
            "first_pass": "RDKit/proxy filters, validity, QED, synthetic accessibility",
            "structure_scoring": "AutoDock Vina or GNINA docking score against the selected target pocket",
            "pose_quality": "PoseBusters pass/fail or pose quality score",
            "high_cost_followup": "OpenFE or higher-cost molecular simulation for top hits only",
            "computational_steps": validation.get("computational", []),
            "experimental_path": validation.get("experimental", []),
            "known_positive_controls": known_drugs if isinstance(known_drugs, list) else [],
        },
        "compliance_boundary": {
            "scope": "computational screening and research planning only",
            "no_claims": [
                "Do not claim therapeutic efficacy.",
                "Do not provide dosing, clinical advice, or experimental human-use claims.",
                "Do not present proxy or docking scores as measured potency.",
                "Do not treat known drugs as newly discovered molecules.",
                "Do not include synthesis procedures unless separately requested for benign feasibility review.",
            ],
        },
        "source_urls": target.get("source_urls", []),
    }


def brief_from_target(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    ensure_project_dirs(project)
    config_path = project / "config.json"
    config = load_config(project) if config_path.exists() else default_config()
    catalog = load_target_catalog(args.catalog)
    target_id = args.target
    if not target_id:
        ranked = sorted_targets(catalog, args.disease)
        if not ranked:
            raise SystemExit(f"No target candidates for disease/query: {args.disease}")
        target_id = str(ranked[0].get("id", ""))
    target = get_target_by_id(catalog, target_id)
    brief = target_to_brief(
        target,
        args.disease,
        args.free_text,
        args.max_heavy_atoms,
        args.max_molecular_weight,
    )

    brief_path = project / "briefs" / "target_brief.json"
    prompt_path = project / "prompts" / "generator_prompt.md"
    if (brief_path.exists() or prompt_path.exists()) and not args.force:
        raise SystemExit("Target brief or prompt already exists. Use --force to overwrite.")

    write_json(brief_path, brief)
    write_text(prompt_path, render_generator_prompt(brief))
    site = brief.get("binding_site", {})
    protein = brief.get("protein", {})
    if not isinstance(site, dict):
        site = {}
    if not isinstance(protein, dict):
        protein = {}
    config["target"] = {
        "name": brief.get("target_name", ""),
        "target_catalog_id": brief.get("target_catalog_id", ""),
        "protein_pdb": "",
        "reference_ligand_sdf": site.get("reference_ligand", ""),
        "pdb_id": protein.get("pdb_id", ""),
        "pocket": {
            "center": site.get("center", []),
            "size": site.get("size", []),
            "source": site.get("source", ""),
            "description": site.get("description", ""),
            "key_residues": site.get("key_residues", []),
            "box_strategy": site.get("box_strategy", ""),
        },
    }
    config["target_brief_file"] = str(brief_path)
    config["generator_prompt_file"] = str(prompt_path)
    write_json(project / "config.json", config)

    print(f"Selected target: {target_id}")
    print(f"Wrote target brief: {brief_path}")
    print(f"Wrote generator prompt: {prompt_path}")
    return prompt_path


def fetch_json(url: str, params: Optional[Dict[str, object]] = None, timeout: int = 20) -> Tuple[Optional[Dict[str, object]], Optional[str]]:
    if params:
        query = urllib.parse.urlencode(params, doseq=True)
        url = f"{url}?{query}"
    request = urllib.request.Request(url, headers={"User-Agent": "ai-mol-loop/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            raw = response.read().decode("utf-8")
        return json.loads(raw), None
    except Exception as exc:  # network and malformed upstream responses should not break offline evidence builds
        return None, str(exc)


def load_evidence_sources(path: Optional[str] = None) -> Dict[str, object]:
    source_path = Path(path).expanduser().resolve() if path else DEFAULT_EVIDENCE_SOURCES
    if not source_path.exists():
        raise SystemExit(f"Missing evidence sources: {source_path}")
    return read_json(source_path)


def catalog_targets_by_id(catalog: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    targets = catalog.get("targets", [])
    if not isinstance(targets, list):
        return {}
    return {
        str(target.get("id")): target
        for target in targets
        if isinstance(target, dict) and target.get("id")
    }


def target_ids_for_evidence(args: argparse.Namespace, catalog: Dict[str, object]) -> List[str]:
    by_id = catalog_targets_by_id(catalog)
    if args.target:
        if args.target not in by_id:
            raise SystemExit(f"Target not found in catalog: {args.target}")
        return [args.target]
    ranked = sorted_targets(catalog, args.disease)
    return [str(target.get("id")) for target in ranked if target.get("id")]


def fetch_rcsb_metadata(pdb_id: str, timeout: int) -> Tuple[Dict[str, object], Optional[str]]:
    url = f"https://data.rcsb.org/rest/v1/core/entry/{urllib.parse.quote(pdb_id)}"
    data, error = fetch_json(url, timeout=timeout)
    if error or not data:
        return {"pdb_id": pdb_id, "fetch_status": "failed"}, error

    citation = data.get("rcsb_primary_citation", {})
    if not isinstance(citation, dict):
        citation = {}
    entry_info = data.get("rcsb_entry_info", {})
    if not isinstance(entry_info, dict):
        entry_info = {}
    accession = data.get("rcsb_accession_info", {})
    if not isinstance(accession, dict):
        accession = {}
    struct = data.get("struct", {})
    if not isinstance(struct, dict):
        struct = {}
    exptl = data.get("exptl", [])
    method = ""
    if isinstance(exptl, list) and exptl and isinstance(exptl[0], dict):
        method = str(exptl[0].get("method", ""))

    resolution = entry_info.get("resolution_combined", [])
    if isinstance(resolution, list):
        resolution_value = resolution[0] if resolution else ""
    else:
        resolution_value = resolution

    return {
        "pdb_id": pdb_id,
        "fetch_status": "ok",
        "title": struct.get("title", ""),
        "experimental_method": method,
        "resolution": resolution_value,
        "release_date": accession.get("initial_release_date", ""),
        "primary_citation_title": citation.get("title", ""),
        "primary_citation_year": citation.get("year", ""),
        "primary_citation_pubmed_id": citation.get("pdbx_database_id_PubMed", ""),
        "rcsb_url": f"https://www.rcsb.org/structure/{pdb_id}",
    }, None


def pubmed_search(query: str, retmax: int, timeout: int) -> Tuple[List[str], Optional[str]]:
    data, error = fetch_json(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": retmax,
            "sort": "relevance",
        },
        timeout=timeout,
    )
    if error or not data:
        return [], error
    result = data.get("esearchresult", {})
    if not isinstance(result, dict):
        return [], "Unexpected PubMed esearch result shape"
    ids = result.get("idlist", [])
    if not isinstance(ids, list):
        return [], "Unexpected PubMed idlist shape"
    return [str(item) for item in ids], None


def pubmed_summary(pmids: Sequence[str], timeout: int) -> Tuple[List[Dict[str, object]], Optional[str]]:
    if not pmids:
        return [], None
    data, error = fetch_json(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "json",
        },
        timeout=timeout,
    )
    if error or not data:
        return [], error
    result = data.get("result", {})
    if not isinstance(result, dict):
        return [], "Unexpected PubMed esummary result shape"
    articles = []
    for pmid in pmids:
        item = result.get(pmid, {})
        if not isinstance(item, dict):
            continue
        authors = item.get("authors", [])
        first_author = ""
        if isinstance(authors, list) and authors and isinstance(authors[0], dict):
            first_author = str(authors[0].get("name", ""))
        articles.append(
            {
                "pmid": pmid,
                "title": item.get("title", ""),
                "journal": item.get("fulljournalname", item.get("source", "")),
                "pubdate": item.get("pubdate", ""),
                "first_author": first_author,
                "doi_or_elocation": item.get("elocationid", ""),
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }
        )
    return articles, None


def collect_pubmed_articles(queries: Sequence[str], retmax: int, timeout: int, offline: bool) -> Tuple[List[Dict[str, object]], List[Dict[str, str]]]:
    articles_by_pmid: Dict[str, Dict[str, object]] = {}
    errors: List[Dict[str, str]] = []
    if offline:
        return [], errors
    for query in queries:
        ids, error = pubmed_search(query, retmax, timeout)
        if error:
            errors.append({"source": "pubmed_esearch", "query": query, "error": error})
            continue
        time.sleep(0.34)
        summaries, summary_error = pubmed_summary(ids, timeout)
        if summary_error:
            errors.append({"source": "pubmed_esummary", "query": query, "error": summary_error})
            continue
        for article in summaries:
            article["query"] = query
            articles_by_pmid[str(article.get("pmid"))] = article
        time.sleep(0.34)
    return list(articles_by_pmid.values()), errors


def rows_for_target(rows: Sequence[Dict[str, str]], target_id: str) -> List[Dict[str, str]]:
    return [row for row in rows if row.get("target_id") == target_id]


def status_rows(rows: Sequence[Dict[str, str]], status: str) -> List[Dict[str, str]]:
    return [row for row in rows if row.get("status_for_workflow") == status]


def first_text(items: object) -> str:
    if isinstance(items, list):
        for item in items:
            if item:
                return str(item)
    if isinstance(items, str):
        return items
    return ""


def first_structure_asset(
    target: Dict[str, object],
    structure_rows: Sequence[Dict[str, str]],
    pdb_entries: Sequence[Dict[str, object]],
) -> Dict[str, object]:
    status_by_pdb = {
        str(entry.get("pdb_id", "")).upper(): entry
        for entry in pdb_entries
        if isinstance(entry, dict) and entry.get("pdb_id")
    }
    preferred_uses = [
        "recommended_demo_structure",
        "secondary_demo_structure",
        "secondary_b_flu_structure",
        "reference_control",
        "candidate_structure_verify_before_use",
        "historical_control_structure",
    ]
    sorted_rows = sorted(
        structure_rows,
        key=lambda row: preferred_uses.index(row.get("use", "zzzz")) if row.get("use") in preferred_uses else 99,
    )
    if sorted_rows:
        row = sorted_rows[0]
        pdb_id = row.get("pdb_id", "")
        metadata = status_by_pdb.get(pdb_id.upper(), {})
        return {
            "pdb_id": pdb_id,
            "ligand": row.get("ligand", ""),
            "use": row.get("use", ""),
            "description": row.get("description", ""),
            "source_url": row.get("source_url", ""),
            "metadata_status": metadata.get("fetch_status", ""),
            "resolution": metadata.get("resolution", ""),
            "rcsb_url": metadata.get("rcsb_url", row.get("source_url", "")),
        }

    structures = target.get("representative_structures", [])
    if isinstance(structures, list) and structures and isinstance(structures[0], dict):
        row = structures[0]
        pdb_id = str(row.get("pdb_id", ""))
        metadata = status_by_pdb.get(pdb_id.upper(), {})
        return {
            "pdb_id": pdb_id,
            "ligand": row.get("ligand", ""),
            "use": row.get("use", ""),
            "description": row.get("description", ""),
            "source_url": metadata.get("rcsb_url", f"https://www.rcsb.org/structure/{pdb_id}" if pdb_id else ""),
            "metadata_status": metadata.get("fetch_status", ""),
            "resolution": metadata.get("resolution", ""),
            "rcsb_url": metadata.get("rcsb_url", ""),
        }
    return {}


def assay_path_from_target(target: Dict[str, object]) -> str:
    validation = target.get("validation_plan", {})
    if not isinstance(validation, dict):
        return ""
    experimental = validation.get("experimental", [])
    if isinstance(experimental, list) and experimental:
        return " | ".join(str(item) for item in experimental)
    return first_text(experimental)


def build_closed_loop_assets(
    target: Dict[str, object],
    known_drug_rows: Sequence[Dict[str, str]],
    structure_rows: Sequence[Dict[str, str]],
    pdb_entries: Sequence[Dict[str, object]],
) -> Dict[str, object]:
    positive = status_rows(known_drug_rows, "positive_control")
    reference = status_rows(known_drug_rows, "reference_control")
    historical = status_rows(known_drug_rows, "historical_control_not_recommended")
    site = target.get("binding_site", {})
    if not isinstance(site, dict):
        site = {}
    validation = target.get("validation_plan", {})
    if not isinstance(validation, dict):
        validation = {}
    return {
        "primary_structure": first_structure_asset(target, structure_rows, pdb_entries),
        "positive_controls": [row.get("drug", "") for row in positive if row.get("drug")],
        "reference_controls": [row.get("drug", "") for row in reference if row.get("drug")],
        "historical_controls": [row.get("drug", "") for row in historical if row.get("drug")],
        "binding_site": {
            "source": site.get("source", ""),
            "description": site.get("description", ""),
            "reference_ligand": site.get("reference_ligand", ""),
            "key_residues": site.get("key_residues", []),
            "box_strategy": site.get("box_strategy", ""),
        },
        "computational_validation": validation.get("computational", []),
        "experimental_validation": validation.get("experimental", []),
        "assay_path": assay_path_from_target(target),
    }


def evidence_quality_score(
    target: Dict[str, object],
    catalog_weights: Dict[str, object],
    package: Dict[str, object],
    known_drug_rows: Sequence[Dict[str, str]],
) -> Dict[str, object]:
    pdb_entries = [entry for entry in package.get("pdb_entries", []) if isinstance(entry, dict)]
    pubmed_articles = [article for article in package.get("pubmed_articles", []) if isinstance(article, dict)]
    open_data_queries = package.get("open_data_queries", [])
    if not isinstance(open_data_queries, list):
        open_data_queries = []

    ok_pdb = sum(1 for entry in pdb_entries if str(entry.get("fetch_status", "")).lower() == "ok")
    queued_pdb = sum(1 for entry in pdb_entries if str(entry.get("fetch_status", "")).lower() == "queued_offline")
    positive_controls = len(status_rows(known_drug_rows, "positive_control"))
    reference_controls = len(status_rows(known_drug_rows, "reference_control"))
    historical_controls = len(status_rows(known_drug_rows, "historical_control_not_recommended"))
    has_assay = 1.0 if "assay" in assay_path_from_target(target).lower() else 0.0

    clinical_component = clamp((positive_controls + 0.5 * reference_controls) / 2.0)
    if historical_controls and not positive_controls:
        clinical_component = min(clinical_component, 0.25)
    structure_component = clamp((ok_pdb + 0.35 * queued_pdb) / 2.0)
    literature_component = clamp(len(pubmed_articles) / 6.0)
    open_data_component = clamp(len(open_data_queries) / 3.0)
    assay_component = has_assay
    catalog_component = target_total_score(target, catalog_weights) if catalog_weights else 0.0

    total = (
        clinical_component * 0.25
        + structure_component * 0.25
        + literature_component * 0.20
        + assay_component * 0.15
        + open_data_component * 0.10
        + catalog_component * 0.05
    )
    if total >= 0.80:
        readiness = "ready_for_mvp_closed_loop"
    elif total >= 0.65:
        readiness = "usable_with_caveats"
    elif total >= 0.45:
        readiness = "secondary_or_needs_manual_review"
    else:
        readiness = "not_recommended_for_mvp"

    return {
        "score": round(total, 4),
        "readiness": readiness,
        "components": {
            "clinical_controls": round(clinical_component, 4),
            "structure_metadata": round(structure_component, 4),
            "literature_metadata": round(literature_component, 4),
            "assay_path": round(assay_component, 4),
            "open_data_entry_points": round(open_data_component, 4),
            "catalog_prior": round(catalog_component, 4),
        },
        "counts": {
            "pdb_entries": len(pdb_entries),
            "ok_pdb_entries": ok_pdb,
            "queued_pdb_entries": queued_pdb,
            "pubmed_articles": len(pubmed_articles),
            "positive_controls": positive_controls,
            "reference_controls": reference_controls,
            "historical_controls": historical_controls,
            "open_data_queries": len(open_data_queries),
        },
        "interpretation": "Evidence score ranks source readiness for the computational closed loop; it is not a measured efficacy score.",
    }


def enrich_evidence_package(
    target: Dict[str, object],
    catalog_weights: Dict[str, object],
    package: Dict[str, object],
    known_drug_rows: Sequence[Dict[str, str]],
    structure_rows: Sequence[Dict[str, str]],
) -> Dict[str, object]:
    package["known_drug_rows"] = list(known_drug_rows)
    package["structure_rows"] = list(structure_rows)
    package["closed_loop_assets"] = build_closed_loop_assets(
        target,
        known_drug_rows,
        structure_rows,
        [entry for entry in package.get("pdb_entries", []) if isinstance(entry, dict)],
    )
    package["evidence_quality"] = evidence_quality_score(target, catalog_weights, package, known_drug_rows)
    package["source_traceability"] = {
        "target_catalog_id": target.get("id", ""),
        "target_catalog_recommendation": target.get("recommendation", ""),
        "known_drugs_file": str(DEFAULT_KNOWN_DRUGS),
        "pdb_structures_file": str(DEFAULT_PDB_STRUCTURES),
        "evidence_sources_file": str(DEFAULT_EVIDENCE_SOURCES),
    }
    return package


def existing_items_by_key(package: Dict[str, object], list_name: str, key: str) -> Dict[str, Dict[str, object]]:
    values = package.get(list_name, [])
    if not isinstance(values, list):
        return {}
    return {
        str(item.get(key, "")): item
        for item in values
        if isinstance(item, dict) and item.get(key)
    }


def load_existing_evidence_package(output_root: Path, target_id: str) -> Dict[str, object]:
    target_dir = output_root / target_id
    evidence_path = target_dir / "evidence.json"
    if evidence_path.exists():
        try:
            data = read_json(evidence_path)
            return data if isinstance(data, dict) else {}
        except (json.JSONDecodeError, OSError):
            pass
    fallback: Dict[str, object] = {}
    pdb_entries = read_csv(target_dir / "pdb_entries.csv")
    pubmed_articles = read_csv(target_dir / "pubmed_articles.csv")
    if pdb_entries:
        fallback["pdb_entries"] = pdb_entries
    if pubmed_articles:
        fallback["pubmed_articles"] = pubmed_articles
    return fallback


def render_evidence_report(
    target: Dict[str, object],
    package: Dict[str, object],
) -> str:
    quality = package.get("evidence_quality", {})
    if not isinstance(quality, dict):
        quality = {}
    components = quality.get("components", {})
    if not isinstance(components, dict):
        components = {}
    counts = quality.get("counts", {})
    if not isinstance(counts, dict):
        counts = {}
    assets = package.get("closed_loop_assets", {})
    if not isinstance(assets, dict):
        assets = {}
    primary_structure = assets.get("primary_structure", {})
    if not isinstance(primary_structure, dict):
        primary_structure = {}
    binding_site = assets.get("binding_site", {})
    if not isinstance(binding_site, dict):
        binding_site = {}

    lines = [
        f"# Evidence Package: {target.get('display_name', target.get('id', 'unknown'))}",
        "",
        f"- Target ID: `{target.get('id', '')}`",
        f"- Recommendation: `{target.get('recommendation', '')}`",
        f"- Reason: {target.get('recommendation_reason', '')}",
        "",
        "## Evidence Readiness",
        "",
        f"- Evidence score: `{quality.get('score', '')}`",
        f"- Readiness: `{quality.get('readiness', '')}`",
        "- Meaning: source readiness for computational closed-loop screening, not measured efficacy.",
        "",
        "| component | score |",
        "|---|---:|",
    ]
    for key in [
        "clinical_controls",
        "structure_metadata",
        "literature_metadata",
        "assay_path",
        "open_data_entry_points",
        "catalog_prior",
    ]:
        lines.append(f"| {key} | {components.get(key, '')} |")

    lines.extend(
        [
            "",
            "## Closed-Loop Assets",
            "",
            f"- Primary structure: `{primary_structure.get('pdb_id', '')}` / ligand `{primary_structure.get('ligand', '')}` / use `{primary_structure.get('use', '')}`",
            f"- Structure metadata: `{primary_structure.get('metadata_status', '')}`, resolution `{primary_structure.get('resolution', '')}`",
            f"- Binding-site source: `{binding_site.get('source', '')}`",
            f"- Reference ligand: `{binding_site.get('reference_ligand', '')}`",
            f"- Docking box strategy: {binding_site.get('box_strategy', '')}",
            f"- Positive controls: {', '.join(str(item) for item in assets.get('positive_controls', []) if item) or 'none'}",
            f"- Reference controls: {', '.join(str(item) for item in assets.get('reference_controls', []) if item) or 'none'}",
            f"- Historical controls: {', '.join(str(item) for item in assets.get('historical_controls', []) if item) or 'none'}",
            f"- Assay path: {assets.get('assay_path', '')}",
            "",
        "## Known Drugs / Controls",
        "",
        ]
    )
    known_drug_rows = package.get("known_drug_rows", [])
    if isinstance(known_drug_rows, list) and known_drug_rows:
        lines.extend(["| drug | mechanism | workflow role | source |", "|---|---|---|---|"])
        for row in known_drug_rows:
            if not isinstance(row, dict):
                continue
            lines.append(
                f"| {row.get('drug', '')} | {row.get('mechanism', '')} | "
                f"{row.get('status_for_workflow', '')} | [source]({row.get('source_url', '')}) |"
            )
    else:
        known_drugs = target.get("known_drugs", [])
        if isinstance(known_drugs, list) and known_drugs:
            lines.extend([f"- {drug}" for drug in known_drugs])
        else:
            lines.append("- Not specified")

    lines.extend(["", "## Structures", "", "| PDB | status | title | method | resolution | PubMed |", "|---|---|---|---|---:|---|"])
    for entry in package.get("pdb_entries", []):
        if not isinstance(entry, dict):
            continue
        lines.append(
            f"| [{entry.get('pdb_id', '')}]({entry.get('rcsb_url', '')}) | {entry.get('fetch_status', '')} | "
            f"{entry.get('title', '')} | {entry.get('experimental_method', '')} | {entry.get('resolution', '')} | "
            f"{entry.get('primary_citation_pubmed_id', '')} |"
        )

    lines.extend(["", "## PubMed Query Results", "", "| PMID | year/date | first author | title |", "|---|---|---|---|"])
    articles = package.get("pubmed_articles", [])
    if isinstance(articles, list) and articles:
        for article in articles[:20]:
            if not isinstance(article, dict):
                continue
            lines.append(
                f"| [{article.get('pmid', '')}]({article.get('pubmed_url', '')}) | "
                f"{article.get('pubdate', '')} | {article.get('first_author', '')} | {article.get('title', '')} |"
            )
    else:
        lines.append("| - | - | - | No PubMed metadata fetched. Use the queued queries or rerun with network access. |")

    lines.extend(["", "## Open Data Entry Points", ""])
    for query in package.get("open_data_queries", []):
        lines.append(f"- {query}")

    lines.extend(["", "## Official / Primary Sources", ""])
    for source in package.get("official_sources", []):
        if not isinstance(source, dict):
            continue
        lines.append(f"- [{source.get('title', source.get('id', 'source'))}]({source.get('url', '')}) - {source.get('evidence_role', '')}")

    errors = package.get("network_errors", [])
    if isinstance(errors, list) and errors:
        lines.extend(["", "## Network / Fetch Warnings", ""])
        for error in errors:
            lines.append(f"- `{error.get('source', 'unknown')}`: {error.get('error', '')}")

    lines.extend(
        [
            "",
            "## Evidence Counts",
            "",
            f"- PDB entries: {counts.get('pdb_entries', 0)}",
            f"- PubMed articles: {counts.get('pubmed_articles', 0)}",
            f"- Positive controls: {counts.get('positive_controls', 0)}",
            f"- Open data queries: {counts.get('open_data_queries', 0)}",
            "",
            "## Interpretation Boundary",
            "",
            "- This is a target-evidence package for computational screening.",
            "- It does not establish therapeutic efficacy.",
            "- Full-text papers are not redistributed; only metadata, links, and local summaries are stored.",
        ]
    )
    return "\n".join(lines) + "\n"


def write_evidence_package(
    output_root: Path,
    target: Dict[str, object],
    package: Dict[str, object],
) -> None:
    target_id = str(target.get("id", "unknown"))
    target_dir = output_root / target_id
    target_dir.mkdir(parents=True, exist_ok=True)
    write_json(target_dir / "evidence.json", package)

    pdb_entries = [entry for entry in package.get("pdb_entries", []) if isinstance(entry, dict)]
    if pdb_entries:
        write_csv(
            target_dir / "pdb_entries.csv",
            pdb_entries,
            [
                "pdb_id",
                "fetch_status",
                "title",
                "experimental_method",
                "resolution",
                "release_date",
                "primary_citation_title",
                "primary_citation_year",
                "primary_citation_pubmed_id",
                "rcsb_url",
            ],
        )

    articles = [article for article in package.get("pubmed_articles", []) if isinstance(article, dict)]
    write_csv(
        target_dir / "pubmed_articles.csv",
        articles,
        ["pmid", "title", "journal", "pubdate", "first_author", "doi_or_elocation", "pubmed_url", "query"],
    )
    write_text(target_dir / "evidence_report.md", render_evidence_report(target, package))


def evidence_refresh(args: argparse.Namespace) -> None:
    catalog = load_target_catalog(args.catalog)
    sources = load_evidence_sources(args.sources)
    by_id = catalog_targets_by_id(catalog)
    target_ids = target_ids_for_evidence(args, catalog)
    output_root = Path(args.output).expanduser().resolve() if args.output else DEFAULT_EVIDENCE_DIR
    output_root.mkdir(parents=True, exist_ok=True)
    write_json(output_root / "evidence_sources_snapshot.json", sources)

    source_queries = sources.get("target_queries", {})
    if not isinstance(source_queries, dict):
        source_queries = {}
    official_sources = sources.get("official_sources", [])
    open_sources = sources.get("open_data_sources", [])
    catalog_weights = catalog.get("scoring_weights", {})
    if not isinstance(catalog_weights, dict):
        catalog_weights = {}
    known_drug_rows_all = read_csv(DEFAULT_KNOWN_DRUGS)
    structure_rows_all = read_csv(DEFAULT_PDB_STRUCTURES)

    summary_rows = []
    for target_id in target_ids:
        target = by_id[target_id]
        existing_package = load_existing_evidence_package(output_root, target_id)
        structures = target.get("representative_structures", [])
        pdb_ids = []
        if isinstance(structures, list):
            pdb_ids = [str(item.get("pdb_id")) for item in structures if isinstance(item, dict) and item.get("pdb_id")]

        network_errors: List[Dict[str, str]] = []
        pdb_entries = []
        existing_pdb = existing_items_by_key(existing_package, "pdb_entries", "pdb_id")
        if args.offline:
            for pdb_id in pdb_ids:
                cached = existing_pdb.get(pdb_id)
                if cached:
                    pdb_entries.append(cached)
                else:
                    pdb_entries.append(
                        {
                            "pdb_id": pdb_id,
                            "fetch_status": "queued_offline",
                            "rcsb_url": f"https://www.rcsb.org/structure/{pdb_id}",
                        }
                    )
        else:
            for pdb_id in pdb_ids:
                entry, error = fetch_rcsb_metadata(pdb_id, args.timeout)
                if error and pdb_id in existing_pdb:
                    entry = dict(existing_pdb[pdb_id])
                    entry["fetch_status"] = entry.get("fetch_status") or "cached_after_fetch_error"
                pdb_entries.append(entry)
                if error:
                    network_errors.append({"source": "rcsb", "query": pdb_id, "error": error})
                time.sleep(0.20)

        target_query_block = source_queries.get(target_id, {})
        if not isinstance(target_query_block, dict):
            target_query_block = {}
        pubmed_queries = target_query_block.get("pubmed", [])
        if not isinstance(pubmed_queries, list):
            pubmed_queries = []
        open_data_queries = target_query_block.get("open_data_queries", [])
        if not isinstance(open_data_queries, list):
            open_data_queries = []

        if args.offline:
            existing_articles = existing_package.get("pubmed_articles", [])
            pubmed_articles = [item for item in existing_articles if isinstance(item, dict)] if isinstance(existing_articles, list) else []
            pubmed_errors = []
        else:
            pubmed_articles, pubmed_errors = collect_pubmed_articles(pubmed_queries, args.retmax, args.timeout, args.offline)
            if not pubmed_articles:
                existing_articles = existing_package.get("pubmed_articles", [])
                if isinstance(existing_articles, list) and existing_articles:
                    pubmed_articles = [item for item in existing_articles if isinstance(item, dict)]
        network_errors.extend(pubmed_errors)

        package = {
            "schema_version": "0.1",
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "target_id": target_id,
            "target_display_name": target.get("display_name", ""),
            "recommendation": target.get("recommendation", ""),
            "known_drugs": target.get("known_drugs", []),
            "official_sources": official_sources,
            "open_data_sources": open_sources,
            "pdb_entries": pdb_entries,
            "pubmed_queries": pubmed_queries,
            "pubmed_articles": pubmed_articles,
            "open_data_queries": open_data_queries,
            "network_errors": network_errors,
        }
        known_drug_rows = rows_for_target(known_drug_rows_all, target_id)
        structure_rows = rows_for_target(structure_rows_all, target_id)
        package = enrich_evidence_package(target, catalog_weights, package, known_drug_rows, structure_rows)
        write_evidence_package(output_root, target, package)
        quality = package.get("evidence_quality", {})
        assets = package.get("closed_loop_assets", {})
        if not isinstance(quality, dict):
            quality = {}
        if not isinstance(assets, dict):
            assets = {}
        primary_structure = assets.get("primary_structure", {})
        if not isinstance(primary_structure, dict):
            primary_structure = {}
        summary_rows.append(
            {
                "target_id": target_id,
                "target": target.get("display_name", ""),
                "pdb_entries": len(pdb_entries),
                "pubmed_articles": len(pubmed_articles),
                "positive_controls": len(assets.get("positive_controls", [])) if isinstance(assets.get("positive_controls", []), list) else 0,
                "primary_pdb": primary_structure.get("pdb_id", ""),
                "evidence_score": quality.get("score", ""),
                "readiness": quality.get("readiness", ""),
                "network_errors": len(network_errors),
                "evidence_dir": str(output_root / target_id),
            }
        )
        print(
            f"{target_id}: pdb={len(pdb_entries)} pubmed={len(pubmed_articles)} "
            f"score={quality.get('score', '')} readiness={quality.get('readiness', '')} "
            f"errors={len(network_errors)}"
        )

    write_csv(
        output_root / "evidence_summary.csv",
        summary_rows,
        [
            "target_id",
            "target",
            "pdb_entries",
            "pubmed_articles",
            "positive_controls",
            "primary_pdb",
            "evidence_score",
            "readiness",
            "network_errors",
            "evidence_dir",
        ],
    )
    write_text(output_root / "README.md", render_evidence_index(summary_rows, sources))
    print(f"Wrote evidence packages: {output_root}")


def render_evidence_index(summary_rows: Sequence[Dict[str, object]], sources: Dict[str, object]) -> str:
    lines = [
        "# Influenza Target Evidence Index",
        "",
        "This directory stores target-evidence metadata used by the second-stage target-source workflow.",
        "",
        "## Target Packages",
        "",
        "| target id | target | evidence score | readiness | primary PDB | PDB entries | PubMed articles | controls | fetch warnings |",
        "|---|---|---:|---|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| `{row.get('target_id', '')}` | {row.get('target', '')} | "
            f"{row.get('evidence_score', '')} | {row.get('readiness', '')} | {row.get('primary_pdb', '')} | "
            f"{row.get('pdb_entries', 0)} | {row.get('pubmed_articles', 0)} | "
            f"{row.get('positive_controls', 0)} | {row.get('network_errors', 0)} |"
        )
    lines.extend(["", "## Source Types", ""])
    for section in ["official_sources", "open_data_sources"]:
        items = sources.get(section, [])
        lines.append(f"### {section}")
        lines.append("")
        if isinstance(items, list):
            for item in items:
                if isinstance(item, dict):
                    lines.append(f"- [{item.get('title', item.get('id', 'source'))}]({item.get('url', '')})")
        lines.append("")
    lines.extend(
        [
            "## Boundary",
            "",
            "The evidence packages store metadata, source links, and local summaries. They do not redistribute copyrighted full text.",
        ]
    )
    return "\n".join(lines) + "\n"


def csv_join(items: object) -> str:
    if isinstance(items, list):
        return "; ".join(str(item) for item in items if item)
    if items:
        return str(items)
    return ""


def stage2_target_row(
    target: Dict[str, object],
    package: Dict[str, object],
    catalog_weights: Dict[str, object],
    evidence_dir: Path,
) -> Dict[str, object]:
    quality = package.get("evidence_quality", {})
    if not isinstance(quality, dict):
        quality = {}
    counts = quality.get("counts", {})
    if not isinstance(counts, dict):
        counts = {}
    assets = package.get("closed_loop_assets", {})
    if not isinstance(assets, dict):
        assets = {}
    primary_structure = assets.get("primary_structure", {})
    if not isinstance(primary_structure, dict):
        primary_structure = {}
    binding_site = assets.get("binding_site", {})
    if not isinstance(binding_site, dict):
        binding_site = {}

    target_id = str(target.get("id", ""))
    return {
        "target_id": target_id,
        "target": target.get("display_name", ""),
        "catalog_score": target_total_score(target, catalog_weights),
        "evidence_score": quality.get("score", ""),
        "readiness": quality.get("readiness", ""),
        "recommendation": target.get("recommendation", ""),
        "primary_pdb": primary_structure.get("pdb_id", ""),
        "reference_ligand": binding_site.get("reference_ligand", ""),
        "positive_controls": csv_join(assets.get("positive_controls", [])),
        "reference_controls": csv_join(assets.get("reference_controls", [])),
        "historical_controls": csv_join(assets.get("historical_controls", [])),
        "pdb_entries": counts.get("pdb_entries", len(package.get("pdb_entries", [])) if isinstance(package.get("pdb_entries", []), list) else 0),
        "pubmed_articles": counts.get("pubmed_articles", len(package.get("pubmed_articles", [])) if isinstance(package.get("pubmed_articles", []), list) else 0),
        "open_data_queries": counts.get("open_data_queries", len(package.get("open_data_queries", [])) if isinstance(package.get("open_data_queries", []), list) else 0),
        "assay_path": assets.get("assay_path", ""),
        "binding_site_source": binding_site.get("source", ""),
        "box_strategy": binding_site.get("box_strategy", ""),
        "evidence_report": str(evidence_dir / target_id / "evidence_report.md"),
    }


def render_stage2_project_report(
    disease: str,
    rows: Sequence[Dict[str, object]],
    sources: Dict[str, object],
    package_root: Path,
) -> str:
    top = rows[0] if rows else {}
    lines = [
        "# Stage 2 Target Evidence Source Report",
        "",
        f"- Disease/query: `{disease}`",
        f"- Evidence package root: `{package_root}`",
        "- Scope: target-source selection for computational screening; not a therapeutic efficacy claim.",
        "",
        "## Decision",
        "",
    ]
    if top:
        lines.extend(
            [
                f"- First target to enter the closed loop: `{top.get('target_id', '')}` ({top.get('target', '')})",
                f"- Evidence score: `{top.get('evidence_score', '')}`",
                f"- Readiness: `{top.get('readiness', '')}`",
                f"- Primary structure: `{top.get('primary_pdb', '')}`",
                f"- Positive controls: {top.get('positive_controls', '')}",
                f"- Assay path: {top.get('assay_path', '')}",
            ]
        )
    else:
        lines.append("- No target packages found.")

    lines.extend(
        [
            "",
            "## Target Source Matrix",
            "",
            "| rank | target id | target | catalog | evidence | readiness | PDB | PubMed | controls | primary PDB |",
            "|---:|---|---|---:|---:|---|---:|---:|---|---|",
        ]
    )
    for idx, row in enumerate(rows, start=1):
        lines.append(
            f"| {idx} | `{row.get('target_id', '')}` | {row.get('target', '')} | "
            f"{row.get('catalog_score', '')} | {row.get('evidence_score', '')} | {row.get('readiness', '')} | "
            f"{row.get('pdb_entries', '')} | {row.get('pubmed_articles', '')} | "
            f"{row.get('positive_controls', '')} | {row.get('primary_pdb', '')} |"
        )

    lines.extend(
        [
            "",
            "## What This Adds To The Closed Loop",
            "",
            "1. 靶点不再靠口头指定，而是由疾病、公开结构、已验证药物、文献元数据和可验证 assay 路径共同筛选。",
            "2. 生成模型的输入可以绑定到明确靶点、参考配体、阳性对照和 docking box 策略。",
            "3. 后续每轮候选分子的打分可以与已知药物和公开结构对照，形成可解释反馈。",
            "4. 当前证据包保存的是元数据和链接，不保存受版权保护的论文全文。",
            "",
            "## Open Data Sources",
            "",
        ]
    )
    for section in ["official_sources", "open_data_sources"]:
        items = sources.get(section, [])
        if not isinstance(items, list):
            continue
        lines.append(f"### {section}")
        lines.append("")
        for item in items:
            if isinstance(item, dict):
                lines.append(
                    f"- [{item.get('title', item.get('id', 'source'))}]({item.get('url', '')}) - "
                    f"{item.get('evidence_role', '')}"
                )
        lines.append("")

    lines.extend(
        [
            "## Boundary",
            "",
            "- Evidence score is a source-readiness score, not activity, potency, safety, or clinical efficacy.",
            "- Docking/proxy results must not be presented as measured drug effect.",
            "- Full validation still requires prepared receptor structures, real docking/pose checks, and preferably biochemical assays.",
        ]
    )
    return "\n".join(lines) + "\n"


def evidence_stage2(args: argparse.Namespace) -> None:
    project = project_path(args.project)
    ensure_project_dirs(project)
    evidence_root = Path(args.evidence_dir).expanduser().resolve() if args.evidence_dir else DEFAULT_EVIDENCE_DIR

    if args.refresh:
        evidence_refresh(
            argparse.Namespace(
                disease=args.disease,
                target=args.target,
                catalog=args.catalog,
                sources=args.sources,
                output=str(evidence_root),
                retmax=args.retmax,
                timeout=args.timeout,
                offline=args.offline,
            )
        )

    catalog = load_target_catalog(args.catalog)
    sources = load_evidence_sources(args.sources)
    by_id = catalog_targets_by_id(catalog)
    target_ids = target_ids_for_evidence(args, catalog)
    if args.top and args.top > 0 and not args.target:
        target_ids = target_ids[: args.top]
    catalog_weights = catalog.get("scoring_weights", {})
    if not isinstance(catalog_weights, dict):
        catalog_weights = {}
    known_drug_rows_all = read_csv(DEFAULT_KNOWN_DRUGS)
    structure_rows_all = read_csv(DEFAULT_PDB_STRUCTURES)

    rows: List[Dict[str, object]] = []
    packages: List[Dict[str, object]] = []
    for target_id in target_ids:
        target = by_id.get(target_id)
        if not target:
            continue
        package = load_existing_evidence_package(evidence_root, target_id)
        if not package:
            package = {
                "schema_version": "0.1",
                "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
                "target_id": target_id,
                "target_display_name": target.get("display_name", ""),
                "recommendation": target.get("recommendation", ""),
                "known_drugs": target.get("known_drugs", []),
                "official_sources": sources.get("official_sources", []),
                "open_data_sources": sources.get("open_data_sources", []),
                "pdb_entries": [],
                "pubmed_queries": [],
                "pubmed_articles": [],
                "open_data_queries": [],
                "network_errors": [{"source": "local_evidence", "error": "Evidence package not found; run evidence-refresh."}],
            }
        package = enrich_evidence_package(
            target,
            catalog_weights,
            package,
            rows_for_target(known_drug_rows_all, target_id),
            rows_for_target(structure_rows_all, target_id),
        )
        write_evidence_package(evidence_root, target, package)
        rows.append(stage2_target_row(target, package, catalog_weights, evidence_root))
        packages.append(
            {
                "target_id": target_id,
                "target": target.get("display_name", ""),
                "evidence_quality": package.get("evidence_quality", {}),
                "closed_loop_assets": package.get("closed_loop_assets", {}),
                "pubmed_article_count": len(package.get("pubmed_articles", [])) if isinstance(package.get("pubmed_articles", []), list) else 0,
                "pdb_entry_count": len(package.get("pdb_entries", [])) if isinstance(package.get("pdb_entries", []), list) else 0,
                "open_data_queries": package.get("open_data_queries", []),
                "evidence_report": str(evidence_root / target_id / "evidence_report.md"),
            }
        )

    rows.sort(
        key=lambda row: (
            float(row.get("evidence_score") or 0.0),
            float(row.get("catalog_score") or 0.0),
        ),
        reverse=True,
    )
    evidence_out = project / "evidence"
    reports_out = project / "reports"
    matrix_path = evidence_out / "stage2_target_sources.csv"
    assets_path = evidence_out / "stage2_closed_loop_assets.json"
    report_path = reports_out / "stage2_evidence_report.md"
    write_csv(
        matrix_path,
        rows,
        [
            "target_id",
            "target",
            "catalog_score",
            "evidence_score",
            "readiness",
            "recommendation",
            "primary_pdb",
            "reference_ligand",
            "positive_controls",
            "reference_controls",
            "historical_controls",
            "pdb_entries",
            "pubmed_articles",
            "open_data_queries",
            "assay_path",
            "binding_site_source",
            "box_strategy",
            "evidence_report",
        ],
    )
    write_json(
        assets_path,
        {
            "schema_version": "0.1",
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "disease": args.disease,
            "evidence_root": str(evidence_root),
            "targets": packages,
        },
    )
    write_text(report_path, render_stage2_project_report(args.disease, rows, sources, evidence_root))

    print(f"Wrote stage 2 target-source matrix: {matrix_path}")
    print(f"Wrote stage 2 closed-loop assets: {assets_path}")
    print(f"Wrote stage 2 report: {report_path}")


def split_items(value: Optional[str]) -> List[str]:
    if not value:
        return []
    parts = re.split(r"[;；\n]+", value)
    if len(parts) == 1:
        parts = re.split(r"[,，]+", value)
    return [part.strip() for part in parts if part.strip()]


def parse_triplet(value: Optional[str]) -> List[float]:
    if not value:
        return []
    parts = [part.strip() for part in re.split(r"[,，\s]+", value) if part.strip()]
    if len(parts) != 3:
        raise SystemExit(f"Expected three numbers, got: {value}")
    try:
        return [float(part) for part in parts]
    except ValueError as exc:
        raise SystemExit(f"Invalid numeric triplet: {value}") from exc


def target_brief_from_args(args: argparse.Namespace) -> Dict[str, object]:
    return {
        "schema_version": "0.1",
        "target_name": args.target_name,
        "disease_context": args.disease,
        "target_rationale": args.rationale,
        "free_text_requirement": args.free_text,
        "protein": {
            "name": args.protein,
            "gene": args.gene,
            "pdb_id": args.pdb_id,
            "structure_file": args.protein_pdb,
        },
        "binding_site": {
            "description": args.pocket,
            "reference_ligand": args.reference_ligand,
            "key_residues": split_items(args.key_residues),
            "center": parse_triplet(args.center),
            "size": parse_triplet(args.size),
            "source": args.pocket_source,
        },
        "design_intent": {
            "style": args.style,
            "primary_goal": args.primary_goal,
            "desired_activity": args.activity,
            "must_have": split_items(args.must_have),
            "avoid": split_items(args.avoid),
            "desired_properties": split_items(args.desired_properties),
            "selectivity_notes": args.selectivity,
        },
        "generation_constraints": {
            "max_heavy_atoms": args.max_heavy_atoms,
            "max_molecular_weight": args.max_molecular_weight,
            "prefer_synthesizable": True,
            "preserve_reference_scaffold": args.preserve_scaffold,
            "diversity_required": args.diversity,
        },
        "validation_plan": {
            "first_pass": "proxy filters, validity, QED, synthetic accessibility",
            "structure_scoring": "AutoDock Vina or GNINA docking score",
            "pose_quality": "PoseBusters pass/fail or pose quality score",
            "high_cost_followup": "OpenFE or higher-cost molecular simulation for top hits only",
        },
        "compliance_boundary": {
            "scope": "computational screening and research planning only",
            "no_claims": [
                "Do not claim therapeutic efficacy.",
                "Do not provide dosing, clinical advice, or experimental human-use claims.",
                "Do not present proxy scores as measured potency.",
                "Do not include synthesis procedures unless separately requested for benign feasibility review.",
            ],
        },
    }


def render_generator_prompt(brief: Dict[str, object]) -> str:
    protein = brief.get("protein", {})
    site = brief.get("binding_site", {})
    intent = brief.get("design_intent", {})
    constraints = brief.get("generation_constraints", {})
    compliance = brief.get("compliance_boundary", {})

    if not isinstance(protein, dict):
        protein = {}
    if not isinstance(site, dict):
        site = {}
    if not isinstance(intent, dict):
        intent = {}
    if not isinstance(constraints, dict):
        constraints = {}
    if not isinstance(compliance, dict):
        compliance = {}

    def list_block(items: object) -> str:
        if not isinstance(items, list) or not items:
            return "- Not specified"
        return "\n".join(f"- {item}" for item in items)

    return f"""# Generator Prompt

You are assisting with computational small-molecule design for virtual screening.
Generate hypothetical candidate molecules only for in silico evaluation. Do not
claim clinical efficacy, measured potency, safety, or human-use relevance.

## Target Brief

- Target name: {brief.get("target_name") or "Not specified"}
- Disease context: {brief.get("disease_context") or "Not specified"}
- Target rationale: {brief.get("target_rationale") or "Not specified"}
- Protein: {protein.get("name") or "Not specified"}
- Gene: {protein.get("gene") or "Not specified"}
- PDB ID: {protein.get("pdb_id") or "Not specified"}
- Structure file: {protein.get("structure_file") or "Not specified"}

## Binding-Site Hypothesis

- Pocket source: {site.get("source") or "Not specified"}
- Pocket description: {site.get("description") or "Not specified"}
- Reference ligand: {site.get("reference_ligand") or "Not specified"}
- Key residues:
{list_block(site.get("key_residues"))}
- Docking center: {site.get("center") or "Not specified"}
- Docking box size: {site.get("size") or "Not specified"}

## Molecule Design Intent

- Design style: {intent.get("style") or "Not specified"}
- Primary goal: {intent.get("primary_goal") or "Not specified"}
- Desired activity: {intent.get("desired_activity") or "Not specified"}
- Selectivity notes: {intent.get("selectivity_notes") or "Not specified"}

Must-have features:
{list_block(intent.get("must_have"))}

Avoid:
{list_block(intent.get("avoid"))}

Desired properties:
{list_block(intent.get("desired_properties"))}

## Hard Constraints

- Max heavy atoms: {constraints.get("max_heavy_atoms")}
- Max molecular weight: {constraints.get("max_molecular_weight")}
- Prefer synthesizable molecules: {constraints.get("prefer_synthesizable")}
- Preserve reference scaffold: {constraints.get("preserve_reference_scaffold")}
- Maintain chemical diversity: {constraints.get("diversity_required")}

## Output Contract

Return candidates as rows with these columns:

```csv
id,smiles,rationale,expected_interaction,design_family,risk_note
```

Requirements:

- Produce valid small-molecule SMILES.
- Keep molecules drug-like and suitable for docking.
- Propose diverse analog families rather than near duplicates.
- Explain the intended pocket interaction at a high level.
- Flag obvious risks such as reactive groups, excessive size, or poor polarity.
- Do not provide synthetic procedures or biological use instructions.

## Compliance Boundary

- Scope: {compliance.get("scope") or "computational screening only"}
{list_block(compliance.get("no_claims"))}

## Free-Text User Requirement

{brief.get("free_text_requirement") or "Not specified"}
"""


def brief_command(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    ensure_project_dirs(project)
    config_path = project / "config.json"
    config = load_config(project) if config_path.exists() else default_config()
    free_text = args.free_text
    if args.from_file:
        free_text = Path(args.from_file).expanduser().resolve().read_text(encoding="utf-8")
    args.free_text = free_text

    brief = target_brief_from_args(args)
    brief_path = project / "briefs" / "target_brief.json"
    prompt_path = project / "prompts" / "generator_prompt.md"
    if (brief_path.exists() or prompt_path.exists()) and not args.force:
        raise SystemExit("Target brief or prompt already exists. Use --force to overwrite.")

    write_json(brief_path, brief)
    write_text(prompt_path, render_generator_prompt(brief))

    target = config.get("target", {})
    if not isinstance(target, dict):
        target = {}
    target.update(
        {
            "name": args.target_name,
            "protein_pdb": args.protein_pdb,
            "reference_ligand_sdf": args.reference_ligand,
            "pocket": {
                "center": parse_triplet(args.center),
                "size": parse_triplet(args.size),
                "source": args.pocket_source,
                "description": args.pocket,
                "key_residues": split_items(args.key_residues),
            },
        }
    )
    config["target"] = target
    config["target_brief_file"] = str(brief_path)
    config["generator_prompt_file"] = str(prompt_path)
    write_json(project / "config.json", config)

    print(f"Wrote target brief: {brief_path}")
    print(f"Wrote generator prompt: {prompt_path}")
    return prompt_path


PROMPT_EXAMPLES = [
    (
        "共晶配体口袋型",
        """为一个已有蛋白-配体共晶结构设计候选小分子。靶点是 SARS-CoV-2 Mpro，设计口袋以共晶配体 N3 所在位置为中心，关键残基包括 His41、Cys145、Gly143、His164、Glu166。希望生成可进入催化口袋、能形成氢键或极性相互作用、同时保持较低反应性和中等分子量的候选分子。输出仅用于虚拟筛选，不声称真实活性。""",
    ),
    (
        "P2Rank 预测口袋型",
        """已知某疾病相关蛋白缺少可靠共晶配体，但 P2Rank 预测出一个体积适中、疏水沟槽和极性边缘并存的潜在口袋。请生成一组适合该口袋的药物样小分子：核心结构可以占据疏水区域，边缘保留 1-3 个氢键供受体用于与极性残基作用，避免过大、过柔性和明显反应性结构。所有候选只进入 docking 和 pose 检查。""",
    ),
    (
        "选择性优化型",
        """围绕 EGFR kinase ATP 结合口袋设计候选分子，但需要尽量降低对同源激酶的非选择性风险。候选应保留一个可与 hinge 区形成氢键的杂芳香核心，同时在 solvent-exposed 区域引入可调极性侧链。避免过强疏水堆叠导致广谱 kinase off-target。输出候选 SMILES、设计理由和可能选择性风险。""",
    ),
    (
        "片段生长型",
        """从一个已知弱结合片段出发做 fragment growing。片段锚定在口袋深部，附近存在一个可扩展的疏水子口袋和一个靠近溶剂的极性出口。请保留片段核心，向疏水子口袋延伸小体积芳香或脂肪基团，同时在溶剂出口方向加入提高溶解性的极性基团。候选应保持合成可及性和结构多样性。""",
    ),
    (
        "ADMET 平衡型",
        """针对一个已知靶点生成候选抑制剂，但重点不是最大化 docking 分数，而是平衡结合力、类药性和早期 ADMET 风险。希望分子量低于 500，氢键供受体数量适中，避免 PAINS、强反应性官能团、过高疏水性和过多芳环。输出适合进入 Vina/GNINA docking 的多样化候选。""",
    ),
    (
        "耐药突变规避型",
        """为一个存在耐药突变的结合口袋设计候选分子。设计时不要过度依赖单一突变热点残基的强相互作用，而应分散相互作用网络：保留与保守残基的氢键，增加与主链或稳定疏水区域的接触，降低因单点突变导致失效的风险。输出候选 SMILES 和每个分子的相互作用假设。""",
    ),
]


def prompt_examples(args: argparse.Namespace) -> None:
    sections = ["# Target Requirement Prompt Examples", ""]
    for idx, (title, text) in enumerate(PROMPT_EXAMPLES, start=1):
        sections.extend([f"## {idx}. {title}", "", text.strip(), ""])
    content = "\n".join(sections)
    if args.project and args.write:
        project = project_path(args.project)
        ensure_project_dirs(project)
        output = project / "prompts" / "target_prompt_examples.md"
        write_text(output, content)
        print(f"Wrote prompt examples: {output}")
    print(content)


def load_config(project: Path) -> Dict[str, object]:
    config_path = project / "config.json"
    if not config_path.exists():
        raise SystemExit(f"Missing config: {config_path}. Run init first.")
    return read_json(config_path)


def seed_file_for_generation(project: Path, round_no: int) -> Path:
    preferred = project / "seeds" / f"round_{round_no - 1}_seeds.csv"
    if preferred.exists():
        return preferred
    return project / "seeds" / "round_0_seeds.csv"


def load_seed_smiles(project: Path, round_no: int) -> List[str]:
    rows = read_csv(seed_file_for_generation(project, round_no))
    smiles = [row.get("smiles", "").strip() for row in rows if row.get("smiles", "").strip()]
    return smiles or [item[1] for item in DEFAULT_SEEDS]


def aromatic_variant(parent: str, substituent: str) -> str:
    if "c1ccccc1" in parent:
        return parent + substituent
    if "c1ccccc1" in substituent:
        return substituent
    if len(parent) <= 10:
        return "c1ccccc1" + substituent
    return parent + substituent


def chain_variant(parent: str, substituent: str) -> str:
    if parent.endswith(")"):
        return substituent + parent
    return parent + substituent


def generate_candidates(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    load_config(project)
    ensure_project_dirs(project)

    if args.source_csv:
        source_rows = read_csv(Path(args.source_csv).expanduser().resolve())
        candidates = []
        for idx, row in enumerate(source_rows, start=1):
            smiles = row.get("smiles", "").strip()
            if not smiles:
                continue
            candidates.append(
                {
                    "round": args.round,
                    "id": row.get("id") or f"r{args.round:02d}_import_{idx:05d}",
                    "smiles": smiles,
                    "parent": row.get("parent", ""),
                    "source": "external_csv",
                }
            )
    else:
        seed_smiles = load_seed_smiles(project, args.round)
        rng = random.Random(args.round * 7919 + args.n)
        pool: List[Tuple[str, str, str]] = []
        for parent in seed_smiles:
            for substituent in SUBSTITUENTS:
                pool.append((aromatic_variant(parent, substituent), parent, "aromatic_or_parent_substitution"))
                pool.append((chain_variant(parent, substituent), parent, "chain_extension"))
        for molecule in TEMPLATE_MOLECULES:
            parent = rng.choice(seed_smiles)
            pool.append((molecule, parent, "template_library"))
        rng.shuffle(pool)

        seen = set()
        candidates = []
        for smiles, parent, source in pool:
            if smiles in seen:
                continue
            seen.add(smiles)
            idx = len(candidates) + 1
            candidates.append(
                {
                    "round": args.round,
                    "id": f"r{args.round:02d}_{idx:05d}",
                    "smiles": smiles,
                    "parent": parent,
                    "source": source,
                }
            )
            if len(candidates) >= args.n:
                break

    output = project / "candidates" / f"round_{args.round}_candidates.csv"
    write_csv(output, candidates, ["round", "id", "smiles", "parent", "source"])
    print(f"Wrote {len(candidates)} candidates: {output}")
    return output


def balanced(text: str, left: str, right: str) -> bool:
    depth = 0
    for char in text:
        if char == left:
            depth += 1
        elif char == right:
            depth -= 1
            if depth < 0:
                return False
    return depth == 0


def ring_digits_ok(smiles: str) -> bool:
    counts = Counter(re.findall(r"\d", smiles))
    return all(count % 2 == 0 for count in counts.values())


def allowed_chars_ok(smiles: str) -> bool:
    return bool(re.fullmatch(r"[A-Za-z0-9@+\-\[\]\(\)=#\\/%.]+", smiles))


def plausible_smiles(smiles: str) -> bool:
    if not smiles or len(smiles) > 180:
        return False
    if not allowed_chars_ok(smiles):
        return False
    if not balanced(smiles, "(", ")") or not balanced(smiles, "[", "]"):
        return False
    if not ring_digits_ok(smiles):
        return False
    if smiles.count("=") + smiles.count("#") > max(1, len(smiles) // 3):
        return False
    return True


def count_atoms(smiles: str) -> Counter:
    tokens = re.findall(r"Cl|Br|[BCNOFPSI]|[cnops]", smiles)
    counts: Counter = Counter()
    for token in tokens:
        key = token.capitalize() if token.islower() else token
        counts[key] += 1
    return counts


def target_score(value: float, target: float, tolerance: float) -> float:
    if tolerance <= 0:
        return 0.0
    return clamp(1.0 - abs(value - target) / tolerance)


def range_score(value: float, low: float, high: float) -> float:
    if low <= value <= high:
        return 1.0
    if value < low:
        return clamp(value / low if low else 0.0)
    span = max(high - low, 1.0)
    return clamp(1.0 - (value - high) / span)


def bigrams(text: str) -> set:
    if len(text) < 2:
        return {text}
    return {text[i : i + 2] for i in range(len(text) - 1)}


def tanimoto_string(a: str, b: str) -> float:
    x = bigrams(a)
    y = bigrams(b)
    if not x and not y:
        return 1.0
    return len(x & y) / max(len(x | y), 1)


def novelty_proxy(smiles: str, seeds: Sequence[str]) -> float:
    if not seeds:
        return 1.0
    max_sim = max(tanimoto_string(smiles, seed) for seed in seeds)
    return clamp(1.0 - max_sim)


def estimate_scores(smiles: str, seeds: Sequence[str]) -> Dict[str, float]:
    valid = 1.0 if plausible_smiles(smiles) else 0.0
    counts = count_atoms(smiles)
    heavy = sum(counts.values())
    hetero = counts["N"] + counts["O"] + counts["S"] + counts["P"]
    halogens = counts["F"] + counts["Cl"] + counts["Br"] + counts["I"]
    aromatic = smiles.count("c") + smiles.count("n") + smiles.count("o") + smiles.count("s")
    branches = smiles.count("(")
    rings = len(re.findall(r"\d", smiles)) / 2.0
    length = len(smiles)

    lipinski = (
        range_score(heavy, 8, 55) * 0.45
        + range_score(hetero, 1, 14) * 0.25
        + range_score(halogens, 0, 6) * 0.10
        + range_score(length, 10, 95) * 0.20
    )
    qed = (
        target_score(heavy, 28, 24) * 0.35
        + target_score(hetero, 5, 8) * 0.25
        + target_score(aromatic, 6, 12) * 0.20
        + range_score(branches, 0, 7) * 0.20
    )
    sa_penalty = 0.03 * branches + 0.04 * rings + 0.02 * halogens + max(0.0, heavy - 48) * 0.015
    sa = clamp(0.95 - sa_penalty)

    hydrophobic = counts["C"] + counts["F"] * 0.4 + counts["Cl"] * 0.5 + aromatic * 0.15
    hbond = min(hetero, 9)
    raw_docking = -(
        2.2
        + hydrophobic * 0.10
        + hbond * 0.26
        + min(aromatic, 18) * 0.04
    )
    raw_docking += max(0.0, 8 - heavy) * 0.18 + max(0.0, heavy - 55) * 0.08
    raw_docking += branches * 0.03
    docking = clamp((-raw_docking - 3.0) / 5.5)

    pose = valid * clamp(0.68 + hetero * 0.035 + aromatic * 0.005 - branches * 0.035 - max(0.0, heavy - 55) * 0.015)
    novelty = novelty_proxy(smiles, seeds)

    if not valid:
        lipinski *= 0.2
        qed *= 0.2
        sa *= 0.4
        docking *= 0.2
        pose = 0.0

    return {
        "validity": round(valid, 4),
        "qed": round(clamp(qed), 4),
        "sa": round(clamp(sa), 4),
        "lipinski": round(clamp(lipinski), 4),
        "docking": round(clamp(docking), 4),
        "pose": round(clamp(pose), 4),
        "novelty": round(clamp(novelty), 4),
        "raw_docking_kcal_mol": round(raw_docking, 3),
        "heavy_atoms_proxy": float(heavy),
        "hetero_atoms_proxy": float(hetero),
    }


def normalized_external_docking(value: str) -> Optional[Tuple[float, float]]:
    if value is None or value == "":
        return None
    try:
        raw = float(value)
    except ValueError:
        return None
    if 0.0 <= raw <= 1.0:
        return raw, raw
    if raw < 0:
        return clamp((-raw - 3.0) / 7.0), raw
    return clamp(raw / 100.0), raw


def parse_pose_value(value: str) -> Optional[float]:
    if value is None or value == "":
        return None
    raw = value.strip().lower()
    if raw in {"true", "pass", "passed", "yes", "y", "1"}:
        return 1.0
    if raw in {"false", "fail", "failed", "no", "n", "0"}:
        return 0.0
    try:
        return clamp(float(raw))
    except ValueError:
        return None


def load_external_scores(path: Optional[str]) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
    if not path:
        return {}, {}
    rows = read_csv(Path(path).expanduser().resolve())
    by_id = {row["id"]: row for row in rows if row.get("id")}
    by_smiles = {row["smiles"]: row for row in rows if row.get("smiles")}
    return by_id, by_smiles


def load_real_descriptors(
    project: Path,
    round_no: int,
    path: Optional[str],
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
    descriptor_path = Path(path).expanduser().resolve() if path else project / "stage4" / f"round_{round_no}_real_descriptors.csv"
    rows = read_csv(descriptor_path)
    by_id = {row["id"]: row for row in rows if row.get("id")}
    by_smiles = {}
    for row in rows:
        for key in ["canonical_smiles", "smiles"]:
            value = row.get(key, "")
            if value:
                by_smiles[value] = row
    return by_id, by_smiles


def apply_real_descriptor_scores(score: Dict[str, float], real: Dict[str, str]) -> bool:
    if not real:
        return False
    try:
        valid = float(real.get("valid", "0") or 0.0)
    except ValueError:
        valid = 0.0
    score["validity"] = clamp(valid)
    if valid <= 0.0:
        score["qed"] = 0.0
        score["lipinski"] = 0.0
        score["pose"] = 0.0
        score["docking"] *= 0.2
        score["heavy_atoms_proxy"] = 0.0
        score["hetero_atoms_proxy"] = 0.0
        return True

    try:
        score["qed"] = clamp(float(real.get("qed", "") or score.get("qed", 0.0)))
    except ValueError:
        pass
    try:
        violations = int(float(real.get("lipinski_violations", "0") or 0))
        score["lipinski"] = clamp(1.0 - violations / 4.0)
    except ValueError:
        pass
    try:
        heavy_atoms = float(real.get("heavy_atoms", "") or score.get("heavy_atoms_proxy", 0.0))
        score["heavy_atoms_proxy"] = heavy_atoms
    except ValueError:
        pass
    return True


def weighted_total(score: Dict[str, float], weights: Dict[str, object]) -> float:
    total_weight = 0.0
    total = 0.0
    for key, raw_weight in weights.items():
        try:
            weight = float(raw_weight)
        except (TypeError, ValueError):
            continue
        total += score.get(key, 0.0) * weight
        total_weight += weight
    if total_weight <= 0:
        return 0.0
    return clamp(total / total_weight)


def score_candidates(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    config = load_config(project)
    candidates_path = project / "candidates" / f"round_{args.round}_candidates.csv"
    rows = read_csv(candidates_path)
    if not rows:
        raise SystemExit(f"No candidates found: {candidates_path}")

    seeds = load_seed_smiles(project, args.round)
    external_by_id, external_by_smiles = load_external_scores(args.external_scores)
    real_by_id, real_by_smiles = load_real_descriptors(project, args.round, getattr(args, "real_descriptors", None))
    weights = config.get("weights", DEFAULT_WEIGHTS)
    if not isinstance(weights, dict):
        weights = DEFAULT_WEIGHTS

    scored = []
    for row in rows:
        smiles = row.get("smiles", "").strip()
        score = estimate_scores(smiles, seeds)
        note = "proxy"
        real = real_by_id.get(row.get("id", "")) or real_by_smiles.get(smiles)
        if real and apply_real_descriptor_scores(score, real):
            note = "proxy+stage4_rdkit"
        external = external_by_id.get(row.get("id", "")) or external_by_smiles.get(smiles)
        if external:
            docking_pair = normalized_external_docking(
                external.get("docking_score")
                or external.get("affinity")
                or external.get("gnina_score")
                or external.get("vina_score")
            )
            if docking_pair:
                score["docking"], score["raw_docking_kcal_mol"] = docking_pair
                note += "+external_docking"
            pose_value = parse_pose_value(
                external.get("pose_pass")
                or external.get("pose_score")
                or external.get("posebusters_pass")
            )
            if pose_value is not None:
                score["pose"] = pose_value
                note += "+external_pose"

        total = weighted_total(score, weights)
        scored.append(
            {
                "round": row.get("round", args.round),
                "id": row.get("id", ""),
                "smiles": smiles,
                "parent": row.get("parent", ""),
                "source": row.get("source", ""),
                "validity_proxy": score["validity"],
                "qed_proxy": score["qed"],
                "sa_proxy": score["sa"],
                "lipinski_proxy": score["lipinski"],
                "docking_proxy": score["docking"],
                "pose_proxy": score["pose"],
                "novelty_proxy": score["novelty"],
                "raw_docking_kcal_mol": score["raw_docking_kcal_mol"],
                "heavy_atoms_proxy": int(score["heavy_atoms_proxy"]),
                "hetero_atoms_proxy": int(score["hetero_atoms_proxy"]),
                "total_proxy": round(total, 4),
                "score_source": note,
            }
        )

    output = project / "scores" / f"round_{args.round}_scores.csv"
    write_csv(
        output,
        scored,
        [
            "round",
            "id",
            "smiles",
            "parent",
            "source",
            "validity_proxy",
            "qed_proxy",
            "sa_proxy",
            "lipinski_proxy",
            "docking_proxy",
            "pose_proxy",
            "novelty_proxy",
            "raw_docking_kcal_mol",
            "heavy_atoms_proxy",
            "hetero_atoms_proxy",
            "total_proxy",
            "score_source",
        ],
    )
    print(f"Wrote scores: {output}")
    return output


def float_field(row: Dict[str, str], field: str, default: float = 0.0) -> float:
    try:
        return float(row.get(field, default))
    except (TypeError, ValueError):
        return default


def rank_candidates(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    config = load_config(project)
    thresholds = config.get("thresholds", {})
    if not isinstance(thresholds, dict):
        thresholds = {}
    advance_total = float(thresholds.get("advance_total", 0.64))
    pose_min = float(thresholds.get("pose_min", 0.50))

    scores_path = project / "scores" / f"round_{args.round}_scores.csv"
    rows = read_csv(scores_path)
    if not rows:
        raise SystemExit(f"No scores found: {scores_path}")

    rows.sort(
        key=lambda row: (
            float_field(row, "total_proxy"),
            float_field(row, "docking_proxy"),
            float_field(row, "pose_proxy"),
        ),
        reverse=True,
    )

    ranked = []
    for idx, row in enumerate(rows, start=1):
        total = float_field(row, "total_proxy")
        pose = float_field(row, "pose_proxy")
        decision = "advance" if idx <= args.top and total >= advance_total and pose >= pose_min else "hold"
        ranked_row = dict(row)
        ranked_row["rank"] = idx
        ranked_row["decision"] = decision
        ranked.append(ranked_row)

    output = project / "ranked" / f"round_{args.round}_ranked.csv"
    fields = ["rank"] + [field for field in rows[0].keys() if field != "rank"] + ["decision"]
    write_csv(output, ranked, fields)

    advanced = [row for row in ranked if row["decision"] == "advance"]
    summary = render_round_summary(args.round, ranked, advanced, advance_total, pose_min)
    write_text(project / "reports" / f"round_{args.round}_summary.md", summary)
    print(f"Wrote ranked list: {output}")
    print(f"Advanced {len(advanced)} / {len(ranked)} candidates")
    return output


def render_round_summary(
    round_no: int,
    ranked: Sequence[Dict[str, object]],
    advanced: Sequence[Dict[str, object]],
    advance_total: float,
    pose_min: float,
) -> str:
    top_rows = ranked[:10]
    lines = [
        f"# Round {round_no} summary",
        "",
        "This report uses proxy scores unless external docking or pose scores were imported.",
        "",
        f"- Total candidates: {len(ranked)}",
        f"- Advanced candidates: {len(advanced)}",
        f"- Advance threshold: total >= {advance_total:.2f}, pose >= {pose_min:.2f}",
        "",
        "## Top candidates",
        "",
        "| rank | id | smiles | total | docking | pose | decision |",
        "|---:|---|---|---:|---:|---:|---|",
    ]
    for row in top_rows:
        lines.append(
            "| {rank} | {id} | `{smiles}` | {total_proxy} | {docking_proxy} | {pose_proxy} | {decision} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Next validation suggestions",
            "",
            "1. Export the top `advance` molecules to SDF and dock them with Vina or GNINA.",
            "2. Run PoseBusters on docked poses before trusting any pose-dependent score.",
            "3. Send only the best few families to OpenFE or a higher-cost physics simulation.",
            "4. Feed the accepted SMILES back into REINVENT4 or DrugEx as seeds, rewards, or fine-tuning data.",
        ]
    )
    return "\n".join(lines) + "\n"


def feedback(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    load_config(project)
    ranked_path = project / "ranked" / f"round_{args.round}_ranked.csv"
    rows = read_csv(ranked_path)
    if not rows:
        raise SystemExit(f"No ranked file found: {ranked_path}")
    rows.sort(key=lambda row: int(float_field(row, "rank", 999999)))

    selected = [
        row
        for row in rows
        if row.get("decision") == "advance"
    ][: args.top]
    if not selected:
        selected = rows[: args.top]

    seed_rows = [
        {
            "id": row.get("id", ""),
            "smiles": row.get("smiles", ""),
            "note": f"round_{args.round}_rank_{row.get('rank', '')}_total_{row.get('total_proxy', '')}",
        }
        for row in selected
    ]
    seed_path = project / "seeds" / f"round_{args.round}_seeds.csv"
    write_csv(seed_path, seed_rows, ["id", "smiles", "note"])

    feedback_doc = {
        "round": args.round,
        "selected_count": len(seed_rows),
        "next_seed_file": str(seed_path),
        "selected_smiles": [row["smiles"] for row in seed_rows],
        "recommended_next_steps": [
            "Use selected_smiles as a seed or fine-tuning set for the next generator round.",
            "Replace proxy docking with Vina/GNINA or DockStream results through score --external-scores.",
            "Use PoseBusters pass/fail as pose_pass in the external scores CSV.",
            "Escalate only diverse top families to OpenFE or wet-lab validation.",
        ],
    }
    output = project / "feedback" / f"round_{args.round}_feedback.json"
    write_json(output, feedback_doc)
    print(f"Wrote next seeds: {seed_path}")
    print(f"Wrote feedback package: {output}")
    return output


def fetch_url_text(url: str, timeout: int = 30, max_bytes: int = 800_000) -> Tuple[str, Optional[str]]:
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme not in {"http", "https"}:
        return "", "Only http/https URLs are supported."
    request = urllib.request.Request(url, headers={"User-Agent": "ai-mol-loop/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            raw = response.read(max_bytes + 1)
            if len(raw) > max_bytes:
                return "", f"URL content exceeds max bytes: {max_bytes}"
            charset = response.headers.get_content_charset() or "utf-8"
        return raw.decode(charset, errors="replace"), None
    except Exception as exc:
        return "", str(exc)


def extract_json_candidates(data: object) -> List[Dict[str, object]]:
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if not isinstance(data, dict):
        return []
    for key in ["candidates", "molecules", "rows", "items", "data", "results"]:
        value = data.get(key)
        if isinstance(value, list):
            return [item for item in value if isinstance(item, dict)]
    if data.get("smiles"):
        return [data]
    return []


def parse_csv_candidates(text: str) -> List[Dict[str, object]]:
    sample = text[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample)
    except csv.Error:
        dialect = csv.excel
    try:
        rows = [dict(row) for row in csv.DictReader(io.StringIO(text), dialect=dialect)]
    except csv.Error:
        return []
    if not rows:
        return []
    header = {field.lower().strip() for field in rows[0].keys() if field}
    if not header & {"smiles", "smile", "canonical_smiles", "mol", "structure"}:
        return []
    return rows


def parse_markdown_table_candidates(text: str) -> List[Dict[str, object]]:
    lines = [line.strip() for line in text.splitlines() if line.strip().startswith("|")]
    if len(lines) < 3:
        return []
    rows: List[Dict[str, object]] = []
    for idx in range(len(lines) - 2):
        header = [part.strip().lower() for part in lines[idx].strip("|").split("|")]
        separator = lines[idx + 1]
        if "smiles" not in header or not re.fullmatch(r"[\|\-\:\s]+", separator):
            continue
        for raw in lines[idx + 2 :]:
            parts = [part.strip() for part in raw.strip("|").split("|")]
            if len(parts) != len(header):
                break
            rows.append(dict(zip(header, parts)))
        return rows
    return []


def parse_smiles_from_text(text: str) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    seen = set()
    token_pattern = re.compile(r"(?<![A-Za-z0-9@+\-\[\]\(\)=#\\/%.])([A-Za-z0-9@+\-\[\]\(\)=#\\/%.]{2,180})(?![A-Za-z0-9@+\-\[\]\(\)=#\\/%.])")
    for match in token_pattern.finditer(text):
        token = match.group(1).strip("` ,.;:")
        if token in seen:
            continue
        if plausible_smiles(token):
            seen.add(token)
            rows.append({"smiles": token, "source": "text_smiles_extraction"})
        if len(rows) >= 500:
            break
    return rows


def parse_candidate_text(text: str) -> List[Dict[str, object]]:
    stripped = text.strip()
    if not stripped:
        return []
    try:
        return extract_json_candidates(json.loads(stripped))
    except json.JSONDecodeError:
        pass

    for match in re.finditer(r"```(?:json|csv)?\s*(.*?)```", stripped, flags=re.DOTALL | re.IGNORECASE):
        block = match.group(1).strip()
        try:
            rows = extract_json_candidates(json.loads(block))
            if rows:
                return rows
        except json.JSONDecodeError:
            pass
        rows = parse_csv_candidates(block)
        if rows:
            return rows

    rows = parse_csv_candidates(stripped)
    if rows:
        return rows
    rows = parse_markdown_table_candidates(stripped)
    if rows:
        return rows
    return parse_smiles_from_text(stripped)


def load_candidate_source(path: Path) -> List[Dict[str, object]]:
    text = path.read_text(encoding="utf-8-sig")
    return parse_candidate_text(text)


def candidate_smiles(row: Dict[str, object]) -> str:
    for key in ["smiles", "SMILES", "smile", "canonical_smiles", "canonical", "mol", "structure"]:
        value = row.get(key)
        if value:
            return str(value).strip()
    return ""


def normalize_candidate_rows(
    rows: Sequence[Dict[str, object]],
    round_no: int,
    source_label: str,
    start_index: int = 1,
) -> List[Dict[str, object]]:
    normalized: List[Dict[str, object]] = []
    for raw_idx, row in enumerate(rows, start=start_index):
        smiles = candidate_smiles(row)
        if not smiles:
            continue
        normalized.append(
            {
                "round": round_no,
                "id": row.get("id") or row.get("name") or f"r{round_no:02d}_stage3_{raw_idx:05d}",
                "smiles": smiles,
                "parent": row.get("parent", ""),
                "source": row.get("source") or source_label,
                "rationale": row.get("rationale", ""),
                "expected_interaction": row.get("expected_interaction", ""),
                "design_family": row.get("design_family", ""),
                "risk_note": row.get("risk_note", ""),
            }
        )
    return normalized


def approximate_molecular_weight(counts: Counter) -> float:
    weights = {
        "B": 10.81,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "F": 18.998,
        "P": 30.974,
        "S": 32.06,
        "Cl": 35.45,
        "Br": 79.904,
        "I": 126.904,
    }
    return sum(counts.get(atom, 0) * weight for atom, weight in weights.items())


def risk_flags_for_smiles(smiles: str) -> List[str]:
    alerts = [
        ("acid_chloride", ["C(=O)Cl", "C(=O)Br"]),
        ("isocyanate_or_isothiocyanate", ["N=C=O", "N=C=S"]),
        ("azide", ["N=[N+]=[N-]", "[N-]=[N+]=N"]),
        ("diazo", ["N=N=N", "[N-]=[N+]=N"]),
        ("peroxide", ["OO", "O-O"]),
        ("nitroso", ["N=O"]),
        ("sulfonyl_halide", ["S(=O)(=O)Cl", "S(=O)(=O)Br"]),
        ("epoxide_like", ["C1OC1"]),
    ]
    flags = []
    for name, patterns in alerts:
        if any(pattern in smiles for pattern in patterns):
            flags.append(name)
    if smiles.count("[") > 8:
        flags.append("many_explicit_atoms")
    if len(smiles) > 140:
        flags.append("very_long_smiles")
    return flags


def rdkit_descriptor_row(smiles: str) -> Optional[Dict[str, object]]:
    try:
        from rdkit import Chem  # type: ignore
        from rdkit import RDLogger  # type: ignore
        from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors  # type: ignore
    except Exception:
        return None
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "valid": 0,
            "canonical_smiles": "",
            "descriptor_source": "rdkit",
            "mw": "",
            "logp": "",
            "tpsa": "",
            "hbd": "",
            "hba": "",
            "rotatable_bonds": "",
            "heavy_atoms": "",
        }
    return {
        "valid": 1,
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "descriptor_source": "rdkit",
        "mw": round(float(Descriptors.MolWt(mol)), 3),
        "logp": round(float(Crippen.MolLogP(mol)), 3),
        "tpsa": round(float(rdMolDescriptors.CalcTPSA(mol)), 3),
        "hbd": int(rdMolDescriptors.CalcNumHBD(mol)),
        "hba": int(rdMolDescriptors.CalcNumHBA(mol)),
        "rotatable_bonds": int(rdMolDescriptors.CalcNumRotatableBonds(mol)),
        "heavy_atoms": int(mol.GetNumHeavyAtoms()),
    }


def fallback_descriptor_row(smiles: str) -> Dict[str, object]:
    valid = 1 if plausible_smiles(smiles) else 0
    counts = count_atoms(smiles)
    heavy = sum(counts.values())
    hetero = counts["N"] + counts["O"] + counts["S"] + counts["P"]
    halogens = counts["F"] + counts["Cl"] + counts["Br"] + counts["I"]
    mw = approximate_molecular_weight(counts)
    hbd = min(counts["N"] + counts["O"], 8)
    hba = min(hetero + counts["F"], 12)
    logp = round(0.34 * counts["C"] + 0.18 * halogens - 0.45 * hetero, 3)
    tpsa = round(17.0 * (counts["N"] + counts["O"]) + 24.0 * (counts["S"] + counts["P"]), 3)
    rotatable = max(0, smiles.count("C") + smiles.count("N") + smiles.count("O") - smiles.count("c") - 4)
    return {
        "valid": valid,
        "canonical_smiles": smiles if valid else "",
        "descriptor_source": "proxy",
        "mw": round(mw, 3),
        "logp": logp,
        "tpsa": tpsa,
        "hbd": int(hbd),
        "hba": int(hba),
        "rotatable_bonds": int(rotatable),
        "heavy_atoms": int(heavy),
    }


def lipinski_violations(desc: Dict[str, object]) -> int:
    violations = 0
    try:
        if float(desc.get("mw", 9999)) > 500:
            violations += 1
        if float(desc.get("logp", 9999)) > 5:
            violations += 1
        if int(float(desc.get("hbd", 999))) > 5:
            violations += 1
        if int(float(desc.get("hba", 999))) > 10:
            violations += 1
    except (TypeError, ValueError):
        return 4
    return violations


def filter_candidate_rows(
    rows: Sequence[Dict[str, object]],
    max_heavy_atoms: int,
    max_molecular_weight: float,
    max_lipinski_violations: int,
    allow_risk: bool,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    filtered: List[Dict[str, object]] = []
    accepted: List[Dict[str, object]] = []
    seen = set()
    for row in rows:
        smiles = str(row.get("smiles", "")).strip()
        desc = rdkit_descriptor_row(smiles) or fallback_descriptor_row(smiles)
        canonical = str(desc.get("canonical_smiles", ""))
        risk_flags = risk_flags_for_smiles(smiles)
        violations = lipinski_violations(desc)
        reasons = []
        if not desc.get("valid"):
            reasons.append("invalid_smiles")
        if canonical and canonical in seen:
            reasons.append("duplicate")
        try:
            if int(float(desc.get("heavy_atoms", 9999))) > max_heavy_atoms:
                reasons.append("too_many_heavy_atoms")
            if float(desc.get("mw", 9999)) > max_molecular_weight:
                reasons.append("mw_above_limit")
        except (TypeError, ValueError):
            reasons.append("descriptor_parse_failed")
        if violations > max_lipinski_violations:
            reasons.append("lipinski_violations")
        if risk_flags and not allow_risk:
            reasons.append("structural_risk_flags")

        passed = not reasons
        out = dict(row)
        out.update(desc)
        out["canonical_smiles"] = canonical
        out["lipinski_violations"] = violations
        out["risk_flags"] = "; ".join(risk_flags)
        out["passed_filter"] = "true" if passed else "false"
        out["filter_reason"] = "; ".join(reasons)
        filtered.append(out)
        if passed:
            seen.add(canonical or smiles)
            accepted.append(out)
    return filtered, accepted


def stage3_candidate_fields(rows: Sequence[Dict[str, object]]) -> List[str]:
    preferred = [
        "round",
        "id",
        "smiles",
        "canonical_smiles",
        "parent",
        "source",
        "valid",
        "passed_filter",
        "filter_reason",
        "descriptor_source",
        "mw",
        "logp",
        "tpsa",
        "hbd",
        "hba",
        "rotatable_bonds",
        "heavy_atoms",
        "lipinski_violations",
        "risk_flags",
        "rationale",
        "expected_interaction",
        "design_family",
        "risk_note",
    ]
    extra = []
    for row in rows:
        for key in row.keys():
            if key not in preferred and key not in extra:
                extra.append(key)
    return preferred + extra


def accepted_candidate_rows(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    return [
        {
            "round": row.get("round", ""),
            "id": row.get("id", ""),
            "smiles": row.get("canonical_smiles") or row.get("smiles", ""),
            "parent": row.get("parent", ""),
            "source": row.get("source", ""),
            "rationale": row.get("rationale", ""),
            "expected_interaction": row.get("expected_interaction", ""),
            "design_family": row.get("design_family", ""),
            "risk_note": row.get("risk_note", ""),
        }
        for row in rows
    ]


def resolve_openai_api_key(args: argparse.Namespace) -> Tuple[Optional[str], str]:
    if getattr(args, "api_key", None):
        return str(args.api_key).strip(), "--api-key"
    if getattr(args, "api_key_file", None):
        key_path = Path(str(args.api_key_file)).expanduser().resolve()
        return key_path.read_text(encoding="utf-8").strip(), "--api-key-file"
    env_name = getattr(args, "api_key_env", "OPENAI_API_KEY") or "OPENAI_API_KEY"
    key = os.environ.get(env_name, "").strip()
    if key:
        return key, f"env:{env_name}"
    return None, "missing"


def extract_openai_output_text(data: Dict[str, object]) -> str:
    direct = data.get("output_text")
    if isinstance(direct, str) and direct.strip():
        return direct
    chunks: List[str] = []
    output = data.get("output", [])
    if isinstance(output, list):
        for item in output:
            if not isinstance(item, dict):
                continue
            content = item.get("content", [])
            if isinstance(content, list):
                for part in content:
                    if isinstance(part, dict):
                        text = part.get("text")
                        if isinstance(text, str):
                            chunks.append(text)
            elif isinstance(content, str):
                chunks.append(content)
    return "\n".join(chunks).strip()


def openai_candidate_schema() -> Dict[str, object]:
    candidate = {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "id": {"type": "string"},
            "smiles": {"type": "string"},
            "rationale": {"type": "string"},
            "expected_interaction": {"type": "string"},
            "design_family": {"type": "string"},
            "risk_note": {"type": "string"},
        },
        "required": ["id", "smiles", "rationale", "expected_interaction", "design_family", "risk_note"],
    }
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "candidates": {
                "type": "array",
                "minItems": 1,
                "items": candidate,
            }
        },
        "required": ["candidates"],
    }


def openai_generate_candidates(
    prompt: str,
    api_key: str,
    model: str,
    timeout: int,
) -> List[Dict[str, object]]:
    payload = {
        "model": model,
        "input": [
            {
                "role": "system",
                "content": [
                    {
                        "type": "input_text",
                        "text": (
                            "You generate hypothetical small-molecule SMILES for computational virtual screening only. "
                            "Do not provide synthesis procedures, dosing, clinical claims, pathogen engineering guidance, or wet-lab instructions."
                        ),
                    }
                ],
            },
            {"role": "user", "content": [{"type": "input_text", "text": prompt}]},
        ],
        "text": {
            "format": {
                "type": "json_schema",
                "name": "molecule_candidates",
                "strict": True,
                "schema": openai_candidate_schema(),
            }
        },
    }
    body = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(
        "https://api.openai.com/v1/responses",
        data=body,
        headers={
            "Content-Type": "application/json",
            "Authorization": f"Bearer {api_key}",
            "User-Agent": "ai-mol-loop/0.1",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            data = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"OpenAI API request failed: HTTP {exc.code}: {detail[:800]}") from exc
    except Exception as exc:
        raise SystemExit(f"OpenAI API request failed: {exc}") from exc

    text = extract_openai_output_text(data)
    if not text:
        raise SystemExit("OpenAI API response did not contain output text.")
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"OpenAI API response was not valid JSON: {exc}") from exc
    rows = extract_json_candidates(parsed)
    if not rows:
        raise SystemExit("OpenAI API response did not contain candidates.")
    return rows


def load_stage3_prompt_context(project: Path, args: argparse.Namespace) -> str:
    config = load_config(project)
    sections: List[str] = []
    prompt_path = config.get("generator_prompt_file")
    if prompt_path and Path(str(prompt_path)).exists():
        sections.append(Path(str(prompt_path)).read_text(encoding="utf-8")[:12000])
    brief_path = config.get("target_brief_file")
    if brief_path and Path(str(brief_path)).exists():
        sections.append(Path(str(brief_path)).read_text(encoding="utf-8")[:8000])
    stage2_assets = project / "evidence" / "stage2_closed_loop_assets.json"
    if stage2_assets.exists():
        sections.append(stage2_assets.read_text(encoding="utf-8")[:12000])
    for url in getattr(args, "context_url", []) or []:
        text, error = fetch_url_text(url, args.timeout, args.max_url_bytes)
        if error:
            sections.append(f"Context URL fetch failed for {url}: {error}")
        else:
            sections.append(f"Context URL: {url}\n{text[:12000]}")
    user_prompt = getattr(args, "prompt", "")
    if user_prompt:
        sections.append(f"User requirement:\n{user_prompt}")
    return "\n\n---\n\n".join(sections)


def stage3_safe_json(path: Path) -> Dict[str, object]:
    if not path.exists():
        return {}
    try:
        data = read_json(path)
    except (json.JSONDecodeError, OSError):
        return {}
    return data if isinstance(data, dict) else {}


def stage3_prior_round_numbers(project: Path, current_round: int, max_rounds: int = 2) -> List[int]:
    available = set()
    for directory in ["ranked", "feedback", "filtered", "stage4", "stage4_5", "stage4_6"]:
        root = project / directory
        if not root.exists():
            continue
        for path in root.glob("round_*"):
            match = re.search(r"round_(\d+)", path.name)
            if not match:
                continue
            round_no = int(match.group(1))
            if 1 <= round_no < current_round:
                available.add(round_no)
    return sorted(available, reverse=True)[:max_rounds]


def stage3_bool_false(value: object) -> bool:
    return str(value or "").strip().lower() in {"0", "false", "fail", "failed", "no", "n"}


def stage3_row_subset(row: Dict[str, object], fields: Sequence[str]) -> Dict[str, object]:
    return {field: row.get(field, "") for field in fields if field in row}


def stage3_limited_rows(rows: Sequence[Dict[str, object]], fields: Sequence[str], limit: int) -> List[Dict[str, object]]:
    return [stage3_row_subset(row, fields) for row in list(rows)[:limit]]


def stage3_confidence_for_score_source(source: str, docking_score: object = "", pose_pass: object = "") -> str:
    text = str(source or "").lower()
    has_docking = docking_score not in {"", None}
    if "external_docking" in text or "vina" in text or "gnina" in text or has_docking:
        if stage3_bool_false(pose_pass):
            return "medium: real docking score present but pose quality failed or is uncertain"
        return "high_screening: real docking/pose evidence when generated by a configured backend"
    if "stage4_rdkit" in text or "rdkit" in text:
        return "medium: real RDKit descriptors plus proxy ranking"
    if "proxy" in text:
        return "low: proxy-only heuristic used for wiring and early triage"
    return "unknown: inspect source files before relying on this signal"


def stage3_prior_round_summary(project: Path, round_no: int, top: int = 8) -> Dict[str, object]:
    ranked_rows = read_csv(project / "ranked" / f"round_{round_no}_ranked.csv")
    ranked_rows.sort(key=lambda row: int(float_field(row, "rank", 999999)))
    filtered_rows = read_csv(project / "filtered" / f"round_{round_no}_filtered.csv")
    failed_rows = [
        row
        for row in filtered_rows
        if stage3_bool_false(row.get("passed_filter", "")) or str(row.get("filter_reason", "")).strip()
    ]
    docking_rows = read_csv(project / "stage4" / f"round_{round_no}_docking_scores_template.csv")
    docking_rows.sort(
        key=lambda row: (
            float(row.get("docking_score", "999999") or 999999)
            if re.fullmatch(r"-?\d+(?:\.\d+)?", str(row.get("docking_score", "")))
            else 999999,
            str(row.get("id", "")),
        )
    )
    feedback_doc = stage3_safe_json(project / "feedback" / f"round_{round_no}_feedback.json")
    stage45 = stage3_safe_json(project / "stage4_5" / f"round_{round_no}_control_validation.json")
    stage46 = stage3_safe_json(project / "stage4_6" / f"round_{round_no}_retrospective_benchmark.json")
    score_sources = sorted({str(row.get("score_source", "")) for row in ranked_rows if row.get("score_source")})
    return {
        "round": round_no,
        "round_label": f"round_{round_no}",
        "ranked_top": stage3_limited_rows(
            ranked_rows,
            ["rank", "id", "smiles", "total_proxy", "docking_proxy", "pose_proxy", "score_source", "decision"],
            top,
        ),
        "failed_filter_examples": stage3_limited_rows(
            failed_rows,
            ["id", "smiles", "passed_filter", "filter_reason", "risk_flags", "mw", "heavy_atoms", "lipinski_violations"],
            top,
        ),
        "feedback_constraints": {
            "selected_count": feedback_doc.get("selected_count", ""),
            "selected_smiles": list(feedback_doc.get("selected_smiles", [])[:top]) if isinstance(feedback_doc.get("selected_smiles"), list) else [],
            "recommended_next_steps": list(feedback_doc.get("recommended_next_steps", [])[:top]) if isinstance(feedback_doc.get("recommended_next_steps"), list) else [],
        },
        "stage4_docking_top": stage3_limited_rows(
            docking_rows,
            ["id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"],
            top,
        ),
        "stage45_control_calibration": {
            "best_known_control": stage45.get("best_known_control", {}),
            "best_candidate": stage45.get("best_candidate", {}),
            "best_decoy": stage45.get("best_decoy", {}),
            "candidate_vs_controls": stage45.get("candidate_vs_controls", {}),
            "docking": stage45.get("docking", {}),
        },
        "stage46_retrospective_benchmark": {
            "metrics": stage46.get("metrics", {}),
            "counts": stage46.get("counts", {}),
        },
        "source_confidence": [
            {
                "score_source": source,
                "confidence": stage3_confidence_for_score_source(source),
            }
            for source in score_sources
        ],
    }


def stage3_prompt_context_audit_policy() -> Dict[str, object]:
    return {
        "claim_boundary": [
            "Do not claim clinical efficacy, potency, safety, toxicity, dosing, or patient usefulness.",
            "Do not call a generated candidate a drug, treatment, cure, or clinically validated hit.",
            "Do not present docking/proxy/RDKit/control-calibration scores as measured biological activity.",
            "Do not infer safety or selectivity from similarity to a known drug.",
            "Do not provide synthesis procedures, experimental protocols, dosing, or wet-lab optimization instructions.",
        ],
        "retry_instruction": (
            "If requested output would violate the boundary, return safer computational-screening candidates "
            "and put the issue in risk_note instead of making the claim."
        ),
    }


def stage3_prior_feedback_context(project: Path, current_round: int, max_rounds: int = 2, top: int = 8) -> str:
    rounds = stage3_prior_round_numbers(project, current_round, max_rounds)
    if not rounds:
        return ""
    payload = {
        "section": "Prior Round Feedback Context",
        "purpose": (
            "Use these prior scores, rejection reasons, controls, decoys, and feedback constraints "
            "to guide the next generation round. Treat all signals as computational screening evidence."
        ),
        "rounds": [stage3_prior_round_summary(project, round_no, top) for round_no in rounds],
        "source_confidence_legend": {
            "low": "proxy-only heuristic; useful for wiring and diversity triage",
            "medium": "real RDKit descriptors or mixed proxy/descriptor ranking",
            "high_screening": "real docking/pose artifacts under configured computational settings; still not experimental evidence",
        },
        "generation_guidance": [
            "Prefer analog families related to advanced seeds while preserving diversity.",
            "Avoid repeating exact selected_smiles or near-duplicate known controls.",
            "Avoid motifs that previously failed filters, especially structural_risk_flags and excessive size.",
            "If prior docking/control calibration exists, use it as a ranking prior only, not as activity proof.",
        ],
        "compliance_audit": stage3_prompt_context_audit_policy(),
    }
    return "## Prior Round Feedback Context\n\n```json\n" + json.dumps(payload, indent=2, ensure_ascii=False)[:18000] + "\n```"


def build_openai_candidate_prompt(project: Path, args: argparse.Namespace) -> str:
    context = load_stage3_prompt_context(project, args)
    prior_context = stage3_prior_feedback_context(
        project,
        int(getattr(args, "round", 1) or 1),
        max_rounds=2,
        top=min(max(int(getattr(args, "top", 8) or 8), 1), 12),
    )
    context_sections = [context] if context else []
    if prior_context:
        context_sections.append(prior_context)
    context = "\n\n---\n\n".join(context_sections)
    return f"""Generate {args.n} diverse candidate molecules for the current computational screening project.

Return only structured JSON matching the schema.

Candidate requirements:
- hypothetical small-molecule SMILES only
- suitable for first-pass in silico screening
- avoid near-duplicate copies of known controls
- avoid obviously reactive, toxicophore-like, or oversized structures
- include high-level rationale and expected pocket interaction
- no synthesis, dosing, clinical, or wet-lab instructions
- use prior-round feedback as constraints when present, especially scores, hold/advance decisions, filter reasons, docking evidence, control calibration, decoys, and benchmark metrics
- label uncertainty in rationale/risk_note when a design is based on low or proxy source confidence
- if a requested statement would imply clinical value, rewrite it as computational prioritization only

Project context:
{context}
"""


def render_stage3_report(
    round_no: int,
    raw_rows: Sequence[Dict[str, object]],
    filtered_rows: Sequence[Dict[str, object]],
    accepted_rows: Sequence[Dict[str, object]],
    ranked_rows: Sequence[Dict[str, str]],
    sources: Sequence[str],
    used_openai: bool,
    model: str,
    api_key_source: str,
) -> str:
    failed = len(filtered_rows) - len(accepted_rows)
    lines = [
        f"# Stage 3 Candidate Input And Screening Report - Round {round_no}",
        "",
        "## Summary",
        "",
        f"- Raw candidates: {len(raw_rows)}",
        f"- Passed filter: {len(accepted_rows)}",
        f"- Failed filter: {failed}",
        f"- Sources: {', '.join(sources) if sources else 'proxy_generator'}",
        f"- OpenAI generation used: {used_openai}",
        f"- OpenAI model: {model if used_openai else 'not_used'}",
        f"- API key source: {api_key_source if used_openai else 'not_used'}",
        "",
        "## Filtered Candidates",
        "",
        "| id | smiles | pass | descriptor | MW | logP | HBD | HBA | reason |",
        "|---|---|---|---|---:|---:|---:|---:|---|",
    ]
    for row in filtered_rows[:20]:
        lines.append(
            f"| {row.get('id', '')} | `{row.get('smiles', '')}` | {row.get('passed_filter', '')} | "
            f"{row.get('descriptor_source', '')} | {row.get('mw', '')} | {row.get('logp', '')} | "
            f"{row.get('hbd', '')} | {row.get('hba', '')} | {row.get('filter_reason', '')} |"
        )
    if ranked_rows:
        lines.extend(
            [
                "",
                "## Top Ranked After Proxy/External Scoring",
                "",
                "| rank | id | smiles | total | docking | pose | decision |",
                "|---:|---|---|---:|---:|---:|---|",
            ]
        )
        for row in ranked_rows[:10]:
            lines.append(
                f"| {row.get('rank', '')} | {row.get('id', '')} | `{row.get('smiles', '')}` | "
                f"{row.get('total_proxy', '')} | {row.get('docking_proxy', '')} | {row.get('pose_proxy', '')} | {row.get('decision', '')} |"
            )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- This stage performs candidate intake, standardization, early filtering, and computational ranking only.",
            "- Proxy scores and docking imports are not measured potency or therapeutic proof.",
            "- OpenAI API keys are never written to project outputs; only the key source label is recorded.",
        ]
    )
    return "\n".join(lines) + "\n"


def stage3_screen(args: argparse.Namespace) -> None:
    project = project_path(args.project)
    ensure_project_dirs(project)
    load_config(project)

    raw_rows: List[Dict[str, object]] = []
    sources: List[str] = []
    openai_prompt_path = ""
    openai_prompt_context = {
        "prior_rounds": stage3_prior_round_numbers(project, int(getattr(args, "round", 1) or 1), 2),
        "uses_feedback_context": False,
        "saved_prompt": "",
        "boundary": "Prompt snapshots omit API keys and store computational-screening context only.",
    }

    if args.source_csv:
        path = Path(args.source_csv).expanduser().resolve()
        rows = read_csv(path)
        raw_rows.extend(normalize_candidate_rows(rows, args.round, f"csv:{path.name}", len(raw_rows) + 1))
        sources.append(f"csv:{path}")
    if args.source_json:
        path = Path(args.source_json).expanduser().resolve()
        rows = load_candidate_source(path)
        raw_rows.extend(normalize_candidate_rows(rows, args.round, f"json:{path.name}", len(raw_rows) + 1))
        sources.append(f"json:{path}")
    for url in args.source_url or []:
        text, error = fetch_url_text(url, args.timeout, args.max_url_bytes)
        if error:
            raise SystemExit(f"Failed to fetch source URL {url}: {error}")
        rows = parse_candidate_text(text)
        if rows:
            raw_rows.extend(normalize_candidate_rows(rows, args.round, f"url:{url}", len(raw_rows) + 1))
            sources.append(f"url:{url}")
        elif args.use_openai:
            sources.append(f"url_context:{url}")
            if not args.context_url:
                args.context_url = []
            args.context_url.append(url)
        else:
            raise SystemExit(f"No candidate rows found in URL: {url}")

    api_key_source = "not_used"
    if args.use_openai:
        api_key, api_key_source = resolve_openai_api_key(args)
        if not api_key:
            raise SystemExit("Missing OpenAI API key. Set OPENAI_API_KEY or pass --api-key-file/--api-key.")
        prompt = build_openai_candidate_prompt(project, args)
        openai_prompt_path = str(project / "prompts" / f"round_{args.round}_openai_candidate_prompt.md")
        write_text(Path(openai_prompt_path), prompt)
        openai_prompt_context["uses_feedback_context"] = bool(openai_prompt_context["prior_rounds"])
        openai_prompt_context["saved_prompt"] = openai_prompt_path
        openai_rows = openai_generate_candidates(prompt, api_key, args.openai_model, args.timeout)
        raw_rows.extend(normalize_candidate_rows(openai_rows, args.round, f"openai:{args.openai_model}", len(raw_rows) + 1))
        sources.append("openai")

    if not raw_rows:
        generate_candidates(argparse.Namespace(project=str(project), round=args.round, n=args.n, source_csv=None))
        raw_rows = normalize_candidate_rows(read_csv(project / "candidates" / f"round_{args.round}_candidates.csv"), args.round, "proxy_generator")
        sources.append("proxy_generator")

    raw_path = project / "stage3" / f"round_{args.round}_raw_candidates.csv"
    write_csv(raw_path, raw_rows, stage3_candidate_fields(raw_rows))

    filtered_rows, accepted_rows = filter_candidate_rows(
        raw_rows,
        args.max_heavy_atoms,
        args.max_molecular_weight,
        args.max_lipinski_violations,
        args.allow_risk,
    )
    filtered_path = project / "filtered" / f"round_{args.round}_filtered.csv"
    write_csv(filtered_path, filtered_rows, stage3_candidate_fields(filtered_rows))

    accepted_for_scoring = accepted_candidate_rows(accepted_rows)
    candidates_path = project / "candidates" / f"round_{args.round}_candidates.csv"
    write_csv(
        candidates_path,
        accepted_for_scoring,
        ["round", "id", "smiles", "parent", "source", "rationale", "expected_interaction", "design_family", "risk_note"],
    )

    ranked_rows: List[Dict[str, str]] = []
    if accepted_for_scoring and not args.no_score:
        score_candidates(argparse.Namespace(project=str(project), round=args.round, external_scores=args.external_scores))
        rank_candidates(argparse.Namespace(project=str(project), round=args.round, top=args.top))
        feedback(argparse.Namespace(project=str(project), round=args.round, top=args.top))
        ranked_rows = read_csv(project / "ranked" / f"round_{args.round}_ranked.csv")

    assets_path = project / "stage3" / f"round_{args.round}_stage3_assets.json"
    write_json(
        assets_path,
        {
            "schema_version": "0.1",
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "round": args.round,
            "sources": sources,
            "raw_candidates": len(raw_rows),
            "passed_filter": len(accepted_rows),
            "failed_filter": len(filtered_rows) - len(accepted_rows),
            "openai": {
                "used": bool(args.use_openai),
                "model": args.openai_model if args.use_openai else "",
                "api_key_source": api_key_source if args.use_openai else "",
                "prompt_file": openai_prompt_path,
                "context": openai_prompt_context,
            },
            "files": {
                "raw_candidates": str(raw_path),
                "filtered_candidates": str(filtered_path),
                "scoring_candidates": str(candidates_path),
                "ranked": str(project / "ranked" / f"round_{args.round}_ranked.csv") if ranked_rows else "",
            },
        },
    )
    report_path = project / "reports" / f"stage3_round_{args.round}_report.md"
    write_text(
        report_path,
        render_stage3_report(
            args.round,
            raw_rows,
            filtered_rows,
            accepted_rows,
            ranked_rows,
            sources,
            bool(args.use_openai),
            args.openai_model,
            api_key_source,
        ),
    )

    print(f"Wrote raw stage 3 candidates: {raw_path}")
    print(f"Wrote filtered candidates: {filtered_path}")
    print(f"Wrote scoring candidates: {candidates_path}")
    print(f"Wrote stage 3 assets: {assets_path}")
    print(f"Wrote stage 3 report: {report_path}")
    if not accepted_for_scoring:
        print("No candidates passed stage 3 filters; scoring was skipped.")


CONTROL_SMILES_BY_NAME: Dict[str, Dict[str, str]] = {
    "oseltamivir": {
        "smiles": "CCC(CC)O[C@@H]1C=C(C[C@@H]([C@H]1NC(=O)C)N)C(=O)OCC",
        "pubchem_cid": "65028",
        "source": "PubChem PUG-REST",
    },
    "zanamivir": {
        "smiles": "CC(=O)N[C@@H]1[C@H](C=C(O[C@H]1[C@@H]([C@@H](CO)O)O)C(=O)O)N=C(N)N",
        "pubchem_cid": "60855",
        "source": "PubChem PUG-REST",
    },
    "peramivir": {
        "smiles": "CCC(CC)[C@@H]([C@H]1[C@@H](C[C@@H]([C@H]1O)C(=O)O)N=C(N)N)NC(=O)C",
        "pubchem_cid": "154234",
        "source": "PubChem PUG-REST",
    },
    "laninamivir": {
        "smiles": "CC(=O)N[C@@H]1[C@H](C=C(O[C@H]1[C@@H]([C@@H](CO)O)OC)C(=O)O)N=C(N)N",
        "pubchem_cid": "502272",
        "source": "PubChem PUG-REST",
    },
    "baloxavir acid": {
        "smiles": "C1COC[C@@H]2N1C(=O)C3=C(C(=O)C=CN3N2[C@H]4C5=C(CSC6=CC=CC=C46)C(=C(C=C5)F)F)O",
        "pubchem_cid": "124081876",
        "source": "PubChem PUG-REST",
    },
    "baloxavir marboxil": {
        "smiles": "COC(=O)OCOC1=C2C(=O)N3CCOC[C@H]3N(N2C=CC1=O)[C@H]4C5=C(CSC6=CC=CC=C46)C(=C(C=C5)F)F",
        "pubchem_cid": "124081896",
        "source": "PubChem PUG-REST",
    },
    "amantadine": {
        "smiles": "C1C2CC3CC1CC(C2)(C3)N",
        "pubchem_cid": "2130",
        "source": "PubChem PUG-REST",
    },
    "rimantadine": {
        "smiles": "CC(C12CC3CC(C1)CC(C3)C2)N",
        "pubchem_cid": "5071",
        "source": "PubChem PUG-REST",
    },
}


STAGE4_DESCRIPTOR_FIELDS = [
    "round",
    "id",
    "smiles",
    "canonical_smiles",
    "parent",
    "source",
    "valid",
    "stage4_status",
    "descriptor_source",
    "mol_formula",
    "exact_mw",
    "mw",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "rotatable_bonds",
    "heavy_atoms",
    "rings",
    "aromatic_rings",
    "formal_charge",
    "fraction_csp3",
    "qed",
    "lipinski_violations",
    "murcko_scaffold",
    "risk_flags",
    "error",
]


STAGE4_SIMILARITY_FIELDS = [
    "candidate_id",
    "candidate_smiles",
    "control_drug",
    "control_role",
    "control_smiles",
    "tanimoto",
    "fingerprint",
    "control_source",
    "pubchem_cid",
]


STAGE4_DIVERSITY_FIELDS = [
    "selection_rank",
    "id",
    "canonical_smiles",
    "qed",
    "max_similarity_to_selected",
    "max_similarity_to_controls",
    "nearest_control",
    "murcko_scaffold",
]


STAGE4_BENCHMARK_FIELDS = [
    "panel_type",
    "id",
    "smiles",
    "canonical_smiles",
    "role",
    "target_id",
    "qed",
    "mw",
    "logp",
    "tpsa",
    "hbd",
    "hba",
    "heavy_atoms",
    "murcko_scaffold",
    "max_similarity_to_controls",
    "nearest_control",
    "source",
]


STAGE4_DOCKING_INPUT_FIELDS = [
    "panel_type",
    "id",
    "canonical_smiles",
    "ligand_sdf",
    "ligand_pdbqt",
    "receptor_pdb",
    "receptor_pdbqt",
    "status",
    "note",
]

STAGE4_DOCKING_SCORE_FIELDS = ["id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"]

STAGE45_DOCKING_SCORE_FIELDS = [
    "panel_type",
    "id",
    "smiles",
    "docking_score",
    "pose_pass",
    "backend",
    "receptor",
    "relative_to_best_control",
    "score_band",
    "notes",
]

STAGE46_RANKING_FIELDS = [
    "rank",
    "id",
    "panel_type",
    "label",
    "docking_score",
    "score_band",
    "pose_pass",
    "backend",
    "receptor",
    "relative_to_best_control",
    "relative_to_best_row",
    "smiles",
    "notes",
]

STAGE46_DEFAULT_POSITIVE_TYPES = ["positive_control", "reference_control", "control"]
STAGE46_DEFAULT_NEGATIVE_TYPES = ["decoy"]
STAGE46_DEFAULT_TOP_K = [1, 3, 5, 10]


DEFAULT_DECOY_SMILES = [
    ("decoy_001", "CCO", "small_polar_fragment"),
    ("decoy_002", "c1ccccc1", "simple_aromatic"),
    ("decoy_003", "CCN(CC)CC", "small_amine"),
    ("decoy_004", "CCOC(=O)c1ccccc1", "ester_aromatic"),
    ("decoy_005", "O=C(NCCO)c1ccccc1", "amide_aromatic"),
    ("decoy_006", "CCS(=O)(=O)Nc1ccccc1", "sulfonamide_aromatic"),
    ("decoy_007", "CC(C)NC(=O)c1ccccc1", "alkyl_amide"),
    ("decoy_008", "COc1ccccc1O", "phenoxy_polar"),
    ("decoy_009", "CC(C)Oc1ccccc1C#N", "nitrile_aromatic"),
    ("decoy_010", "c1ccncc1", "pyridine_fragment"),
    ("decoy_011", "CC(C)(C)OC(=O)N1CCOCC1", "protected_morpholine"),
    ("decoy_012", "CC(C)C(=O)Nc1ccccc1", "branched_anilide"),
]


def python_module_status(module_name: str) -> Dict[str, object]:
    spec = importlib.util.find_spec(module_name)
    status: Dict[str, object] = {
        "name": module_name,
        "status": "available" if spec else "missing",
        "origin": spec.origin if spec and spec.origin else "",
    }
    if spec:
        try:
            module = __import__(module_name)
            status["version"] = str(getattr(module, "__version__", ""))
        except Exception as exc:
            status["status"] = "import_error"
            status["error"] = str(exc)
    return status


def executable_search_paths() -> List[Path]:
    paths: List[Path] = []
    for raw in os.environ.get("PATH", "").split(os.pathsep):
        if raw:
            paths.append(Path(raw).expanduser())
    for raw in [
        sysconfig.get_path("scripts"),
        sysconfig.get_path("scripts", scheme="osx_framework_user"),
        str(Path(site.getuserbase()) / "bin"),
        str(DEFAULT_EVAL_TOOLS / "micromamba-root" / "envs" / "docking" / "bin"),
        str(DEFAULT_EVAL_TOOLS / "bin"),
    ]:
        if raw:
            paths.append(Path(raw).expanduser())
    unique: List[Path] = []
    seen = set()
    for path in paths:
        key = str(path)
        if key not in seen:
            seen.add(key)
            unique.append(path)
    return unique


def validate_executable_status(executable: str, status: Dict[str, str]) -> Dict[str, str]:
    if status.get("status") != "found" or executable != "vina":
        return status
    path = Path(status.get("path", ""))
    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return status
    if (
        text.startswith("#!/bin/sh")
        and "-7.4" in text
        and "--out" not in text
        and "printf" not in text
    ):
        updated = dict(status)
        updated["status"] = "invalid"
        updated["validation"] = "placeholder_no_pose_output"
        updated["message"] = "This vina executable is a placeholder script that prints a fixed score and cannot write docking pose files."
        return updated
    return status


def executable_status(executable: str) -> Dict[str, str]:
    path = shutil.which(executable)
    source = "PATH" if path else ""
    if not path:
        for directory in executable_search_paths():
            candidate = directory / executable
            if candidate.exists() and os.access(candidate, os.X_OK):
                path = str(candidate)
                source = str(directory)
                break
    return validate_executable_status(executable, {"name": executable, "status": "found" if path else "not_found", "path": path or "", "source": source})


def stage4_capabilities() -> Dict[str, object]:
    modules = {
        name: python_module_status(name)
        for name in ["rdkit", "numpy", "pandas", "openbabel", "meeko", "vina", "posebusters", "gemmi"]
    }
    executables = {
        name: executable_status(name)
        for name in [
            "vina",
            "gnina",
            "obabel",
            "mk_prepare_ligand.py",
            "mk_prepare_receptor.py",
            "bust",
            "openfe",
        ]
    }
    docking_backends = [
        name
        for name in ["vina", "gnina"]
        if executables.get(name, {}).get("status") == "found" or modules.get(name, {}).get("status") == "available"
    ]
    if docking_backends:
        docking_status = "available"
        docking_message = "At least one docking backend is installed and can be wired as the next scoring step."
    else:
        docking_status = "not_available"
        docking_message = "No Vina/GNINA backend found; Stage 4 produced docking-ready ligand assets but did not run docking."
    return {
        "modules": modules,
        "executables": executables,
        "rdkit": modules["rdkit"],
        "docking_backend": {
            "status": docking_status,
            "available_backends": docking_backends,
            "message": docking_message,
        },
        "ligand_preparation": {
            "openbabel": executables["obabel"],
            "meeko": modules["meeko"],
            "message": "Install OpenBabel or Meeko to convert SDF into PDBQT for AutoDock Vina.",
        },
        "install_options": {
            "pip_user": [
                "python3 -m pip install --user meeko posebusters openbabel-wheel gemmi prody",
                "python3 -m pip install --user vina",
            ],
            "binary_tools": [
                "Install AutoDock Vina or GNINA executable and make it visible on PATH.",
                f"Alternatively place executables under {DEFAULT_EVAL_TOOLS / 'bin'}.",
            ],
            "docker_tools": [
                "GNINA upstream release binaries are Linux/CUDA oriented; on macOS use Docker Desktop or a Linux workstation.",
                "Example Docker path: mount the project stage4 directory and run a gnina image against receptor/ligand files, then import the score CSV.",
            ],
            "notes": [
                "Vina Python wheels may be unavailable on macOS arm64 and can require Boost to build from source.",
                "GNINA usually ships as platform-specific binaries; macOS support depends on the released artifact.",
                "This macOS product path treats GNINA as external/Docker-backed unless a native gnina executable is explicitly found.",
            ],
        },
        "search_paths": [str(path) for path in executable_search_paths()],
        "boundary": [
            "Installed Python packages do not by themselves prove docking has run.",
            "Docking status is completed only after a backend writes score rows.",
        ],
    }


def require_rdkit() -> Dict[str, object]:
    try:
        from rdkit import Chem, DataStructs, RDLogger, rdBase  # type: ignore
        from rdkit.Chem import AllChem, Crippen, Descriptors, QED, rdMolDescriptors  # type: ignore
        from rdkit.Chem.Scaffolds import MurckoScaffold  # type: ignore
    except Exception as exc:
        raise SystemExit(
            "Stage 4 requires RDKit for the real-library path. "
            "Install RDKit or run earlier proxy-only stages."
        ) from exc
    RDLogger.DisableLog("rdApp.*")
    return {
        "Chem": Chem,
        "DataStructs": DataStructs,
        "AllChem": AllChem,
        "Crippen": Crippen,
        "Descriptors": Descriptors,
        "QED": QED,
        "rdMolDescriptors": rdMolDescriptors,
        "MurckoScaffold": MurckoScaffold,
        "rdBase": rdBase,
    }


def load_stage4_candidates(project: Path, round_no: int, input_csv: Optional[str]) -> List[Dict[str, str]]:
    if input_csv:
        path = Path(input_csv).expanduser().resolve()
    else:
        path = project / "candidates" / f"round_{round_no}_candidates.csv"
    rows = read_csv(path)
    if not rows:
        raise SystemExit(f"No candidate rows found for Stage 4: {path}")
    return rows


def load_stage4_controls(target_id: str, controls_csv: Optional[str]) -> List[Dict[str, object]]:
    rows = read_csv(Path(controls_csv).expanduser().resolve()) if controls_csv else read_csv(DEFAULT_KNOWN_DRUGS)
    controls: List[Dict[str, object]] = []
    seen = set()
    for row in rows:
        if target_id and row.get("target_id") and row.get("target_id") != target_id:
            continue
        drug = str(row.get("drug", "")).strip()
        if not drug:
            continue
        smiles = row.get("smiles", "").strip() if isinstance(row.get("smiles"), str) else ""
        meta = CONTROL_SMILES_BY_NAME.get(drug.lower())
        if not smiles and meta:
            smiles = meta["smiles"]
        if not smiles:
            continue
        key = (drug.lower(), smiles)
        if key in seen:
            continue
        seen.add(key)
        controls.append(
            {
                "drug": drug,
                "target_id": row.get("target_id", target_id),
                "role": row.get("status_for_workflow", row.get("role", "")),
                "mechanism": row.get("mechanism", ""),
                "smiles": smiles,
                "pubchem_cid": meta.get("pubchem_cid", "") if meta else row.get("pubchem_cid", ""),
                "source": meta.get("source", "controls_csv") if meta else row.get("source", "controls_csv"),
            }
        )
    return controls


def stage4_target_context(project: Path, target_id: str, config: Dict[str, object]) -> Dict[str, object]:
    context: Dict[str, object] = {"target_id": target_id}
    if target_id:
        try:
            catalog = load_target_catalog(None)
            target = get_target_by_id(catalog, target_id)
            context["target_catalog"] = target
        except SystemExit:
            context["target_catalog"] = {}
    stage2_assets = project / "evidence" / "stage2_closed_loop_assets.json"
    if stage2_assets.exists():
        try:
            assets = read_json(stage2_assets)
            targets = assets.get("targets", [])
            if isinstance(targets, list):
                for item in targets:
                    if isinstance(item, dict) and item.get("target_id") == target_id:
                        context["stage2_target"] = item
                        break
        except (json.JSONDecodeError, OSError):
            pass
    target_config = config.get("target", {})
    if isinstance(target_config, dict):
        context["config_target"] = target_config
    return context


def nested_get(data: Dict[str, object], keys: Sequence[str]) -> object:
    current: object = data
    for key in keys:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def stage4_primary_structure(context: Dict[str, object]) -> Dict[str, object]:
    stage2_primary = nested_get(context, ["stage2_target", "closed_loop_assets", "primary_structure"])
    if isinstance(stage2_primary, dict) and stage2_primary.get("pdb_id"):
        return dict(stage2_primary)
    catalog = context.get("target_catalog", {})
    if isinstance(catalog, dict):
        structures = catalog.get("representative_structures", [])
        if isinstance(structures, list) and structures and isinstance(structures[0], dict):
            return dict(structures[0])
    config_target = context.get("config_target", {})
    if isinstance(config_target, dict):
        pdb_id = config_target.get("pdb_id") or config_target.get("protein_pdb")
        if pdb_id:
            return {"pdb_id": pdb_id, "source": "project_config"}
    return {}


def stage4_binding_site(context: Dict[str, object]) -> Dict[str, object]:
    site: Dict[str, object] = {}
    stage2_site = nested_get(context, ["stage2_target", "closed_loop_assets", "binding_site"])
    if isinstance(stage2_site, dict) and stage2_site:
        site.update(stage2_site)
    if not site:
        catalog_site = nested_get(context, ["target_catalog", "binding_site"])
        if isinstance(catalog_site, dict) and catalog_site:
            site.update(catalog_site)
    config_site = nested_get(context, ["config_target", "pocket"])
    if isinstance(config_site, dict):
        site.update(config_site)
    return site


def fetch_pdb_file(pdb_id: str, output: Path, timeout: int = 30) -> Tuple[bool, str]:
    if not pdb_id:
        return False, "missing_pdb_id"
    url = f"https://files.rcsb.org/download/{urllib.parse.quote(pdb_id.upper())}.pdb"
    request = urllib.request.Request(url, headers={"User-Agent": "ai-mol-loop/0.1"})
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            raw = response.read()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(raw)
        return True, "fetched_from_rcsb"
    except Exception as exc:
        return False, str(exc)


PDB_NON_BINDING_HETERO_RESNAMES = {
    "HOH",
    "WAT",
    "DOD",
    "H2O",
    "NA",
    "CL",
    "K",
    "CA",
    "MG",
    "MN",
    "ZN",
    "FE",
    "CU",
    "NI",
    "CO",
    "CD",
    "HG",
    "LI",
    "RB",
    "CS",
    "SR",
    "BA",
    "AL",
    "F",
    "BR",
    "I",
    "SO4",
    "PO4",
    "EDO",
    "PEG",
    "PG4",
    "GOL",
    "DMS",
    "ACT",
    "ACY",
    "ACE",
    "TRS",
    "MES",
    "HEP",
    "BME",
    "MPD",
    "NAG",
    "BMA",
    "MAN",
    "FUC",
}

PDB_REFERENCE_LIGAND_ALIASES = {
    "oseltamivir": {"G39", "OTV", "OSE", "OSV", "GS4"},
    "zanamivir": {"ZMR", "ZAN", "G20"},
    "laninamivir": {"G28", "LAN"},
    "peramivir": {"PRV", "PRM"},
    "baloxavir acid": {"BXA", "BXM", "BAV"},
}

PDB_STANDARD_PROTEIN_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "MSE",
}


def stage4_has_binding_box(binding_site: Dict[str, object]) -> bool:
    center = binding_site.get("center") or binding_site.get("box_center")
    size = binding_site.get("size") or binding_site.get("box_size")
    if isinstance(center, str):
        center = parse_triplet(center)
    if isinstance(size, str):
        size = parse_triplet(size)
    return isinstance(center, list) and isinstance(size, list) and len(center) == 3 and len(size) == 3


def stage4_reference_ligand_aliases(reference_ligand: str) -> set:
    text = str(reference_ligand or "").strip().lower()
    aliases = set()
    if text:
        aliases.add(text.upper())
        aliases.update(PDB_REFERENCE_LIGAND_ALIASES.get(text, set()))
    return {str(item).upper() for item in aliases if str(item).strip()}


def parse_pdb_hetero_ligands(path: Path, reference_ligand: str = "") -> List[Dict[str, object]]:
    aliases = stage4_reference_ligand_aliases(reference_ligand)
    groups: Dict[Tuple[str, str, int, str], Dict[str, object]] = {}
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []

    for line in lines:
        if not line.startswith("HETATM"):
            continue
        alt_loc = line[16:17].strip() if len(line) > 16 else ""
        if alt_loc and alt_loc not in {"A", "1"}:
            continue
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        if not resname:
            continue
        chain = line[21:22].strip() if len(line) >= 22 else ""
        icode = line[26:27].strip() if len(line) >= 27 else ""
        try:
            resseq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except Exception:
            continue
        atom_name = line[12:16].strip() if len(line) >= 16 else ""
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element and atom_name:
            element = re.sub(r"[^A-Za-z]", "", atom_name)[:1].upper()
        if element == "H":
            continue
        key = (resname, chain, resseq, icode)
        group = groups.setdefault(
            key,
            {
                "resname": resname,
                "chain": chain,
                "resseq": resseq,
                "icode": icode,
                "atoms": [],
                "elements": Counter(),
            },
        )
        group["atoms"].append((x, y, z))
        elements = group.get("elements")
        if isinstance(elements, Counter):
            elements[element or "?"] += 1

    ligands: List[Dict[str, object]] = []
    for group in groups.values():
        resname = str(group.get("resname", "")).upper()
        atoms = group.get("atoms", [])
        if not isinstance(atoms, list) or not atoms:
            continue
        is_reference_alias = resname in aliases
        if not is_reference_alias and resname in PDB_NON_BINDING_HETERO_RESNAMES:
            continue
        if len(atoms) < 4 and not is_reference_alias:
            continue
        xs = [float(item[0]) for item in atoms]
        ys = [float(item[1]) for item in atoms]
        zs = [float(item[2]) for item in atoms]
        center = [round(sum(xs) / len(xs), 3), round(sum(ys) / len(ys), 3), round(sum(zs) / len(zs), 3)]
        bounds_min = [round(min(xs), 3), round(min(ys), 3), round(min(zs), 3)]
        bounds_max = [round(max(xs), 3), round(max(ys), 3), round(max(zs), 3)]
        extent = [round(bounds_max[i] - bounds_min[i], 3) for i in range(3)]
        atom_count = len(atoms)
        score = atom_count + (1000 if is_reference_alias else 0)
        ligands.append(
            {
                "resname": resname,
                "chain": group.get("chain", ""),
                "resseq": group.get("resseq", ""),
                "icode": group.get("icode", ""),
                "atom_count": atom_count,
                "center": center,
                "bounds_min": bounds_min,
                "bounds_max": bounds_max,
                "extent": extent,
                "reference_alias_match": is_reference_alias,
                "selection_score": score,
                "elements": dict(group.get("elements", {})),
            }
        )
    ligands.sort(key=lambda item: (-int(item.get("selection_score", 0)), str(item.get("resname", "")), str(item.get("chain", "")), int(item.get("resseq") or 0)))
    return ligands


def extract_stage4_cocrystal_pocket(
    receptor_pdb: Path,
    reference_ligand: str = "",
    min_size: float = 20.0,
    padding: float = 8.0,
) -> Dict[str, object]:
    ligands = parse_pdb_hetero_ligands(receptor_pdb, reference_ligand)
    if not ligands:
        return {
            "status": "no_ligand_detected",
            "source": "co_crystal_ligand",
            "message": "No non-solvent HETATM ligand was detected in the receptor PDB.",
        }
    selected = dict(ligands[0])
    extent = selected.get("extent", [])
    if not isinstance(extent, list) or len(extent) != 3:
        extent = [0.0, 0.0, 0.0]
    box_size = [round(max(float(min_size), float(value) + 2.0 * float(padding)), 3) for value in extent]
    return {
        "status": "extracted_from_receptor",
        "source": "co_crystal_ligand",
        "reference_ligand": reference_ligand,
        "detected_ligand": selected,
        "center": selected.get("center", []),
        "size": box_size,
        "box_padding_angstrom": float(padding),
        "box_min_size_angstrom": float(min_size),
        "candidate_ligands": ligands[:10],
        "message": "Docking box was computed from the selected co-crystallized ligand centroid and ligand bounds.",
    }


def clean_stage4_receptor_pdb_for_docking(receptor_pdb: Path, output_pdb: Path) -> Dict[str, object]:
    atom_count = 0
    residue_keys = set()
    kept_lines: List[str] = []
    try:
        lines = receptor_pdb.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        return {"status": "failed", "path": "", "atom_count": 0, "residue_count": 0, "note": str(exc)}

    for line in lines:
        record = line[:6].strip()
        if record == "ATOM":
            resname = line[17:20].strip().upper() if len(line) >= 20 else ""
            if resname and resname not in PDB_STANDARD_PROTEIN_RESNAMES:
                continue
            alt_loc = line[16:17].strip() if len(line) > 16 else ""
            if alt_loc and alt_loc not in {"A", "1"}:
                continue
            kept_lines.append(line)
            atom_count += 1
            residue_keys.add((line[21:22].strip(), line[22:26].strip(), line[26:27].strip(), resname))
        elif record == "TER" and kept_lines and not kept_lines[-1].startswith("TER"):
            kept_lines.append(line)

    if not kept_lines:
        return {"status": "failed", "path": "", "atom_count": 0, "residue_count": 0, "note": "No protein ATOM records found after receptor cleanup."}
    if not kept_lines[-1].startswith("END"):
        kept_lines.append("END")
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    output_pdb.write_text("\n".join(kept_lines) + "\n", encoding="utf-8")
    return {
        "status": "ready",
        "path": str(output_pdb),
        "atom_count": atom_count,
        "residue_count": len(residue_keys),
        "note": "Protein-only receptor PDB written for docking preparation; HETATM records were removed after pocket extraction.",
    }


def write_stage4_receptor_package(
    path: Path,
    project: Path,
    round_no: int,
    target_id: str,
    config: Dict[str, object],
    args: argparse.Namespace,
) -> Dict[str, object]:
    context = stage4_target_context(project, target_id, config)
    primary = stage4_primary_structure(context)
    binding_site = stage4_binding_site(context)
    requested_pdb_id = getattr(args, "pdb_id", "") or str(primary.get("pdb_id", "") or "")
    receptor_pdb = Path(str(args.receptor_pdb)).expanduser().resolve() if getattr(args, "receptor_pdb", None) else None
    local_receptor_path = ""
    receptor_pdbqt = ""
    fetch_status = "not_requested"
    if receptor_pdb:
        local_receptor_path = str(receptor_pdb)
        fetch_status = "local_present" if receptor_pdb.exists() else "local_missing"
    elif getattr(args, "fetch_receptor", False) and requested_pdb_id:
        receptor_pdb = project / "stage4" / "receptors" / f"{requested_pdb_id.upper()}.pdb"
        ok, message = fetch_pdb_file(requested_pdb_id, receptor_pdb)
        local_receptor_path = str(receptor_pdb) if ok else ""
        fetch_status = message if ok else f"fetch_failed: {message}"
        if not ok:
            fallback_receptor_path, fallback_receptor_pdbqt = find_stage4_project_receptor(project, requested_pdb_id)
            if fallback_receptor_path or fallback_receptor_pdbqt:
                local_receptor_path = fallback_receptor_path
                receptor_pdbqt = fallback_receptor_pdbqt
                fetch_status = f"{fetch_status}; local_project_receptor_found"
    elif requested_pdb_id:
        local_receptor_path, receptor_pdbqt = find_stage4_project_receptor(project, requested_pdb_id)
        if local_receptor_path or receptor_pdbqt:
            fetch_status = "local_project_receptor_found"

    preparation_status = "metadata_only"
    if local_receptor_path and Path(local_receptor_path).exists():
        preparation_status = "receptor_pdb_available"
    elif receptor_pdbqt and Path(receptor_pdbqt).exists():
        preparation_status = "receptor_pdbqt_available"
    elif requested_pdb_id:
        preparation_status = "pdb_id_known_receptor_file_not_prepared"

    if local_receptor_path and not receptor_pdbqt:
        candidate = Path(local_receptor_path).with_suffix(".pdbqt")
        if candidate.exists():
            receptor_pdbqt = str(candidate)

    co_crystal_ligand: Dict[str, object] = {}
    if local_receptor_path and Path(local_receptor_path).exists() and not stage4_has_binding_box(binding_site):
        co_crystal_ligand = extract_stage4_cocrystal_pocket(
            Path(local_receptor_path),
            str(binding_site.get("reference_ligand", "")),
        )
        if co_crystal_ligand.get("status") == "extracted_from_receptor":
            binding_site = dict(binding_site)
            original_site = dict(binding_site)
            binding_site.update(
                {
                    "status": "extracted_from_receptor",
                    "source": "co_crystal_ligand",
                    "center": co_crystal_ligand.get("center", []),
                    "size": co_crystal_ligand.get("size", []),
                    "detected_ligand": co_crystal_ligand.get("detected_ligand", {}),
                    "box_padding_angstrom": co_crystal_ligand.get("box_padding_angstrom", ""),
                    "box_min_size_angstrom": co_crystal_ligand.get("box_min_size_angstrom", ""),
                    "extraction_message": co_crystal_ligand.get("message", ""),
                    "catalog_binding_site": original_site,
                }
            )
        else:
            binding_site = dict(binding_site)
            binding_site.setdefault("status", "needs_curated_coordinates")
            binding_site["extraction_message"] = co_crystal_ligand.get("message", "")

    package = {
        "schema_version": "0.1",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "round": round_no,
        "target_id": target_id,
        "pdb_id": requested_pdb_id,
        "primary_structure": primary,
        "binding_site": binding_site,
        "co_crystal_ligand": co_crystal_ligand,
        "local_receptor_pdb": local_receptor_path,
        "local_receptor_pdbqt": receptor_pdbqt,
        "fetch_status": fetch_status,
        "preparation_status": preparation_status,
        "recommended_preparation": [
            "Remove crystallographic waters unless they mediate known binding-site interactions.",
            "Preserve catalytic metals/cofactors when the target requires them.",
            "Add hydrogens and assign protonation states before docking.",
            "Generate receptor PDBQT with Meeko or AutoDockTools before Vina.",
            "Use the co-crystallized ligand centroid or curated binding-site center for the docking box.",
        ],
        "boundary": "This package records receptor-readiness metadata and preparation instructions; it does not claim a prepared receptor unless a local receptor file exists.",
    }
    write_json(path, package)
    return package


def find_stage4_project_receptor(project: Path, pdb_id: str) -> Tuple[str, str]:
    if not pdb_id:
        return "", ""
    receptor_dir = project / "stage4" / "receptors"
    pdb_key = pdb_id.upper()
    pdb_candidates = [
        receptor_dir / f"{pdb_key}.pdb",
        receptor_dir / f"{pdb_key}_protein_only.pdb",
        receptor_dir / f"{pdb_key.lower()}.pdb",
    ]
    local_pdb = next((p.resolve() for p in pdb_candidates if p.exists() and p.is_file()), None)
    pdbqt_candidates: List[Path] = []
    if local_pdb:
        pdbqt_candidates.append(local_pdb.with_suffix(".pdbqt"))
        pdbqt_candidates.append(local_pdb.with_name(f"{local_pdb.stem}_protein_only_obabel.pdbqt"))
        pdbqt_candidates.append(local_pdb.with_name(f"{local_pdb.stem}_obabel.pdbqt"))
    pdbqt_candidates.extend(
        [
            receptor_dir / f"{pdb_key}_protein_only_obabel.pdbqt",
            receptor_dir / f"{pdb_key}_protein_only.pdbqt",
            receptor_dir / f"{pdb_key}.pdbqt",
        ]
    )
    local_pdbqt = next((p.resolve() for p in pdbqt_candidates if p.exists() and p.is_file()), None)
    return (str(local_pdb) if local_pdb else "", str(local_pdbqt) if local_pdbqt else "")


def rdkit_mol_and_fingerprint(smiles: str, libs: Dict[str, object]) -> Tuple[object, object]:
    chem = libs["Chem"]
    all_chem = libs["AllChem"]
    mol = chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    fingerprint = all_chem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return mol, fingerprint


def rdkit_stage4_descriptor(row: Dict[str, str], round_no: int, libs: Dict[str, object]) -> Dict[str, object]:
    chem = libs["Chem"]
    crippen = libs["Crippen"]
    descriptors = libs["Descriptors"]
    qed_lib = libs["QED"]
    rd_mol_descriptors = libs["rdMolDescriptors"]
    murcko = libs["MurckoScaffold"]

    smiles = row.get("smiles", "").strip()
    mol = chem.MolFromSmiles(smiles)
    base = {
        "round": row.get("round", round_no),
        "id": row.get("id", ""),
        "smiles": smiles,
        "parent": row.get("parent", ""),
        "source": row.get("source", ""),
        "descriptor_source": "rdkit",
    }
    if mol is None:
        out = dict(base)
        out.update(
            {
                "canonical_smiles": "",
                "valid": 0,
                "stage4_status": "invalid_smiles",
                "mol_formula": "",
                "exact_mw": "",
                "mw": "",
                "logp": "",
                "tpsa": "",
                "hbd": "",
                "hba": "",
                "rotatable_bonds": "",
                "heavy_atoms": "",
                "rings": "",
                "aromatic_rings": "",
                "formal_charge": "",
                "fraction_csp3": "",
                "qed": "",
                "lipinski_violations": 4,
                "murcko_scaffold": "",
                "risk_flags": "; ".join(risk_flags_for_smiles(smiles)),
                "error": "RDKit MolFromSmiles failed",
            }
        )
        return out

    canonical = chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    scaffold_mol = murcko.GetScaffoldForMol(mol)
    scaffold = chem.MolToSmiles(scaffold_mol, canonical=True) if scaffold_mol and scaffold_mol.GetNumAtoms() else ""
    desc = {
        "mw": round(float(descriptors.MolWt(mol)), 3),
        "logp": round(float(crippen.MolLogP(mol)), 3),
        "tpsa": round(float(rd_mol_descriptors.CalcTPSA(mol)), 3),
        "hbd": int(rd_mol_descriptors.CalcNumHBD(mol)),
        "hba": int(rd_mol_descriptors.CalcNumHBA(mol)),
    }
    lipinski = lipinski_violations(desc)
    out = dict(base)
    out.update(
        {
            "canonical_smiles": canonical,
            "valid": 1,
            "stage4_status": "ready_for_similarity_and_sdf",
            "mol_formula": rd_mol_descriptors.CalcMolFormula(mol),
            "exact_mw": round(float(descriptors.ExactMolWt(mol)), 5),
            "mw": desc["mw"],
            "logp": desc["logp"],
            "tpsa": desc["tpsa"],
            "hbd": desc["hbd"],
            "hba": desc["hba"],
            "rotatable_bonds": int(rd_mol_descriptors.CalcNumRotatableBonds(mol)),
            "heavy_atoms": int(mol.GetNumHeavyAtoms()),
            "rings": int(rd_mol_descriptors.CalcNumRings(mol)),
            "aromatic_rings": int(rd_mol_descriptors.CalcNumAromaticRings(mol)),
            "formal_charge": int(sum(atom.GetFormalCharge() for atom in mol.GetAtoms())),
            "fraction_csp3": round(float(rd_mol_descriptors.CalcFractionCSP3(mol)), 4),
            "qed": round(float(qed_lib.qed(mol)), 4),
            "lipinski_violations": lipinski,
            "murcko_scaffold": scaffold,
            "risk_flags": "; ".join(risk_flags_for_smiles(smiles)),
            "error": "",
        }
    )
    return out


def write_stage4_sdf(
    path: Path,
    descriptor_rows: Sequence[Dict[str, object]],
    libs: Dict[str, object],
    seed: int,
    max_conformers: int,
) -> Tuple[int, List[Dict[str, str]]]:
    chem = libs["Chem"]
    all_chem = libs["AllChem"]
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = chem.SDWriter(str(path))
    written = 0
    warnings: List[Dict[str, str]] = []
    try:
        for idx, row in enumerate(descriptor_rows, start=1):
            if str(row.get("valid", "")) != "1":
                continue
            smiles = str(row.get("canonical_smiles") or row.get("smiles", ""))
            mol = chem.MolFromSmiles(smiles)
            if mol is None:
                warnings.append({"id": str(row.get("id", "")), "warning": "MolFromSmiles failed during SDF export"})
                continue
            mol = chem.AddHs(mol)
            embedded = False
            if max_conformers > 0:
                status = all_chem.EmbedMolecule(mol, randomSeed=int(seed + idx), useRandomCoords=True)
                if status == 0:
                    embedded = True
                    try:
                        if all_chem.MMFFHasAllMoleculeParams(mol):
                            all_chem.MMFFOptimizeMolecule(mol, maxIters=250)
                        else:
                            all_chem.UFFOptimizeMolecule(mol, maxIters=250)
                    except Exception as exc:
                        warnings.append({"id": str(row.get("id", "")), "warning": f"force-field optimization failed: {exc}"})
                else:
                    warnings.append({"id": str(row.get("id", "")), "warning": "3D embedding failed; writing 2D/graph SDF"})
            mol.SetProp("_Name", str(row.get("id", "")) or f"ligand_{idx:05d}")
            for prop in ["id", "canonical_smiles", "qed", "mw", "logp", "tpsa", "hbd", "hba", "rotatable_bonds"]:
                mol.SetProp(prop, str(row.get(prop, "")))
            mol.SetProp("stage4_coordinate_status", "3d_embedded" if embedded else "not_3d_embedded")
            writer.write(mol)
            written += 1
    finally:
        writer.close()
    return written, warnings


def stage4_similarity_rows(
    descriptor_rows: Sequence[Dict[str, object]],
    controls: Sequence[Dict[str, object]],
    libs: Dict[str, object],
) -> Tuple[List[Dict[str, object]], Dict[str, Dict[str, object]]]:
    data_structs = libs["DataStructs"]
    control_fps: List[Tuple[Dict[str, object], object]] = []
    for control in controls:
        _, fp = rdkit_mol_and_fingerprint(str(control.get("smiles", "")), libs)
        if fp is not None:
            control_fps.append((dict(control), fp))

    rows: List[Dict[str, object]] = []
    best_by_candidate: Dict[str, Dict[str, object]] = {}
    for row in descriptor_rows:
        if str(row.get("valid", "")) != "1":
            continue
        _, fp = rdkit_mol_and_fingerprint(str(row.get("canonical_smiles") or row.get("smiles", "")), libs)
        if fp is None:
            continue
        candidate_id = str(row.get("id", ""))
        best: Optional[Dict[str, object]] = None
        for control, control_fp in control_fps:
            similarity = round(float(data_structs.TanimotoSimilarity(fp, control_fp)), 4)
            sim_row = {
                "candidate_id": candidate_id,
                "candidate_smiles": row.get("canonical_smiles", ""),
                "control_drug": control.get("drug", ""),
                "control_role": control.get("role", ""),
                "control_smiles": control.get("smiles", ""),
                "tanimoto": similarity,
                "fingerprint": "rdkit_morgan_2048_r2",
                "control_source": control.get("source", ""),
                "pubchem_cid": control.get("pubchem_cid", ""),
            }
            rows.append(sim_row)
            if best is None or similarity > float(best.get("tanimoto", 0.0)):
                best = sim_row
        if best:
            best_by_candidate[candidate_id] = best

    rows.sort(key=lambda item: (str(item.get("candidate_id", "")), -float(item.get("tanimoto", 0.0)), str(item.get("control_drug", ""))))
    return rows, best_by_candidate


def stage4_diversity_rows(
    descriptor_rows: Sequence[Dict[str, object]],
    best_by_candidate: Dict[str, Dict[str, object]],
    libs: Dict[str, object],
    top: int,
) -> List[Dict[str, object]]:
    data_structs = libs["DataStructs"]
    valid_rows = [row for row in descriptor_rows if str(row.get("valid", "")) == "1"]
    valid_rows.sort(key=lambda row: (float(row.get("qed") or 0.0), -float(row.get("lipinski_violations") or 0.0)), reverse=True)
    selected: List[Tuple[Dict[str, object], object]] = []
    output: List[Dict[str, object]] = []
    for row in valid_rows:
        if len(output) >= top:
            break
        _, fp = rdkit_mol_and_fingerprint(str(row.get("canonical_smiles") or row.get("smiles", "")), libs)
        if fp is None:
            continue
        max_selected = 0.0
        if selected:
            max_selected = max(float(data_structs.TanimotoSimilarity(fp, selected_fp)) for _, selected_fp in selected)
        if selected and max_selected > 0.82:
            continue
        best_control = best_by_candidate.get(str(row.get("id", "")), {})
        output.append(
            {
                "selection_rank": len(output) + 1,
                "id": row.get("id", ""),
                "canonical_smiles": row.get("canonical_smiles", ""),
                "qed": row.get("qed", ""),
                "max_similarity_to_selected": round(max_selected, 4),
                "max_similarity_to_controls": best_control.get("tanimoto", ""),
                "nearest_control": best_control.get("control_drug", ""),
                "murcko_scaffold": row.get("murcko_scaffold", ""),
            }
        )
        selected.append((row, fp))
    return output


def stage4_descriptor_for_smiles(
    panel_type: str,
    item_id: str,
    smiles: str,
    role: str,
    target_id: str,
    source: str,
    libs: Dict[str, object],
) -> Dict[str, object]:
    row = {
        "round": "",
        "id": item_id,
        "smiles": smiles,
        "parent": "",
        "source": source,
    }
    desc = rdkit_stage4_descriptor(row, 0, libs)
    desc["panel_type"] = panel_type
    desc["role"] = role
    desc["target_id"] = target_id
    return desc


def build_stage4_decoys(count: int, libs: Dict[str, object]) -> List[Dict[str, object]]:
    decoys: List[Dict[str, object]] = []
    for item_id, smiles, note in DEFAULT_DECOY_SMILES[: max(0, count)]:
        desc = stage4_descriptor_for_smiles("decoy", item_id, smiles, "decoy", "", note, libs)
        decoys.append(desc)
    return decoys


def nearest_control_summary(
    item_smiles: str,
    controls: Sequence[Dict[str, object]],
    libs: Dict[str, object],
) -> Tuple[object, object]:
    data_structs = libs["DataStructs"]
    _, fp = rdkit_mol_and_fingerprint(item_smiles, libs)
    if fp is None:
        return "", ""
    best_name = ""
    best_value = 0.0
    for control in controls:
        _, control_fp = rdkit_mol_and_fingerprint(str(control.get("smiles", "")), libs)
        if control_fp is None:
            continue
        value = float(data_structs.TanimotoSimilarity(fp, control_fp))
        if value > best_value:
            best_value = value
            best_name = str(control.get("drug", ""))
    return round(best_value, 4), best_name


def stage4_benchmark_panel_rows(
    descriptor_rows: Sequence[Dict[str, object]],
    controls: Sequence[Dict[str, object]],
    decoy_rows: Sequence[Dict[str, object]],
    best_by_candidate: Dict[str, Dict[str, object]],
    target_id: str,
    libs: Dict[str, object],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for row in descriptor_rows:
        if str(row.get("valid", "")) != "1":
            continue
        best = best_by_candidate.get(str(row.get("id", "")), {})
        rows.append(
            {
                "panel_type": "candidate",
                "id": row.get("id", ""),
                "smiles": row.get("smiles", ""),
                "canonical_smiles": row.get("canonical_smiles", ""),
                "role": "candidate",
                "target_id": target_id,
                "qed": row.get("qed", ""),
                "mw": row.get("mw", ""),
                "logp": row.get("logp", ""),
                "tpsa": row.get("tpsa", ""),
                "hbd": row.get("hbd", ""),
                "hba": row.get("hba", ""),
                "heavy_atoms": row.get("heavy_atoms", ""),
                "murcko_scaffold": row.get("murcko_scaffold", ""),
                "max_similarity_to_controls": best.get("tanimoto", ""),
                "nearest_control": best.get("control_drug", ""),
                "source": row.get("source", ""),
            }
        )

    for control in controls:
        desc = stage4_descriptor_for_smiles(
            str(control.get("role", "control")) or "control",
            str(control.get("drug", "")),
            str(control.get("smiles", "")),
            str(control.get("role", "")),
            target_id,
            str(control.get("source", "")),
            libs,
        )
        if str(desc.get("valid", "")) != "1":
            continue
        rows.append(
            {
                "panel_type": str(control.get("role", "control")) or "control",
                "id": control.get("drug", ""),
                "smiles": control.get("smiles", ""),
                "canonical_smiles": desc.get("canonical_smiles", ""),
                "role": control.get("role", ""),
                "target_id": target_id,
                "qed": desc.get("qed", ""),
                "mw": desc.get("mw", ""),
                "logp": desc.get("logp", ""),
                "tpsa": desc.get("tpsa", ""),
                "hbd": desc.get("hbd", ""),
                "hba": desc.get("hba", ""),
                "heavy_atoms": desc.get("heavy_atoms", ""),
                "murcko_scaffold": desc.get("murcko_scaffold", ""),
                "max_similarity_to_controls": 1.0,
                "nearest_control": control.get("drug", ""),
                "source": control.get("source", ""),
            }
        )

    for decoy in decoy_rows:
        if str(decoy.get("valid", "")) != "1":
            continue
        max_sim, nearest = nearest_control_summary(str(decoy.get("canonical_smiles") or decoy.get("smiles", "")), controls, libs)
        rows.append(
            {
                "panel_type": "decoy",
                "id": decoy.get("id", ""),
                "smiles": decoy.get("smiles", ""),
                "canonical_smiles": decoy.get("canonical_smiles", ""),
                "role": "decoy",
                "target_id": target_id,
                "qed": decoy.get("qed", ""),
                "mw": decoy.get("mw", ""),
                "logp": decoy.get("logp", ""),
                "tpsa": decoy.get("tpsa", ""),
                "hbd": decoy.get("hbd", ""),
                "hba": decoy.get("hba", ""),
                "heavy_atoms": decoy.get("heavy_atoms", ""),
                "murcko_scaffold": decoy.get("murcko_scaffold", ""),
                "max_similarity_to_controls": max_sim,
                "nearest_control": nearest,
                "source": decoy.get("source", ""),
            }
        )
    return rows


def write_stage4_decoys(path: Path, decoy_rows: Sequence[Dict[str, object]]) -> None:
    write_csv(path, decoy_rows, STAGE4_DESCRIPTOR_FIELDS)


def render_stage4_2d_images(
    output_dir: Path,
    descriptor_rows: Sequence[Dict[str, object]],
    controls: Sequence[Dict[str, object]],
    decoy_rows: Sequence[Dict[str, object]],
    libs: Dict[str, object],
    top: int,
    enabled: bool,
) -> List[Dict[str, str]]:
    if not enabled:
        return []
    try:
        from rdkit.Chem import Draw  # type: ignore
    except Exception:
        return []
    chem = libs["Chem"]
    output_dir.mkdir(parents=True, exist_ok=True)
    items: List[Tuple[str, str, str]] = []
    for row in descriptor_rows:
        if str(row.get("valid", "")) == "1":
            items.append(("candidate", str(row.get("id", "")), str(row.get("canonical_smiles") or row.get("smiles", ""))))
        if len([item for item in items if item[0] == "candidate"]) >= top:
            break
    for control in controls[:4]:
        items.append(("control", str(control.get("drug", "")), str(control.get("smiles", ""))))
    for decoy in decoy_rows[:2]:
        items.append(("decoy", str(decoy.get("id", "")), str(decoy.get("canonical_smiles") or decoy.get("smiles", ""))))

    rendered: List[Dict[str, str]] = []
    for panel_type, item_id, smiles in items:
        mol = chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", item_id).strip("_") or "molecule"
        path = output_dir / f"{panel_type}_{safe_id}.png"
        Draw.MolToFile(mol, str(path), size=(360, 260), legend=f"{panel_type}: {item_id}")
        rendered.append({"panel_type": panel_type, "id": item_id, "path": str(path)})
    return rendered


def stage4_docking_inputs(
    benchmark_rows: Sequence[Dict[str, object]],
    receptor_package: Dict[str, object],
    sdf_path: Path,
    ligand_records: Optional[Dict[str, Dict[str, str]]] = None,
    dockable_panel_types: Optional[Sequence[str]] = None,
) -> List[Dict[str, object]]:
    receptor_pdb = str(receptor_package.get("local_receptor_pdb", ""))
    receptor_pdbqt = str(receptor_package.get("local_receptor_pdbqt", ""))
    ligand_records = ligand_records or {}
    panel_types = set(dockable_panel_types) if dockable_panel_types is not None else None
    rows = []
    for row in benchmark_rows:
        panel_type = str(row.get("panel_type", ""))
        if panel_types is None and panel_type == "decoy":
            continue
        if panel_types is not None and panel_type not in panel_types:
            continue
        item_id = str(row.get("id", ""))
        prepared = ligand_records.get(item_id, {})
        ligand_sdf = prepared.get("ligand_sdf") or (str(sdf_path) if panel_type == "candidate" else "")
        ligand_pdbqt = prepared.get("ligand_pdbqt", "")
        status = prepared.get("status") or ("ready_for_ligand_preparation" if Path(sdf_path).exists() else "missing_sdf")
        note = prepared.get("note") or "Candidate ligands share the round SDF; controls require separate SDF/PDBQT export before docking."
        rows.append(
            {
                "panel_type": panel_type,
                "id": item_id,
                "canonical_smiles": row.get("canonical_smiles", ""),
                "ligand_sdf": ligand_sdf,
                "ligand_pdbqt": ligand_pdbqt,
                "receptor_pdb": receptor_pdb,
                "receptor_pdbqt": receptor_pdbqt,
                "status": status,
                "note": note,
            }
        )
    return rows


def write_stage4_benchmark_sdf(
    path: Path,
    benchmark_rows: Sequence[Dict[str, object]],
    libs: Dict[str, object],
    seed: int,
) -> int:
    chem = libs["Chem"]
    all_chem = libs["AllChem"]
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = chem.SDWriter(str(path))
    written = 0
    try:
        for idx, row in enumerate(benchmark_rows, start=1):
            smiles = str(row.get("canonical_smiles") or row.get("smiles", ""))
            mol = chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mol = chem.AddHs(mol)
            status = all_chem.EmbedMolecule(mol, randomSeed=int(seed + 10_000 + idx), useRandomCoords=True)
            if status == 0:
                try:
                    if all_chem.MMFFHasAllMoleculeParams(mol):
                        all_chem.MMFFOptimizeMolecule(mol, maxIters=200)
                    else:
                        all_chem.UFFOptimizeMolecule(mol, maxIters=200)
                except Exception:
                    pass
            item_id = str(row.get("id", "")) or f"panel_{idx:05d}"
            mol.SetProp("_Name", item_id)
            for prop in ["panel_type", "id", "role", "target_id", "canonical_smiles", "qed", "mw", "source"]:
                mol.SetProp(prop, str(row.get(prop, "")))
            writer.write(mol)
            written += 1
    finally:
        writer.close()
    return written


def stage4_box_from_receptor_package(receptor_package: Dict[str, object]) -> Dict[str, object]:
    site = receptor_package.get("binding_site", {})
    if not isinstance(site, dict):
        site = {}
    center = site.get("center") or site.get("box_center")
    size = site.get("size") or site.get("box_size")
    if isinstance(center, str):
        center = parse_triplet(center)
    if isinstance(size, str):
        size = parse_triplet(size)
    if isinstance(center, list) and isinstance(size, list) and len(center) == 3 and len(size) == 3:
        try:
            return {"status": "ready", "center": [float(x) for x in center], "size": [float(x) for x in size], "source": "binding_site"}
        except Exception:
            pass
    return {
        "status": "missing",
        "center": [],
        "size": [],
        "source": site.get("source", ""),
        "message": "Docking box center/size is missing; use curated co-crystal coordinates before real docking.",
    }


def stage4_prepare_receptor_pdbqt(
    receptor_package: Dict[str, object],
    capabilities: Dict[str, object],
    docking_box: Dict[str, object],
    timeout: int,
) -> Dict[str, object]:
    existing = str(receptor_package.get("local_receptor_pdbqt", ""))
    if existing and Path(existing).exists():
        return {"status": "ready", "receptor_pdbqt": existing, "tool": "existing", "note": "Existing receptor PDBQT found."}
    receptor_pdb = str(receptor_package.get("local_receptor_pdb", ""))
    if not receptor_pdb or not Path(receptor_pdb).exists():
        return {"status": "missing_receptor_pdb", "receptor_pdbqt": "", "tool": "", "note": "No local receptor PDB available."}
    executables = capabilities.get("executables", {})
    exe = executables.get("mk_prepare_receptor.py", {}) if isinstance(executables, dict) else {}
    path = str(exe.get("path", "")) if isinstance(exe, dict) else ""
    obabel_exe = executables.get("obabel", {}) if isinstance(executables, dict) else {}
    obabel_path = str(obabel_exe.get("path", "")) if isinstance(obabel_exe, dict) else ""
    if not path and not obabel_path:
        return {"status": "missing_receptor_preparation_tool", "receptor_pdbqt": "", "tool": "", "note": "Neither mk_prepare_receptor.py nor obabel is available."}
    if docking_box.get("status") != "ready":
        return {"status": "missing_docking_box", "receptor_pdbqt": "", "tool": path, "note": str(docking_box.get("message", "missing docking box"))}
    center = docking_box.get("center", [])
    size = docking_box.get("size", [])

    def run_prepare(input_pdb: Path, output_base: Path) -> Tuple[List[str], Optional[subprocess.CompletedProcess], Path, str]:
        output_pdbqt = Path(str(output_base) + ".pdbqt")
        command = [
            path,
            "--read_pdb",
            str(input_pdb),
            "-o",
            str(output_base),
            "-p",
            str(output_pdbqt),
            "--box_center",
            str(center[0]),
            str(center[1]),
            str(center[2]),
            "--box_size",
            str(size[0]),
            str(size[1]),
            str(size[2]),
            "-a",
        ]
        try:
            proc = subprocess.run(command, check=False, capture_output=True, text=True, timeout=max(1, timeout))
            return command, proc, output_pdbqt, ""
        except Exception as exc:
            return command, None, output_pdbqt, str(exc)

    def run_obabel_prepare(input_pdb: Path, output_pdbqt: Path) -> Tuple[List[str], Optional[subprocess.CompletedProcess], str]:
        command = [obabel_path, str(input_pdb), "-O", str(output_pdbqt), "-xr", "-h"]
        try:
            proc = subprocess.run(command, check=False, capture_output=True, text=True, timeout=max(1, timeout))
            return command, proc, ""
        except Exception as exc:
            return command, None, str(exc)

    raw_input = Path(receptor_pdb)
    out_base = raw_input.with_suffix("")
    command: List[str] = []
    proc: Optional[subprocess.CompletedProcess] = None
    out_pdbqt = Path(str(out_base) + ".pdbqt")
    exception_note = ""
    if path:
        command, proc, out_pdbqt, exception_note = run_prepare(raw_input, out_base)
        if proc is not None and proc.returncode == 0 and out_pdbqt.exists():
            receptor_package["local_receptor_pdbqt"] = str(out_pdbqt)
            receptor_package["preparation_status"] = "receptor_pdbqt_available"
            return {"status": "ready", "receptor_pdbqt": str(out_pdbqt), "tool": path, "note": "Prepared receptor PDBQT with Meeko.", "command": command}

    raw_note = exception_note or ((proc.stderr or proc.stdout or "mk_prepare_receptor.py failed") if proc is not None else "mk_prepare_receptor.py failed")
    cleaned_pdb = raw_input.with_name(f"{raw_input.stem}_protein_only.pdb")
    cleanup = clean_stage4_receptor_pdb_for_docking(raw_input, cleaned_pdb)
    if cleanup.get("status") == "ready":
        cleaned_base = cleaned_pdb.with_suffix("")
        fallback_command: List[str] = []
        fallback_proc: Optional[subprocess.CompletedProcess] = None
        fallback_pdbqt = Path(str(cleaned_base) + ".pdbqt")
        fallback_exception = ""
        fallback_note = ""
        if path:
            fallback_command, fallback_proc, fallback_pdbqt, fallback_exception = run_prepare(cleaned_pdb, cleaned_base)
            if fallback_proc is not None and fallback_proc.returncode == 0 and fallback_pdbqt.exists():
                receptor_package["local_receptor_pdbqt"] = str(fallback_pdbqt)
                receptor_package["prepared_receptor_pdb"] = str(cleaned_pdb)
                receptor_package["preparation_status"] = "receptor_pdbqt_available"
                receptor_package["preparation_note"] = "Prepared receptor PDBQT from protein-only cleaned PDB after original receptor preparation failed."
                return {
                    "status": "ready",
                    "receptor_pdbqt": str(fallback_pdbqt),
                    "tool": path,
                    "note": "Prepared receptor PDBQT with Meeko from protein-only cleaned PDB.",
                    "command": fallback_command,
                    "original_command": command,
                    "original_note": raw_note[-800:],
                    "fallback": "protein_only_clean_pdb",
                    "cleaned_receptor_pdb": str(cleaned_pdb),
                    "cleanup": cleanup,
                }
            fallback_note = fallback_exception or ((fallback_proc.stderr or fallback_proc.stdout or "mk_prepare_receptor.py failed on cleaned receptor") if fallback_proc is not None else "mk_prepare_receptor.py failed on cleaned receptor")
        if obabel_path:
            obabel_pdbqt = cleaned_pdb.with_name(f"{cleaned_pdb.stem}_obabel.pdbqt")
            obabel_command, obabel_proc, obabel_exception = run_obabel_prepare(cleaned_pdb, obabel_pdbqt)
            if obabel_proc is not None and obabel_proc.returncode == 0 and obabel_pdbqt.exists():
                receptor_package["local_receptor_pdbqt"] = str(obabel_pdbqt)
                receptor_package["prepared_receptor_pdb"] = str(cleaned_pdb)
                receptor_package["preparation_status"] = "receptor_pdbqt_available"
                receptor_package["preparation_note"] = "Prepared receptor PDBQT with OpenBabel from protein-only cleaned PDB after Meeko receptor preparation failed."
                return {
                    "status": "ready",
                    "receptor_pdbqt": str(obabel_pdbqt),
                    "tool": "openbabel",
                    "note": "Prepared receptor PDBQT with OpenBabel from protein-only cleaned PDB.",
                    "command": obabel_command,
                    "original_command": command,
                    "original_note": raw_note[-800:],
                    "meeko_cleaned_command": fallback_command,
                    "meeko_cleaned_note": fallback_note[-800:],
                    "fallback": "openbabel_protein_only_pdbqt",
                    "cleaned_receptor_pdb": str(cleaned_pdb),
                    "cleanup": cleanup,
                }
            obabel_note = obabel_exception or ((obabel_proc.stderr or obabel_proc.stdout or "OpenBabel receptor PDBQT conversion failed") if obabel_proc is not None else "OpenBabel receptor PDBQT conversion failed")
            fallback_note = (fallback_note + "\nOpenBabel: " + obabel_note).strip()
        return {
            "status": "prepare_failed",
            "receptor_pdbqt": "",
            "tool": path,
            "note": fallback_note[-800:],
            "command": fallback_command,
            "original_command": command,
            "original_note": raw_note[-800:],
            "fallback": "protein_only_clean_pdb",
            "cleaned_receptor_pdb": str(cleaned_pdb),
            "cleanup": cleanup,
        }
    return {
        "status": "prepare_failed",
        "receptor_pdbqt": "",
        "tool": path,
        "note": raw_note[-800:],
        "command": command,
        "cleanup": cleanup,
    }


def stage4_write_single_ligand_sdf(
    path: Path,
    row: Dict[str, object],
    libs: Dict[str, object],
    seed: int,
    idx: int,
) -> bool:
    chem = libs["Chem"]
    all_chem = libs["AllChem"]
    smiles = str(row.get("canonical_smiles") or row.get("smiles", ""))
    mol = chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    mol = chem.AddHs(mol)
    status = all_chem.EmbedMolecule(mol, randomSeed=int(seed + 20_000 + idx), useRandomCoords=True)
    if status == 0:
        try:
            if all_chem.MMFFHasAllMoleculeParams(mol):
                all_chem.MMFFOptimizeMolecule(mol, maxIters=200)
            else:
                all_chem.UFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
    mol.SetProp("_Name", str(row.get("id", "")) or f"ligand_{idx:05d}")
    mol.SetProp("id", str(row.get("id", "")))
    mol.SetProp("canonical_smiles", smiles)
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = chem.SDWriter(str(path))
    try:
        writer.write(mol)
    finally:
        writer.close()
    return True


def stage4_prepare_ligand_records(
    benchmark_rows: Sequence[Dict[str, object]],
    libs: Dict[str, object],
    capabilities: Dict[str, object],
    output_dir: Path,
    seed: int,
    timeout: int,
    enabled: bool,
    dockable_panel_types: Optional[Sequence[str]] = None,
) -> Tuple[Dict[str, Dict[str, str]], List[Dict[str, object]]]:
    records: Dict[str, Dict[str, str]] = {}
    attempts: List[Dict[str, object]] = []
    if not enabled:
        return records, attempts
    executables = capabilities.get("executables", {})
    exe = executables.get("mk_prepare_ligand.py", {}) if isinstance(executables, dict) else {}
    meeko_tool = str(exe.get("path", "")) if isinstance(exe, dict) else ""
    obabel_exe = executables.get("obabel", {}) if isinstance(executables, dict) else {}
    obabel_tool = str(obabel_exe.get("path", "")) if isinstance(obabel_exe, dict) else ""
    output_dir.mkdir(parents=True, exist_ok=True)
    panel_types = set(dockable_panel_types or ["candidate"])
    for idx, row in enumerate(benchmark_rows, start=1):
        if str(row.get("panel_type", "")) not in panel_types:
            continue
        item_id = str(row.get("id", "")) or f"ligand_{idx:05d}"
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", item_id).strip("_") or f"ligand_{idx:05d}"
        ligand_sdf = output_dir / f"{safe_id}.sdf"
        ligand_pdbqt = output_dir / f"{safe_id}.pdbqt"
        wrote = stage4_write_single_ligand_sdf(ligand_sdf, row, libs, seed, idx)
        if not wrote:
            records[item_id] = {"ligand_sdf": "", "ligand_pdbqt": "", "status": "invalid_ligand", "note": "Could not write ligand SDF."}
            continue
        prep_errors: List[str] = []
        if meeko_tool:
            command = [meeko_tool, "-i", str(ligand_sdf), "-o", str(ligand_pdbqt)]
            try:
                proc = subprocess.run(command, check=False, capture_output=True, text=True, timeout=max(1, timeout))
                ok = proc.returncode == 0 and ligand_pdbqt.exists()
                status = "ready_for_docking" if ok else "ligand_prepare_failed"
                note = "Prepared ligand PDBQT with Meeko." if ok else (proc.stderr or proc.stdout or "mk_prepare_ligand.py failed")[-500:]
                attempts.append({"id": item_id, "status": status, "returncode": proc.returncode, "command": command, "tool": "meeko"})
                if ok:
                    records[item_id] = {
                        "ligand_sdf": str(ligand_sdf),
                        "ligand_pdbqt": str(ligand_pdbqt),
                        "status": status,
                        "note": note,
                        "tool": "meeko",
                    }
                    continue
                prep_errors.append(f"meeko: {note}")
            except Exception as exc:
                attempts.append({"id": item_id, "status": "ligand_prepare_failed", "error": str(exc), "command": command, "tool": "meeko"})
                prep_errors.append(f"meeko: {exc}")
        if obabel_tool:
            command = [obabel_tool, str(ligand_sdf), "-O", str(ligand_pdbqt), "-h"]
            try:
                proc = subprocess.run(command, check=False, capture_output=True, text=True, timeout=max(1, timeout))
                ok = proc.returncode == 0 and ligand_pdbqt.exists()
                status = "ready_for_docking" if ok else "ligand_prepare_failed"
                note = "Prepared ligand PDBQT with OpenBabel." if ok else (proc.stderr or proc.stdout or "obabel failed")[-500:]
                attempts.append({"id": item_id, "status": status, "returncode": proc.returncode, "command": command, "tool": "openbabel"})
                records[item_id] = {
                    "ligand_sdf": str(ligand_sdf),
                    "ligand_pdbqt": str(ligand_pdbqt) if ok else "",
                    "status": status,
                    "note": note if ok else "; ".join(prep_errors + [f"openbabel: {note}"]),
                    "tool": "openbabel",
                }
                continue
            except Exception as exc:
                attempts.append({"id": item_id, "status": "ligand_prepare_failed", "error": str(exc), "command": command, "tool": "openbabel"})
                prep_errors.append(f"openbabel: {exc}")
        if not meeko_tool and not obabel_tool:
            records[item_id] = {"ligand_sdf": str(ligand_sdf), "ligand_pdbqt": "", "status": "missing_ligand_preparation_tool", "note": "Neither mk_prepare_ligand.py nor obabel is available."}
            continue
        records[item_id] = {"ligand_sdf": str(ligand_sdf), "ligand_pdbqt": "", "status": "ligand_prepare_failed", "note": "; ".join(prep_errors) or "Ligand preparation failed."}
    return records, attempts


def write_stage4_docking_plan(
    path: Path,
    backend: str,
    run_docking: bool,
    capabilities: Dict[str, object],
    receptor_package: Dict[str, object],
    docking_inputs_path: Path,
    docking_scores_path: Path,
    docking_box: Optional[Dict[str, object]] = None,
    receptor_preparation: Optional[Dict[str, object]] = None,
    ligand_preparation_count: int = 0,
) -> Dict[str, object]:
    executables = capabilities.get("executables", {})
    if not isinstance(executables, dict):
        executables = {}
    available = []
    for name in ["vina", "gnina"]:
        status = executables.get(name, {})
        if isinstance(status, dict) and status.get("status") == "found":
            available.append(name)
    selected = ""
    if backend in {"vina", "gnina"}:
        selected = backend if backend in available else ""
    elif available:
        selected = available[0]

    docking_box = docking_box or {}
    receptor_preparation = receptor_preparation or {}
    receptor_ready = bool(receptor_package.get("local_receptor_pdbqt") or receptor_package.get("local_receptor_pdb"))
    box_ready = docking_box.get("status") == "ready"
    if run_docking and selected and receptor_ready and box_ready:
        status = "planned"
    elif run_docking and not selected:
        status = "not_available"
    elif run_docking and not receptor_ready:
        status = "blocked_missing_receptor"
    elif run_docking and not box_ready:
        status = "blocked_missing_box"
    else:
        status = "skipped"

    commands = []
    if selected == "vina":
        commands.append(
            "vina --receptor <receptor.pdbqt> --ligand <ligand.pdbqt> --center_x <x> --center_y <y> --center_z <z> --size_x <sx> --size_y <sy> --size_z <sz> --out <pose.pdbqt>"
        )
    elif selected == "gnina":
        commands.append(
            "gnina --receptor <receptor.pdb> --ligand <ligand.sdf> --center_x <x> --center_y <y> --center_z <z> --size_x <sx> --size_y <sy> --size_z <sz> --out <poses.sdf>"
        )
    commands.extend(
        [
            "Convert ligand SDF to PDBQT with Meeko/OpenBabel before Vina.",
            "After docking, write id,smiles,docking_score,pose_pass to the external scores CSV and run score --external-scores.",
        ]
    )
    plan = {
        "schema_version": "0.1",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "requested_backend": backend,
        "selected_backend": selected,
        "status": status,
        "run_docking_requested": bool(run_docking),
        "required_tools": {
            "ligand_preparation": ["meeko", "obabel"],
            "docking": ["vina", "gnina"],
            "pose_quality": ["posebusters"],
        },
        "tool_status": {
            "vina": executables.get("vina", {}),
            "gnina": executables.get("gnina", {}),
            "obabel": executables.get("obabel", {}),
            "mk_prepare_ligand.py": executables.get("mk_prepare_ligand.py", {}),
            "mk_prepare_receptor.py": executables.get("mk_prepare_receptor.py", {}),
            "posebusters_cli": executables.get("bust", {}),
        },
        "receptor_status": receptor_package.get("preparation_status", ""),
        "docking_box": docking_box,
        "receptor_preparation": receptor_preparation,
        "ligand_preparation_count": ligand_preparation_count,
        "docking_inputs": str(docking_inputs_path),
        "expected_scores_csv": str(docking_scores_path),
        "commands": commands,
        "message": "Docking is not executed unless --run-docking is set and required executables/receptor files are available.",
    }
    write_json(path, plan)
    return plan


def stage4_run_posebusters_if_possible(
    item_id: str,
    pose_path: Path,
    receptor_pdb: str,
    run_dir: Path,
    executables: Dict[str, object],
    timeout: int,
) -> Dict[str, object]:
    bust_status = executables.get("posebusters_cli", {})
    bust_path = str(bust_status.get("path", "")) if isinstance(bust_status, dict) else ""
    if not bust_path:
        return {"status": "not_available", "note": "PoseBusters CLI is not available."}
    if not pose_path.exists():
        return {"status": "missing_pose", "note": f"Pose file not found: {pose_path}"}

    mol_pred = pose_path
    conversion_command: List[str] = []
    if pose_path.suffix.lower() != ".sdf":
        obabel_status = executables.get("obabel", {})
        obabel_path = str(obabel_status.get("path", "")) if isinstance(obabel_status, dict) else ""
        if not obabel_path:
            return {"status": "missing_obabel", "note": "PoseBusters needs SDF pose conversion, but obabel is not available."}
        mol_pred = run_dir / f"{item_id}_pose.sdf"
        conversion_command = [obabel_path, str(pose_path), "-O", str(mol_pred)]
        try:
            convert = subprocess.run(conversion_command, check=False, capture_output=True, text=True, timeout=max(1, timeout))
        except Exception as exc:
            return {"status": "conversion_failed", "note": str(exc), "command": conversion_command}
        if convert.returncode != 0 or not mol_pred.exists():
            return {
                "status": "conversion_failed",
                "note": (convert.stderr or convert.stdout or "obabel pose conversion failed")[-500:],
                "command": conversion_command,
            }

    report_path = run_dir / f"{item_id}_posebusters.csv"
    command = [bust_path, "--outfmt", "csv", "--output", str(report_path), str(mol_pred)]
    if receptor_pdb and Path(receptor_pdb).exists():
        command.extend(["-p", receptor_pdb])
    try:
        proc = subprocess.run(command, check=False, capture_output=True, text=True, timeout=max(1, timeout))
    except Exception as exc:
        return {"status": "failed", "note": str(exc), "command": command, "conversion_command": conversion_command, "report": str(report_path)}
    report_text = report_path.read_text(encoding="utf-8", errors="replace") if report_path.exists() else (proc.stdout or "")
    lowered = report_text.lower()
    passed = proc.returncode == 0 and "false" not in lowered and "failed" not in lowered
    return {
        "status": "passed" if passed else "failed",
        "note": (proc.stderr or proc.stdout or "")[-500:],
        "command": command,
        "conversion_command": conversion_command,
        "report": str(report_path),
        "mol_pred": str(mol_pred),
        "returncode": proc.returncode,
    }


def run_stage4_docking_if_possible(
    docking_plan: Dict[str, object],
    docking_input_rows: Sequence[Dict[str, object]],
    descriptor_rows: Sequence[Dict[str, object]],
    output_csv: Path,
    timeout: int,
    dockable_panel_types: Optional[Sequence[str]] = None,
    output_fields: Optional[Sequence[str]] = None,
    run_dir_name: str = "docking_runs",
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    if docking_plan.get("status") != "planned":
        return docking_plan, []
    backend = str(docking_plan.get("selected_backend", ""))
    if backend not in {"vina", "gnina"}:
        docking_plan["status"] = "not_available"
        docking_plan["message"] = "No supported docking backend selected."
        return docking_plan, []

    descriptor_by_id = {str(row.get("id", "")): row for row in descriptor_rows}
    panel_types = set(dockable_panel_types or ["candidate"])
    executables = docking_plan.get("tool_status", {})
    if not isinstance(executables, dict):
        executables = {}
    backend_status = executables.get(backend, {})
    backend_exe = str(backend_status.get("path", "")) if isinstance(backend_status, dict) else ""
    if not backend_exe:
        backend_exe = backend
    docking_box = docking_plan.get("docking_box", {})
    if not isinstance(docking_box, dict):
        docking_box = {}
    center = docking_box.get("center") or []
    size = docking_box.get("size") or []
    if len(center) != 3 or len(size) != 3:
        docking_plan["status"] = "blocked_missing_box"
        docking_plan["message"] = "Docking box center/size is missing; docking command was not run."
        return docking_plan, []
    results: List[Dict[str, object]] = []
    run_dir = output_csv.parent / run_dir_name
    run_dir.mkdir(parents=True, exist_ok=True)

    for item in docking_input_rows:
        panel_type = str(item.get("panel_type", ""))
        if panel_type not in panel_types:
            continue
        ligand_pdbqt = str(item.get("ligand_pdbqt", ""))
        receptor_pdbqt = str(item.get("receptor_pdbqt", ""))
        receptor_pdb = str(item.get("receptor_pdb", ""))
        item_id = str(item.get("id", ""))
        desc = descriptor_by_id.get(item_id, {})
        smiles = str(desc.get("canonical_smiles") or desc.get("smiles") or item.get("canonical_smiles", ""))
        log_path = run_dir / f"{item_id}.log"
        if backend == "vina":
            if not ligand_pdbqt or not receptor_pdbqt:
                results.append({"panel_type": panel_type, "id": item_id, "smiles": smiles, "docking_score": "", "pose_pass": "", "backend": backend, "receptor": receptor_pdbqt, "notes": "missing_ligand_or_receptor_pdbqt"})
                continue
            pose_output = run_dir / f"{item_id}_pose.pdbqt"
            command = [
                backend_exe,
                "--receptor",
                receptor_pdbqt,
                "--ligand",
                ligand_pdbqt,
                "--center_x",
                str(center[0]),
                "--center_y",
                str(center[1]),
                "--center_z",
                str(center[2]),
                "--size_x",
                str(size[0]),
                "--size_y",
                str(size[1]),
                "--size_z",
                str(size[2]),
                "--out",
                str(pose_output),
            ]
        else:
            if not receptor_pdb:
                results.append({"panel_type": panel_type, "id": item_id, "smiles": smiles, "docking_score": "", "pose_pass": "", "backend": backend, "receptor": receptor_pdb, "notes": "missing_receptor_pdb"})
                continue
            pose_output = run_dir / f"{item_id}_poses.sdf"
            command = [
                backend_exe,
                "--receptor",
                receptor_pdb,
                "--ligand",
                str(item.get("ligand_sdf", "")),
                "--center_x",
                str(center[0]),
                "--center_y",
                str(center[1]),
                "--center_z",
                str(center[2]),
                "--size_x",
                str(size[0]),
                "--size_y",
                str(size[1]),
                "--size_z",
                str(size[2]),
                "--out",
                str(pose_output),
            ]
        try:
            proc = subprocess.run(command, check=False, capture_output=True, text=True, timeout=timeout)
            log_path.write_text(
                "\n".join(
                    [
                        "$ " + " ".join(command),
                        "",
                        "[stdout]",
                        proc.stdout or "",
                        "",
                        "[stderr]",
                        proc.stderr or "",
                    ]
                ),
                encoding="utf-8",
            )
            score = parse_docking_score_from_text(proc.stdout + "\n" + proc.stderr)
            pose_quality = {}
            pose_missing_note = ""
            if score != "" and proc.returncode == 0 and not pose_output.exists():
                pose_missing_note = f"missing_pose_output={pose_output}"
            if score != "" and proc.returncode == 0 and not pose_missing_note:
                pose_quality = stage4_run_posebusters_if_possible(
                    item_id,
                    pose_output,
                    receptor_pdb,
                    run_dir,
                    executables,
                    timeout,
                )
            pose_quality_status = str(pose_quality.get("status", "")) if pose_quality else ""
            notes = f"returncode={proc.returncode};log={log_path}"
            if pose_missing_note:
                notes += f";{pose_missing_note}"
            if pose_quality_status:
                notes += f";posebusters={pose_quality_status}"
            if pose_quality.get("report"):
                notes += f";posebusters_report={pose_quality.get('report')}"
            if pose_missing_note:
                docking_plan.setdefault("rejected_docking_outputs", [])
                rejected = docking_plan["rejected_docking_outputs"]
                if isinstance(rejected, list):
                    rejected.append({"id": item_id, "reason": "missing_pose_output", "score": score, "log": str(log_path), "pose": str(pose_output)})
                continue
            results.append(
                {
                    "panel_type": panel_type,
                    "id": item_id,
                    "smiles": smiles,
                    "docking_score": score,
                    "pose_pass": "true" if score != "" and proc.returncode == 0 and pose_quality_status == "passed" else "",
                    "backend": backend,
                    "receptor": receptor_pdbqt or receptor_pdb,
                    "notes": notes,
                }
            )
        except Exception as exc:
            results.append({"panel_type": panel_type, "id": item_id, "smiles": smiles, "docking_score": "", "pose_pass": "", "backend": backend, "receptor": receptor_pdbqt or receptor_pdb, "notes": str(exc)})

    if results:
        write_csv(output_csv, results, output_fields or STAGE4_DOCKING_SCORE_FIELDS)
    docking_plan["status"] = "completed" if any(row.get("docking_score") not in {"", None} for row in results) else "attempted_no_scores"
    rejected = docking_plan.get("rejected_docking_outputs", [])
    if isinstance(rejected, list) and rejected:
        docking_plan["message"] = f"Docking command attempted but {len(rejected)} scored outputs were rejected because pose files were missing (missing_pose_output); inspect run logs."
    else:
        docking_plan["message"] = "Docking command attempted; inspect docking_scores_template/results CSV and run logs."
    return docking_plan, results


def import_stage4_external_docking_scores(
    external_scores: Optional[str],
    output_csv: Path,
    docking_plan: Dict[str, object],
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    if not external_scores:
        return docking_plan, []
    source = Path(external_scores).expanduser()
    if not source.exists():
        raise SystemExit(f"External docking scores not found: {source}")
    rows = read_csv(source)
    imported: List[Dict[str, object]] = []
    for row in rows:
        score = (
            row.get("docking_score")
            or row.get("affinity")
            or row.get("vina_score")
            or row.get("gnina_score")
            or row.get("CNNaffinity")
            or row.get("cnn_affinity")
            or ""
        )
        score = str(score).strip()
        if not score:
            continue
        imported.append(
            {
                "id": row.get("id") or row.get("candidate_id") or row.get("name") or "",
                "smiles": row.get("smiles") or row.get("canonical_smiles") or "",
                "docking_score": score,
                "pose_pass": row.get("pose_pass") or row.get("posebusters_pass") or row.get("pose_score") or "",
                "backend": row.get("backend") or row.get("tool") or row.get("source") or "external",
                "receptor": row.get("receptor") or row.get("receptor_pdb") or row.get("receptor_pdbqt") or "",
                "notes": row.get("notes") or row.get("note") or f"imported_external_scores={source}",
            }
        )
    write_csv(output_csv, imported, STAGE4_DOCKING_SCORE_FIELDS)
    if imported:
        docking_plan["status"] = "imported_external_scores"
        docking_plan["message"] = f"Imported {len(imported)} external docking score rows from {source}."
        docking_plan["external_scores_csv"] = str(source)
    else:
        docking_plan["status"] = "external_scores_empty"
        docking_plan["message"] = f"External docking score CSV had no usable score rows: {source}."
        docking_plan["external_scores_csv"] = str(source)
    return docking_plan, imported


def parse_docking_score_from_text(text: str) -> str:
    vina_match = re.search(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+", text, flags=re.MULTILINE)
    if vina_match:
        return vina_match.group(1)
    affinity_match = re.search(r"affinity[:\s]+(-?\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)
    if affinity_match:
        return affinity_match.group(1)
    cnn_match = re.search(r"CNNaffinity[:\s]+(-?\d+(?:\.\d+)?)", text, flags=re.IGNORECASE)
    if cnn_match:
        return cnn_match.group(1)
    return ""


def write_stage4_validation_metrics(
    path: Path,
    benchmark_rows: Sequence[Dict[str, object]],
    similarity_rows: Sequence[Dict[str, object]],
    docking_plan: Dict[str, object],
) -> Dict[str, object]:
    counts = Counter(str(row.get("panel_type", "")) for row in benchmark_rows)
    candidate_sims = [float(row.get("tanimoto", 0.0)) for row in similarity_rows if row.get("candidate_id")]
    metrics = {
        "schema_version": "0.1",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "panel_counts": dict(counts),
        "control_similarity": {
            "candidate_control_pairs": len(candidate_sims),
            "max_candidate_to_control_tanimoto": round(max(candidate_sims), 4) if candidate_sims else "",
            "mean_candidate_to_control_tanimoto": round(sum(candidate_sims) / len(candidate_sims), 4) if candidate_sims else "",
            "interpretation": "Fingerprint similarity is a benchmark for chemical relatedness to controls, not evidence of antiviral activity.",
        },
        "docking_status": docking_plan.get("status", ""),
        "docking_backend": docking_plan.get("selected_backend", ""),
        "readiness": {
            "rdkit_descriptors": "ready" if counts.get("candidate", 0) else "missing_candidates",
            "control_panel": "ready" if counts.get("positive_control", 0) or counts.get("reference_control", 0) else "missing_controls",
            "decoy_panel": "ready" if counts.get("decoy", 0) else "missing_decoys",
            "docking": docking_plan.get("status", ""),
        },
    }
    write_json(path, metrics)
    return metrics


def render_stage4_report(
    round_no: int,
    target_id: str,
    descriptor_rows: Sequence[Dict[str, object]],
    controls: Sequence[Dict[str, object]],
    similarity_rows: Sequence[Dict[str, object]],
    diversity_rows: Sequence[Dict[str, object]],
    capabilities: Dict[str, object],
    files: Dict[str, str],
    sdf_written: int,
    sdf_warnings: Sequence[Dict[str, str]],
    receptor_package: Optional[Dict[str, object]] = None,
    docking_plan: Optional[Dict[str, object]] = None,
    validation_metrics: Optional[Dict[str, object]] = None,
    rendered_2d_count: int = 0,
) -> str:
    valid_count = sum(1 for row in descriptor_rows if str(row.get("valid", "")) == "1")
    invalid_count = len(descriptor_rows) - valid_count
    rdkit_status = capabilities.get("rdkit", {})
    if not isinstance(rdkit_status, dict):
        rdkit_status = {}
    docking = capabilities.get("docking_backend", {})
    if not isinstance(docking, dict):
        docking = {}
    receptor_package = receptor_package or {}
    docking_plan = docking_plan or {}
    validation_metrics = validation_metrics or {}
    panel_counts = validation_metrics.get("panel_counts", {})
    if not isinstance(panel_counts, dict):
        panel_counts = {}
    best_similarity = sorted(similarity_rows, key=lambda row: float(row.get("tanimoto", 0.0)), reverse=True)[:10]
    lines = [
        f"# Stage 4 Real Library Validation - Round {round_no}",
        "",
        "## Summary",
        "",
        f"- Target: `{target_id or 'not_specified'}`",
        f"- Candidate rows: {len(descriptor_rows)}",
        f"- RDKit-valid molecules: {valid_count}",
        f"- Invalid molecules: {invalid_count}",
        f"- Controls with structures: {len(controls)}",
        f"- Benchmark panel: candidates `{panel_counts.get('candidate', 0)}`, controls `{int(panel_counts.get('positive_control', 0) or 0) + int(panel_counts.get('reference_control', 0) or 0)}`, decoys `{panel_counts.get('decoy', 0)}`",
        f"- SDF ligands written: {sdf_written}",
        f"- 2D images rendered: {rendered_2d_count}",
        f"- Receptor preparation status: `{receptor_package.get('preparation_status', '')}`",
        f"- RDKit status: `{rdkit_status.get('status', '')}` version `{rdkit_status.get('version', '')}`",
        f"- Docking backend status: `{docking_plan.get('status', docking.get('status', ''))}`",
        f"- Docking backend message: {docking_plan.get('message', docking.get('message', ''))}",
        "",
        "## Output Files",
        "",
    ]
    for label, path in files.items():
        lines.append(f"- {label}: `{path}`")
    lines.extend(
        [
            "",
            "## Top Candidate-Control Similarities",
            "",
            "| candidate | control | Tanimoto | role |",
            "|---|---|---:|---|",
        ]
    )
    for row in best_similarity:
        lines.append(
            f"| {row.get('candidate_id', '')} | {row.get('control_drug', '')} | "
            f"{row.get('tanimoto', '')} | {row.get('control_role', '')} |"
        )
    lines.extend(
        [
            "",
            "## Diverse RDKit-Ready Selection",
            "",
            "| selection | id | QED | max control sim | nearest control | scaffold |",
            "|---:|---|---:|---:|---|---|",
        ]
    )
    for row in diversity_rows:
        lines.append(
            f"| {row.get('selection_rank', '')} | {row.get('id', '')} | {row.get('qed', '')} | "
            f"{row.get('max_similarity_to_controls', '')} | {row.get('nearest_control', '')} | `{row.get('murcko_scaffold', '')}` |"
        )
    if sdf_warnings:
        lines.extend(["", "## SDF Warnings", ""])
        for warning in sdf_warnings[:20]:
            lines.append(f"- {warning.get('id', '')}: {warning.get('warning', '')}")
    lines.extend(
        [
            "",
            "## Receptor And Docking Readiness",
            "",
            f"- PDB ID: `{receptor_package.get('pdb_id', '')}`",
            f"- Local receptor PDB: `{receptor_package.get('local_receptor_pdb', '')}`",
            f"- Receptor PDBQT: `{receptor_package.get('local_receptor_pdbqt', '')}`",
            f"- Docking plan status: `{docking_plan.get('status', '')}`",
            f"- Selected backend: `{docking_plan.get('selected_backend', '')}`",
            f"- Expected docking score CSV: `{docking_plan.get('expected_scores_csv', '')}`",
        ]
    )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- Stage 4 uses real RDKit chemistry objects, descriptors, fingerprints, scaffolds, and SDF export.",
            "- Similarity to a known drug is a computational benchmark, not proof of activity or safety.",
            "- SDF export prepares ligands for docking, but no docking score is claimed unless Vina/GNINA is actually installed and run.",
            "- Next product step: install/configure Meeko or OpenBabel plus Vina/GNINA, then feed real docking CSV back through `score --external-scores`.",
        ]
    )
    return "\n".join(lines) + "\n"


def stage4_docking_score_band(score: object) -> str:
    try:
        value = float(score)
    except Exception:
        return "missing"
    if value <= -7.5:
        return "strong"
    if value <= -6.0:
        return "moderate"
    return "weak"


def stage45_score_value(row: Dict[str, object]) -> Optional[float]:
    try:
        raw = row.get("docking_score", "")
        if raw in {"", None}:
            return None
        return float(raw)
    except Exception:
        return None


def stage45_best_row(rows: Sequence[Dict[str, object]], panel_types: Sequence[str]) -> Dict[str, object]:
    allowed = set(panel_types)
    scored = [row for row in rows if str(row.get("panel_type", "")) in allowed and stage45_score_value(row) is not None]
    if not scored:
        return {}
    scored.sort(key=lambda row: (float(stage45_score_value(row) or 0.0), str(row.get("id", ""))))
    best = dict(scored[0])
    value = stage45_score_value(best)
    if value is not None:
        best["docking_score"] = round(value, 4)
    return best


def stage45_panel_counts(rows: Sequence[Dict[str, object]]) -> Dict[str, int]:
    raw = Counter(str(row.get("panel_type", "")) for row in rows)
    known = sum(count for panel, count in raw.items() if panel not in {"candidate", "decoy", ""})
    return {
        "candidate": int(raw.get("candidate", 0)),
        "known_control": int(known),
        "decoy": int(raw.get("decoy", 0)),
        "by_panel_type": dict(raw),
    }


def stage45_score_summary(rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    by_panel: Dict[str, List[float]] = defaultdict(list)
    for row in rows:
        value = stage45_score_value(row)
        if value is None:
            continue
        by_panel[str(row.get("panel_type", ""))].append(value)
    summary: Dict[str, object] = {}
    for panel, values in by_panel.items():
        summary[panel] = {
            "count": len(values),
            "best": round(min(values), 4),
            "worst": round(max(values), 4),
            "mean": round(sum(values) / len(values), 4),
        }
    return summary


def stage45_load_candidate_panel(project: Path, round_no: int, top: int) -> List[Dict[str, object]]:
    benchmark = read_csv(project / "stage4" / f"round_{round_no}_benchmark_panel.csv")
    candidates_by_id = {row.get("id", ""): row for row in benchmark if row.get("panel_type") == "candidate"}
    ranked = read_csv(project / "ranked" / f"round_{round_no}_ranked.csv")
    selected: List[Dict[str, object]] = []
    if ranked:
        for row in ranked:
            item_id = row.get("id", "")
            source = candidates_by_id.get(item_id)
            if source:
                selected.append(dict(source))
            elif row.get("smiles"):
                selected.append(
                    {
                        "panel_type": "candidate",
                        "id": item_id,
                        "smiles": row.get("smiles", ""),
                        "canonical_smiles": row.get("smiles", ""),
                        "role": "candidate",
                        "target_id": "",
                        "source": "ranked",
                    }
                )
            if len(selected) >= top:
                break
    if not selected:
        selected = [dict(row) for row in benchmark if row.get("panel_type") == "candidate"][:top]
    if not selected:
        descriptors = read_csv(project / "stage4" / f"round_{round_no}_real_descriptors.csv")
        for row in descriptors:
            if str(row.get("valid", "")) != "1":
                continue
            selected.append(
                {
                    "panel_type": "candidate",
                    "id": row.get("id", ""),
                    "smiles": row.get("smiles", ""),
                    "canonical_smiles": row.get("canonical_smiles", ""),
                    "role": "candidate",
                    "target_id": "",
                    "qed": row.get("qed", ""),
                    "mw": row.get("mw", ""),
                    "logp": row.get("logp", ""),
                    "tpsa": row.get("tpsa", ""),
                    "hbd": row.get("hbd", ""),
                    "hba": row.get("hba", ""),
                    "heavy_atoms": row.get("heavy_atoms", ""),
                    "murcko_scaffold": row.get("murcko_scaffold", ""),
                    "source": row.get("source", "stage4_descriptors"),
                }
            )
            if len(selected) >= top:
                break
    return selected[:top]


def stage45_control_panel_rows(
    controls: Sequence[Dict[str, object]],
    target_id: str,
    libs: Dict[str, object],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for control in controls:
        desc = stage4_descriptor_for_smiles(
            str(control.get("role", "control")) or "control",
            str(control.get("drug", "")),
            str(control.get("smiles", "")),
            str(control.get("role", "")),
            target_id,
            str(control.get("source", "")),
            libs,
        )
        if str(desc.get("valid", "")) != "1":
            continue
        rows.append(
            {
                "panel_type": str(control.get("role", "control")) or "control",
                "id": control.get("drug", ""),
                "smiles": control.get("smiles", ""),
                "canonical_smiles": desc.get("canonical_smiles", ""),
                "role": control.get("role", ""),
                "target_id": target_id,
                "qed": desc.get("qed", ""),
                "mw": desc.get("mw", ""),
                "logp": desc.get("logp", ""),
                "tpsa": desc.get("tpsa", ""),
                "hbd": desc.get("hbd", ""),
                "hba": desc.get("hba", ""),
                "heavy_atoms": desc.get("heavy_atoms", ""),
                "murcko_scaffold": desc.get("murcko_scaffold", ""),
                "max_similarity_to_controls": 1.0,
                "nearest_control": control.get("drug", ""),
                "source": control.get("source", ""),
            }
        )
    return rows


def stage45_decoy_panel_rows(
    decoy_rows: Sequence[Dict[str, object]],
    controls: Sequence[Dict[str, object]],
    target_id: str,
    libs: Dict[str, object],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for decoy in decoy_rows:
        if str(decoy.get("valid", "")) != "1":
            continue
        max_sim, nearest = nearest_control_summary(str(decoy.get("canonical_smiles") or decoy.get("smiles", "")), controls, libs)
        rows.append(
            {
                "panel_type": "decoy",
                "id": decoy.get("id", ""),
                "smiles": decoy.get("smiles", ""),
                "canonical_smiles": decoy.get("canonical_smiles", ""),
                "role": "decoy",
                "target_id": target_id,
                "qed": decoy.get("qed", ""),
                "mw": decoy.get("mw", ""),
                "logp": decoy.get("logp", ""),
                "tpsa": decoy.get("tpsa", ""),
                "hbd": decoy.get("hbd", ""),
                "hba": decoy.get("hba", ""),
                "heavy_atoms": decoy.get("heavy_atoms", ""),
                "murcko_scaffold": decoy.get("murcko_scaffold", ""),
                "max_similarity_to_controls": max_sim,
                "nearest_control": nearest,
                "source": decoy.get("source", ""),
            }
        )
    return rows


def stage45_export_reference_ligand_pdb(receptor_package: Dict[str, object], output_dir: Path) -> Dict[str, object]:
    receptor_pdb = Path(str(receptor_package.get("local_receptor_pdb", "")))
    site = receptor_package.get("binding_site", {})
    if not receptor_pdb.exists() or not isinstance(site, dict):
        return {"status": "unavailable", "path": "", "note": "No local receptor PDB or binding-site metadata was available."}
    ligand = site.get("detected_ligand", {})
    if not isinstance(ligand, dict):
        return {"status": "unavailable", "path": "", "note": "No co-crystal ligand was recorded in the receptor package."}
    resname = str(ligand.get("resname", "")).upper()
    chain = str(ligand.get("chain", ""))
    resseq = str(ligand.get("resseq", ""))
    icode = str(ligand.get("icode", ""))
    if not resname or not resseq:
        return {"status": "unavailable", "path": "", "note": "Co-crystal ligand metadata lacks residue identity."}
    kept: List[str] = []
    for line in receptor_pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("HETATM"):
            continue
        line_resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        line_chain = line[21:22].strip() if len(line) >= 22 else ""
        line_resseq = line[22:26].strip() if len(line) >= 26 else ""
        line_icode = line[26:27].strip() if len(line) >= 27 else ""
        alt_loc = line[16:17].strip() if len(line) > 16 else ""
        if alt_loc and alt_loc not in {"A", "1"}:
            continue
        if line_resname == resname and line_chain == chain and line_resseq == resseq and line_icode == icode:
            kept.append(line)
    if not kept:
        return {"status": "unavailable", "path": "", "note": "No matching HETATM records were found for the detected ligand."}
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{resname}_{chain}_{resseq}{icode}").strip("_")
    path = output_dir / f"reference_ligand_{safe}.pdb"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(kept + ["END"]) + "\n", encoding="utf-8")
    return {
        "status": "reference_pose_exported",
        "path": str(path),
        "resname": resname,
        "chain": chain,
        "resseq": resseq,
        "icode": icode,
        "atom_count": len(kept),
        "rmsd_status": "unavailable",
        "rmsd_note": "RMSD is not reported unless atom mapping between the crystallographic ligand and the docked pose is robustly established.",
    }


def stage45_redock_reference_ligand(
    reference_export: Dict[str, object],
    receptor_package: Dict[str, object],
    docking_plan: Dict[str, object],
    capabilities: Dict[str, object],
    output_dir: Path,
    timeout: int,
) -> Dict[str, object]:
    result = dict(reference_export)
    if reference_export.get("status") != "reference_pose_exported":
        return result
    if docking_plan.get("status") not in {"planned", "completed"}:
        result.update({"redocking_status": "not_run", "redocking_note": "Docking was not planned for Stage 4.5."})
        return result
    executables = capabilities.get("executables", {})
    if not isinstance(executables, dict):
        executables = {}
    obabel = executables.get("obabel", {}) if isinstance(executables.get("obabel", {}), dict) else {}
    obabel_path = str(obabel.get("path", ""))
    backend = str(docking_plan.get("selected_backend", ""))
    backend_status = executables.get(backend, {}) if isinstance(executables.get(backend, {}), dict) else {}
    backend_path = str(backend_status.get("path", "")) or backend
    receptor_pdbqt = str(receptor_package.get("local_receptor_pdbqt", ""))
    receptor_pdb = str(receptor_package.get("local_receptor_pdb", ""))
    ligand_pdb = Path(str(reference_export.get("path", "")))
    if backend != "vina":
        result.update({"redocking_status": "not_run", "redocking_note": "Reference redocking is currently implemented for Vina PDBQT inputs."})
        return result
    if not obabel_path or not receptor_pdbqt or not Path(receptor_pdbqt).exists():
        result.update({"redocking_status": "not_run", "redocking_note": "Missing OpenBabel or receptor PDBQT for reference ligand redocking."})
        return result
    ligand_pdbqt = output_dir / f"{ligand_pdb.stem}.pdbqt"
    ligand_sdf = output_dir / f"{ligand_pdb.stem}.sdf"
    try:
        to_pdbqt = subprocess.run([obabel_path, str(ligand_pdb), "-O", str(ligand_pdbqt), "-h"], check=False, capture_output=True, text=True, timeout=max(1, timeout))
        to_sdf = subprocess.run([obabel_path, str(ligand_pdb), "-O", str(ligand_sdf), "-h"], check=False, capture_output=True, text=True, timeout=max(1, timeout))
    except Exception as exc:
        result.update({"redocking_status": "conversion_failed", "redocking_note": str(exc)})
        return result
    if to_pdbqt.returncode != 0 or not ligand_pdbqt.exists():
        result.update({"redocking_status": "conversion_failed", "redocking_note": (to_pdbqt.stderr or to_pdbqt.stdout or "OpenBabel PDBQT conversion failed")[-500:]})
        return result
    docking_box = docking_plan.get("docking_box", {})
    if not isinstance(docking_box, dict):
        docking_box = {}
    center = docking_box.get("center", [])
    size = docking_box.get("size", [])
    if not isinstance(center, list) or not isinstance(size, list) or len(center) != 3 or len(size) != 3:
        result.update({"redocking_status": "not_run", "redocking_note": "Docking box was missing."})
        return result
    run_dir = output_dir / "reference_redocking"
    run_dir.mkdir(parents=True, exist_ok=True)
    pose = run_dir / f"{ligand_pdb.stem}_pose.pdbqt"
    log = run_dir / f"{ligand_pdb.stem}.log"
    command = [
        backend_path,
        "--receptor",
        receptor_pdbqt,
        "--ligand",
        str(ligand_pdbqt),
        "--center_x",
        str(center[0]),
        "--center_y",
        str(center[1]),
        "--center_z",
        str(center[2]),
        "--size_x",
        str(size[0]),
        "--size_y",
        str(size[1]),
        "--size_z",
        str(size[2]),
        "--out",
        str(pose),
    ]
    try:
        proc = subprocess.run(command, check=False, capture_output=True, text=True, timeout=max(1, timeout))
        log.write_text("$ " + " ".join(command) + "\n\n[stdout]\n" + (proc.stdout or "") + "\n\n[stderr]\n" + (proc.stderr or ""), encoding="utf-8")
    except Exception as exc:
        result.update({"redocking_status": "failed", "redocking_note": str(exc), "command": command})
        return result
    score = parse_docking_score_from_text((proc.stdout or "") + "\n" + (proc.stderr or ""))
    posebusters = {}
    if score and proc.returncode == 0:
        posebusters = stage4_run_posebusters_if_possible(ligand_pdb.stem, pose, receptor_pdb, run_dir, executables, timeout)
    result.update(
        {
            "redocking_status": "completed" if score else "attempted_no_score",
            "docking_score": score,
            "score_band": stage4_docking_score_band(score),
            "pose_pass": "true" if score and posebusters.get("status") == "passed" else "",
            "ligand_pdbqt": str(ligand_pdbqt),
            "ligand_sdf": str(ligand_sdf) if ligand_sdf.exists() or to_sdf.returncode == 0 else "",
            "pose": str(pose),
            "log": str(log),
            "command": command,
            "posebusters": posebusters,
        }
    )
    return result


def stage45_enrich_docking_results(rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    known_types = {str(row.get("panel_type", "")) for row in rows if str(row.get("panel_type", "")) not in {"candidate", "decoy", ""}}
    best_control = stage45_best_row(rows, sorted(known_types))
    best_control_score = stage45_score_value(best_control) if best_control else None
    enriched: List[Dict[str, object]] = []
    for row in rows:
        out = dict(row)
        value = stage45_score_value(out)
        out["score_band"] = stage4_docking_score_band(out.get("docking_score", ""))
        if value is not None and best_control_score is not None:
            out["relative_to_best_control"] = round(value - best_control_score, 4)
        else:
            out["relative_to_best_control"] = ""
        enriched.append(out)
    return enriched


def render_stage45_report(validation: Dict[str, object], rows: Sequence[Dict[str, object]]) -> str:
    counts = validation.get("counts", {})
    if not isinstance(counts, dict):
        counts = {}
    best_control = validation.get("best_known_control", {})
    best_candidate = validation.get("best_candidate", {})
    best_decoy = validation.get("best_decoy", {})
    candidate_vs = validation.get("candidate_vs_controls", {})
    redocking = validation.get("redocking", {})
    if not isinstance(best_control, dict):
        best_control = {}
    if not isinstance(best_candidate, dict):
        best_candidate = {}
    if not isinstance(best_decoy, dict):
        best_decoy = {}
    if not isinstance(candidate_vs, dict):
        candidate_vs = {}
    if not isinstance(redocking, dict):
        redocking = {}
    lines = [
        f"# Stage 4.5 Control Calibration - Round {validation.get('round', '')}",
        "",
        "## Summary",
        "",
        f"- Target: `{validation.get('target_id', '')}`",
        f"- Candidates docked: {counts.get('candidate', 0)}",
        f"- Known controls docked: {counts.get('known_control', 0)}",
        f"- Decoys docked: {counts.get('decoy', 0)}",
        f"- Docking status: `{nested_get(validation, ['docking', 'status']) or ''}`",
        f"- Best known control: `{best_control.get('id', '')}` score `{best_control.get('docking_score', '')}`",
        f"- Best candidate: `{best_candidate.get('id', '')}` score `{best_candidate.get('docking_score', '')}`",
        f"- Best decoy: `{best_decoy.get('id', '')}` score `{best_decoy.get('docking_score', '')}`",
        f"- Best candidate delta to best control: `{candidate_vs.get('best_candidate_delta_kcal_mol', '')}` kcal/mol",
        "",
        "known controls calibrate the docking workflow; they do not make generated candidates active by analogy.",
        "",
        "## Reference Ligand",
        "",
        f"- Export status: `{redocking.get('status', '')}`",
        f"- Redocking status: `{redocking.get('redocking_status', '')}`",
        f"- Reference ligand file: `{redocking.get('path', '')}`",
        f"- RMSD status: `{redocking.get('rmsd_status', 'unavailable')}`",
        "",
        "## Docking Pool",
        "",
        "| panel | id | score | band | pose pass | delta to best control |",
        "|---|---|---:|---|---|---:|",
    ]
    for row in rows[:80]:
        lines.append(
            f"| {row.get('panel_type', '')} | {row.get('id', '')} | {row.get('docking_score', '')} | "
            f"{row.get('score_band', '')} | {row.get('pose_pass', '')} | {row.get('relative_to_best_control', '')} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- This control layer checks whether candidates, known controls, and decoys behave plausibly under the same computational setup.",
            "- A candidate near or better than a known control is a prioritization signal only; it does not prove antiviral efficacy.",
            "- Docking and PoseBusters are computational filters and must be followed by assay design, orthogonal modeling, and experimental validation.",
        ]
    )
    return "\n".join(lines) + "\n"


def stage45_validate_controls(args: argparse.Namespace) -> None:
    project = project_path(args.project)
    ensure_project_dirs(project)
    config = load_config(project)
    target = config.get("target", {})
    if not isinstance(target, dict):
        target = {}
    target_id = args.target or str(target.get("target_catalog_id", "") or "")
    if not target_id:
        stage4_assets = read_json(project / "stage4" / f"round_{args.round}_stage4_assets.json") if (project / "stage4" / f"round_{args.round}_stage4_assets.json").exists() else {}
        target_id = str(stage4_assets.get("target_id", ""))
    if not target_id:
        raise SystemExit("Stage 4.5 requires --target or a Stage 4 assets target_id.")

    stage4_dir = project / "stage4"
    stage45_dir = project / "stage4_5"
    stage45_dir.mkdir(parents=True, exist_ok=True)
    receptor_package_path = stage4_dir / f"round_{args.round}_receptor_package.json"
    if not receptor_package_path.exists():
        raise SystemExit(f"Stage 4 receptor package is required before Stage 4.5: {receptor_package_path}")
    receptor_package = read_json(receptor_package_path)

    libs = require_rdkit()
    capabilities = stage4_capabilities()
    rdkit_status = capabilities.get("rdkit", {})
    if isinstance(rdkit_status, dict):
        rdkit_status["status"] = "available"
        rdkit_status["version"] = str(getattr(libs["rdBase"], "rdkitVersion", rdkit_status.get("version", "")))

    candidates = stage45_load_candidate_panel(project, args.round, max(1, int(getattr(args, "top_candidates", 12) or 12)))
    controls = load_stage4_controls(target_id, getattr(args, "controls_csv", None))
    control_rows = stage45_control_panel_rows(controls, target_id, libs)
    decoy_rows = build_stage4_decoys(int(getattr(args, "decoys", 8) or 0), libs)
    decoy_panel = stage45_decoy_panel_rows(decoy_rows, controls, target_id, libs)
    panel_rows = candidates + control_rows + decoy_panel
    if not panel_rows:
        raise SystemExit("Stage 4.5 found no candidate/control/decoy rows to validate.")
    for row in panel_rows:
        row["target_id"] = row.get("target_id") or target_id

    panel_path = stage45_dir / f"round_{args.round}_control_panel.csv"
    docking_inputs_path = stage45_dir / f"round_{args.round}_control_docking_inputs.csv"
    docking_plan_path = stage45_dir / f"round_{args.round}_control_docking_plan.json"
    scores_path = stage45_dir / f"round_{args.round}_control_docking_scores.csv"
    validation_path = stage45_dir / f"round_{args.round}_control_validation.json"
    report_path = project / "reports" / f"stage4_5_round_{args.round}_control_validation.md"
    reference_dir = stage45_dir

    write_csv(panel_path, panel_rows, STAGE4_BENCHMARK_FIELDS)
    docking_timeout = int(getattr(args, "docking_timeout", 600) or 600)
    docking_box = stage4_box_from_receptor_package(receptor_package)
    run_docking = not bool(getattr(args, "no_docking", False))
    receptor_preparation = stage4_prepare_receptor_pdbqt(receptor_package, capabilities, docking_box, docking_timeout) if run_docking else {}
    if receptor_preparation.get("status") == "ready" and receptor_preparation.get("receptor_pdbqt"):
        receptor_package["local_receptor_pdbqt"] = receptor_preparation["receptor_pdbqt"]
        write_json(receptor_package_path, receptor_package)

    dockable_types = ["candidate", "positive_control", "reference_control", "historical_control_not_recommended", "control", "decoy"]
    ligand_records, ligand_attempts = stage4_prepare_ligand_records(
        panel_rows,
        libs,
        capabilities,
        stage45_dir / f"round_{args.round}_prepared_ligands",
        int(getattr(args, "seed", 61453) or 61453),
        docking_timeout,
        run_docking,
        dockable_panel_types=dockable_types,
    )
    docking_input_rows = stage4_docking_inputs(
        panel_rows,
        receptor_package,
        stage4_dir / f"round_{args.round}_ligands.sdf",
        ligand_records,
        dockable_panel_types=dockable_types,
    )
    write_csv(docking_inputs_path, docking_input_rows, STAGE4_DOCKING_INPUT_FIELDS)
    write_csv(scores_path, [], STAGE45_DOCKING_SCORE_FIELDS)
    docking_plan = write_stage4_docking_plan(
        docking_plan_path,
        str(getattr(args, "docking_backend", "auto") or "auto"),
        run_docking,
        capabilities,
        receptor_package,
        docking_inputs_path,
        scores_path,
        docking_box,
        receptor_preparation,
        sum(1 for item in ligand_records.values() if item.get("ligand_pdbqt")),
    )
    docking_plan, docking_results = run_stage4_docking_if_possible(
        docking_plan,
        docking_input_rows,
        panel_rows,
        scores_path,
        docking_timeout,
        dockable_panel_types=dockable_types,
        output_fields=STAGE45_DOCKING_SCORE_FIELDS,
        run_dir_name="control_docking_runs",
    )
    docking_plan["ligand_preparation_attempts"] = ligand_attempts[:80]
    write_json(docking_plan_path, docking_plan)
    enriched_results = stage45_enrich_docking_results(docking_results)
    write_csv(scores_path, enriched_results, STAGE45_DOCKING_SCORE_FIELDS)

    reference_export = stage45_export_reference_ligand_pdb(receptor_package, reference_dir)
    redocking = stage45_redock_reference_ligand(reference_export, receptor_package, docking_plan, capabilities, reference_dir, docking_timeout) if run_docking else reference_export
    counts = stage45_panel_counts(panel_rows)
    scored_count = sum(1 for row in enriched_results if stage45_score_value(row) is not None)
    pose_pass_count = sum(1 for row in enriched_results if str(row.get("pose_pass", "")).lower() == "true")
    known_types = sorted({str(row.get("panel_type", "")) for row in panel_rows if str(row.get("panel_type", "")) not in {"candidate", "decoy", ""}})
    best_control = stage45_best_row(enriched_results, known_types)
    best_candidate = stage45_best_row(enriched_results, ["candidate"])
    best_decoy = stage45_best_row(enriched_results, ["decoy"])
    control_score = stage45_score_value(best_control) if best_control else None
    candidate_score = stage45_score_value(best_candidate) if best_candidate else None
    delta = round(candidate_score - control_score, 4) if candidate_score is not None and control_score is not None else ""
    validation = {
        "schema_version": "0.1",
        "stage": 4.5,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "round": args.round,
        "target_id": target_id,
        "counts": counts,
        "docking": {
            "status": docking_plan.get("status", ""),
            "backend": docking_plan.get("selected_backend", ""),
            "scored_count": scored_count,
            "pose_pass_count": pose_pass_count,
            "pose_pass_rate": round(pose_pass_count / scored_count, 4) if scored_count else "",
            "docking_box": docking_box,
        },
        "score_summary": stage45_score_summary(enriched_results),
        "best_known_control": best_control,
        "best_candidate": best_candidate,
        "best_decoy": best_decoy,
        "candidate_vs_controls": {
            "best_candidate_delta_kcal_mol": delta,
            "candidate_beats_or_matches_best_control": bool(candidate_score is not None and control_score is not None and candidate_score <= control_score),
            "interpretation": "Negative or near-zero delta is a computational prioritization signal only; it is not a measured activity claim.",
        },
        "redocking": redocking,
        "files": {
            "control_panel": str(panel_path),
            "docking_inputs": str(docking_inputs_path),
            "docking_plan": str(docking_plan_path),
            "docking_scores": str(scores_path),
            "validation": str(validation_path),
            "report": str(report_path),
        },
        "boundary": [
            "Stage 4.5 is a computational control-calibration layer, not a wet-lab validation result.",
            "Known controls calibrate the scoring workflow and reveal whether decoys/candidates behave plausibly under the same setup.",
            "Docking scores, pose checks, and reference redocking do not prove antiviral efficacy, selectivity, toxicity, or synthesizability.",
        ],
    }
    write_json(validation_path, validation)
    write_text(report_path, render_stage45_report(validation, enriched_results))
    print(f"Wrote Stage 4.5 control panel: {panel_path}")
    print(f"Wrote Stage 4.5 docking scores: {scores_path}")
    print(f"Wrote Stage 4.5 validation summary: {validation_path}")
    print(f"Wrote Stage 4.5 report: {report_path}")


def stage46_parse_types(value: object, default: Sequence[str]) -> List[str]:
    raw = str(value or "").strip()
    if not raw:
        return list(default)
    items = [item.strip() for item in re.split(r"[,;\s]+", raw) if item.strip()]
    seen = set()
    parsed: List[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            parsed.append(item)
    return parsed or list(default)


def stage46_parse_top_k(value: object) -> List[int]:
    raw = str(value or "").strip()
    if not raw:
        return list(STAGE46_DEFAULT_TOP_K)
    parsed: List[int] = []
    for item in re.split(r"[,;\s]+", raw):
        if not item.strip():
            continue
        try:
            k = int(item)
        except ValueError:
            continue
        if k > 0 and k not in parsed:
            parsed.append(k)
    return parsed or list(STAGE46_DEFAULT_TOP_K)


def stage46_label_for_panel(panel_type: str, positive_types: Sequence[str], negative_types: Sequence[str]) -> str:
    if panel_type in set(positive_types):
        return "positive"
    if panel_type in set(negative_types):
        return "negative"
    if panel_type == "candidate":
        return "candidate"
    return "ignored"


def stage46_score_summary(values: Sequence[float]) -> Dict[str, object]:
    if not values:
        return {"count": 0, "best": "", "worst": "", "mean": ""}
    return {
        "count": len(values),
        "best": round(min(values), 4),
        "worst": round(max(values), 4),
        "mean": round(sum(values) / len(values), 4),
    }


def stage46_auc(positive_scores: Sequence[float], negative_scores: Sequence[float]) -> object:
    if not positive_scores or not negative_scores:
        return ""
    wins = 0.0
    pairs = 0
    for positive in positive_scores:
        for negative in negative_scores:
            pairs += 1
            if positive < negative:
                wins += 1.0
            elif positive == negative:
                wins += 0.5
    return round(wins / pairs, 4) if pairs else ""


def stage46_rank_rows(
    source_rows: Sequence[Dict[str, object]],
    positive_types: Sequence[str],
    negative_types: Sequence[str],
) -> List[Dict[str, object]]:
    scored: List[Tuple[float, Dict[str, object]]] = []
    for row in source_rows:
        score = stage45_score_value(row)
        if score is None:
            continue
        scored.append((score, dict(row)))
    scored.sort(key=lambda item: (item[0], str(item[1].get("id", ""))))
    best_score = scored[0][0] if scored else None
    positive_scores = [score for score, row in scored if stage46_label_for_panel(str(row.get("panel_type", "")), positive_types, negative_types) == "positive"]
    best_control_score = min(positive_scores) if positive_scores else None
    ranking: List[Dict[str, object]] = []
    for idx, (score, row) in enumerate(scored, start=1):
        panel_type = str(row.get("panel_type", ""))
        relative_to_best_control = row.get("relative_to_best_control", "")
        if (relative_to_best_control in {"", None}) and best_control_score is not None:
            relative_to_best_control = round(score - best_control_score, 4)
        ranking.append(
            {
                "rank": idx,
                "id": row.get("id", ""),
                "panel_type": panel_type,
                "label": stage46_label_for_panel(panel_type, positive_types, negative_types),
                "docking_score": round(score, 4),
                "score_band": row.get("score_band", "") or stage4_docking_score_band(score),
                "pose_pass": row.get("pose_pass", ""),
                "backend": row.get("backend", ""),
                "receptor": row.get("receptor", ""),
                "relative_to_best_control": relative_to_best_control,
                "relative_to_best_row": round(score - best_score, 4) if best_score is not None else "",
                "smiles": row.get("smiles", ""),
                "notes": row.get("notes", ""),
            }
        )
    return ranking


def stage46_top_k_metrics(ranking: Sequence[Dict[str, object]], top_k: Sequence[int], positives: int) -> Dict[str, object]:
    metrics: Dict[str, object] = {}
    total = len(ranking)
    background_positive_rate = (positives / total) if total else 0.0
    for k in top_k:
        top_rows = list(ranking[: min(k, total)])
        denominator = len(top_rows)
        control_hits = sum(1 for row in top_rows if row.get("label") == "positive")
        negative_hits = sum(1 for row in top_rows if row.get("label") == "negative")
        candidate_hits = sum(1 for row in top_rows if row.get("label") == "candidate")
        hit_rate = round(control_hits / positives, 4) if positives else ""
        top_fraction = (control_hits / denominator) if denominator else 0.0
        enrichment = round(top_fraction / background_positive_rate, 4) if background_positive_rate else ""
        metrics[f"top_{k}"] = {
            "k": k,
            "evaluated": denominator,
            "control_hit_count": control_hits,
            "negative_count": negative_hits,
            "candidate_count": candidate_hits,
            "control_hit_rate": hit_rate,
            "enrichment_factor": enrichment,
        }
    return metrics


def stage46_best_candidate(ranking: Sequence[Dict[str, object]]) -> Dict[str, object]:
    for row in ranking:
        if row.get("label") == "candidate":
            return {
                "rank": row.get("rank", ""),
                "id": row.get("id", ""),
                "docking_score": row.get("docking_score", ""),
                "score_band": row.get("score_band", ""),
                "relative_to_best_control": row.get("relative_to_best_control", ""),
            }
    return {}


def render_stage46_report(benchmark: Dict[str, object], ranking: Sequence[Dict[str, object]]) -> str:
    counts = benchmark.get("counts", {})
    metrics = benchmark.get("metrics", {})
    top_k = metrics.get("top_k", {}) if isinstance(metrics, dict) else {}
    best_candidate = metrics.get("best_candidate", {}) if isinstance(metrics, dict) else {}
    if not isinstance(counts, dict):
        counts = {}
    if not isinstance(top_k, dict):
        top_k = {}
    if not isinstance(best_candidate, dict):
        best_candidate = {}
    lines = [
        f"# Stage 4.6 Retrospective Benchmark - Round {benchmark.get('round', '')}",
        "",
        "## Summary",
        "",
        f"- Source scores: `{benchmark.get('source_scores', '')}`",
        f"- Positives: {counts.get('positives', 0)}",
        f"- Negatives: {counts.get('negatives', 0)}",
        f"- Candidates carried as context: {counts.get('candidates', 0)}",
        f"- ROC-AUC: `{metrics.get('roc_auc', '') if isinstance(metrics, dict) else ''}`",
        f"- Best candidate: `{best_candidate.get('id', '')}` rank `{best_candidate.get('rank', '')}` score `{best_candidate.get('docking_score', '')}`",
        "",
        "## Top-K Control Recovery",
        "",
        "| cutoff | controls in top K | control hit rate | decoys in top K | candidates in top K | enrichment |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for key in sorted(top_k, key=lambda item: int(str(item).replace("top_", "") or 0)):
        item = top_k.get(key, {})
        if not isinstance(item, dict):
            continue
        lines.append(
            f"| {item.get('k', '')} | {item.get('control_hit_count', '')} | {item.get('control_hit_rate', '')} | "
            f"{item.get('negative_count', '')} | {item.get('candidate_count', '')} | {item.get('enrichment_factor', '')} |"
        )
    lines.extend(
        [
            "",
            "## Ranking Preview",
            "",
            "| rank | label | panel | id | docking score | pose pass | backend |",
            "|---:|---|---|---|---:|---|---|",
        ]
    )
    for row in ranking[:80]:
        lines.append(
            f"| {row.get('rank', '')} | {row.get('label', '')} | {row.get('panel_type', '')} | "
            f"{row.get('id', '')} | {row.get('docking_score', '')} | {row.get('pose_pass', '')} | {row.get('backend', '')} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "- This is a retrospective computational benchmark against known controls and decoys, not wet-lab validation.",
            "- A high ROC-AUC means the current docking setup ranks configured positive controls above configured decoys in this run.",
            "- Candidate ranks are prioritization signals only and must not be read as measured potency, antiviral efficacy, safety, or clinical value.",
        ]
    )
    return "\n".join(lines) + "\n"


def stage46_retrospective_benchmark(args: argparse.Namespace) -> None:
    project = project_path(args.project)
    ensure_project_dirs(project)
    round_no = int(getattr(args, "round", 1) or 1)
    positive_types = stage46_parse_types(getattr(args, "positive_types", ""), STAGE46_DEFAULT_POSITIVE_TYPES)
    negative_types = stage46_parse_types(getattr(args, "negative_types", ""), STAGE46_DEFAULT_NEGATIVE_TYPES)
    top_k = stage46_parse_top_k(getattr(args, "top_k", ""))

    source_scores = project / "stage4_5" / f"round_{round_no}_control_docking_scores.csv"
    if not source_scores.exists():
        raise SystemExit(f"Stage 4.6 requires Stage 4.5 docking scores first: {source_scores}")
    source_rows = read_csv(source_scores)
    ranking = stage46_rank_rows(source_rows, positive_types, negative_types)
    if not ranking:
        raise SystemExit(f"Stage 4.6 found no scored rows in: {source_scores}")

    positive_scores = [float(row["docking_score"]) for row in ranking if row.get("label") == "positive"]
    negative_scores = [float(row["docking_score"]) for row in ranking if row.get("label") == "negative"]
    candidate_count = sum(1 for row in ranking if row.get("label") == "candidate")
    ignored_count = sum(1 for row in ranking if row.get("label") == "ignored")
    panel_counts = Counter(str(row.get("panel_type", "")) for row in ranking)
    stage46_dir = project / "stage4_6"
    ranking_path = stage46_dir / f"round_{round_no}_retrospective_ranking.csv"
    benchmark_path = stage46_dir / f"round_{round_no}_retrospective_benchmark.json"
    report_path = project / "reports" / f"stage4_6_round_{round_no}_retrospective_benchmark.md"

    benchmark = {
        "schema_version": "0.1",
        "stage": 4.6,
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "round": round_no,
        "source_scores": str(source_scores),
        "positive_types": positive_types,
        "negative_types": negative_types,
        "counts": {
            "source_rows": len(source_rows),
            "scored_rows": len(ranking),
            "positives": len(positive_scores),
            "negatives": len(negative_scores),
            "candidates": candidate_count,
            "ignored": ignored_count,
            "by_panel_type": dict(panel_counts),
        },
        "metrics": {
            "roc_auc": stage46_auc(positive_scores, negative_scores),
            "positive_score_summary": stage46_score_summary(positive_scores),
            "negative_score_summary": stage46_score_summary(negative_scores),
            "top_k": stage46_top_k_metrics(ranking, top_k, len(positive_scores)),
            "best_candidate": stage46_best_candidate(ranking),
        },
        "ranking_preview": ranking[:20],
        "files": {
            "benchmark": str(benchmark_path),
            "ranking": str(ranking_path),
            "source_scores": str(source_scores),
            "report": str(report_path),
        },
        "boundary": [
            "Stage 4.6 is a retrospective computational benchmark, not wet-lab validation.",
            "ROC-AUC and Top-K recovery only measure whether configured controls outrank configured decoys under the same docking setup.",
            "Candidate ranks remain prioritization evidence and do not prove efficacy, potency, toxicity, selectivity, or clinical value.",
        ],
    }
    write_csv(ranking_path, ranking, STAGE46_RANKING_FIELDS)
    write_json(benchmark_path, benchmark)
    write_text(report_path, render_stage46_report(benchmark, ranking))
    print(f"Wrote Stage 4.6 retrospective ranking: {ranking_path}")
    print(f"Wrote Stage 4.6 retrospective benchmark: {benchmark_path}")
    print(f"Wrote Stage 4.6 report: {report_path}")


def stage4_real(args: argparse.Namespace) -> None:
    project = project_path(args.project)
    ensure_project_dirs(project)
    config = load_config(project)
    target = config.get("target", {})
    if not isinstance(target, dict):
        target = {}
    cli_center = parse_triplet(getattr(args, "pocket_center", "")) if str(getattr(args, "pocket_center", "") or "").strip() else []
    cli_size = parse_triplet(getattr(args, "pocket_size", "")) if str(getattr(args, "pocket_size", "") or "").strip() else []
    if cli_center or cli_size or str(getattr(args, "pocket_source", "") or "").strip():
        pocket = target.get("pocket", {})
        if not isinstance(pocket, dict):
            pocket = {}
        if cli_center:
            pocket["center"] = cli_center
        if cli_size:
            pocket["size"] = cli_size
        if str(getattr(args, "pocket_source", "") or "").strip():
            pocket["source"] = str(getattr(args, "pocket_source"))
        target["pocket"] = pocket
        config["target"] = target
        write_json(project / "config.json", config)
    target_id = args.target or str(target.get("target_catalog_id", "") or "")
    libs = require_rdkit()
    capabilities = stage4_capabilities()
    rdkit_status = capabilities.get("rdkit", {})
    if isinstance(rdkit_status, dict):
        rdkit_status["status"] = "available"
        rdkit_status["version"] = str(getattr(libs["rdBase"], "rdkitVersion", rdkit_status.get("version", "")))

    candidate_rows = load_stage4_candidates(project, args.round, args.input_csv)
    controls = load_stage4_controls(target_id, args.controls_csv)
    descriptor_rows = [rdkit_stage4_descriptor(row, args.round, libs) for row in candidate_rows]
    similarity_rows, best_by_candidate = stage4_similarity_rows(descriptor_rows, controls, libs)
    diversity_rows = stage4_diversity_rows(descriptor_rows, best_by_candidate, libs, args.top)

    stage4_dir = project / "stage4"
    descriptor_path = stage4_dir / f"round_{args.round}_real_descriptors.csv"
    similarity_path = stage4_dir / f"round_{args.round}_similarity_to_controls.csv"
    diversity_path = stage4_dir / f"round_{args.round}_diverse_selection.csv"
    decoy_path = stage4_dir / f"round_{args.round}_decoys.csv"
    benchmark_path = stage4_dir / f"round_{args.round}_benchmark_panel.csv"
    benchmark_sdf_path = stage4_dir / f"round_{args.round}_benchmark_panel.sdf"
    receptor_package_path = stage4_dir / f"round_{args.round}_receptor_package.json"
    docking_inputs_path = stage4_dir / f"round_{args.round}_docking_inputs.csv"
    docking_plan_path = stage4_dir / f"round_{args.round}_docking_plan.json"
    docking_scores_path = stage4_dir / f"round_{args.round}_docking_scores_template.csv"
    validation_metrics_path = stage4_dir / f"round_{args.round}_validation_metrics.json"
    images_dir = stage4_dir / f"round_{args.round}_2d"
    sdf_path = stage4_dir / f"round_{args.round}_ligands.sdf"
    assets_path = stage4_dir / f"round_{args.round}_stage4_assets.json"
    report_path = project / "reports" / f"stage4_round_{args.round}_report.md"

    receptor_package = write_stage4_receptor_package(
        receptor_package_path,
        project,
        args.round,
        target_id,
        config,
        args,
    )
    decoy_rows = build_stage4_decoys(int(getattr(args, "decoys", 0) or 0), libs)
    benchmark_rows = stage4_benchmark_panel_rows(
        descriptor_rows,
        controls,
        decoy_rows,
        best_by_candidate,
        target_id,
        libs,
    )
    image_rows = render_stage4_2d_images(
        images_dir,
        descriptor_rows,
        controls,
        decoy_rows,
        libs,
        min(args.top, 12),
        bool(getattr(args, "render_2d", True)),
    )
    write_csv(descriptor_path, descriptor_rows, STAGE4_DESCRIPTOR_FIELDS)
    write_csv(similarity_path, similarity_rows, STAGE4_SIMILARITY_FIELDS)
    write_csv(diversity_path, diversity_rows, STAGE4_DIVERSITY_FIELDS)
    write_stage4_decoys(decoy_path, decoy_rows)
    write_csv(benchmark_path, benchmark_rows, STAGE4_BENCHMARK_FIELDS)
    benchmark_sdf_written = write_stage4_benchmark_sdf(benchmark_sdf_path, benchmark_rows, libs, args.seed)
    sdf_written = 0
    sdf_warnings: List[Dict[str, str]] = []
    if not args.no_sdf:
        sdf_written, sdf_warnings = write_stage4_sdf(sdf_path, descriptor_rows, libs, args.seed, args.max_conformers)
    else:
        write_text(sdf_path, "")
    run_docking_requested = bool(getattr(args, "run_docking", False))
    docking_timeout = int(getattr(args, "docking_timeout", 600) or 600)
    docking_box = stage4_box_from_receptor_package(receptor_package)
    receptor_preparation = stage4_prepare_receptor_pdbqt(receptor_package, capabilities, docking_box, docking_timeout) if run_docking_requested else {}
    if receptor_preparation.get("status") == "ready" and receptor_preparation.get("receptor_pdbqt"):
        receptor_package["local_receptor_pdbqt"] = receptor_preparation["receptor_pdbqt"]
    ligand_records, ligand_preparation_attempts = stage4_prepare_ligand_records(
        benchmark_rows,
        libs,
        capabilities,
        stage4_dir / f"round_{args.round}_prepared_ligands",
        args.seed,
        docking_timeout,
        run_docking_requested,
    )
    docking_input_rows = stage4_docking_inputs(benchmark_rows, receptor_package, sdf_path, ligand_records)
    write_csv(docking_inputs_path, docking_input_rows, STAGE4_DOCKING_INPUT_FIELDS)
    write_csv(docking_scores_path, [], STAGE4_DOCKING_SCORE_FIELDS)
    docking_plan = write_stage4_docking_plan(
        docking_plan_path,
        str(getattr(args, "docking_backend", "auto") or "auto"),
        bool(getattr(args, "run_docking", False)),
        capabilities,
        receptor_package,
        docking_inputs_path,
        docking_scores_path,
        docking_box,
        receptor_preparation,
        sum(1 for item in ligand_records.values() if item.get("ligand_pdbqt")),
    )
    docking_plan, docking_results = run_stage4_docking_if_possible(
        docking_plan,
        docking_input_rows,
        descriptor_rows,
        docking_scores_path,
        docking_timeout,
    )
    docking_plan, imported_docking_results = import_stage4_external_docking_scores(
        getattr(args, "external_scores", None),
        docking_scores_path,
        docking_plan,
    )
    if imported_docking_results:
        docking_results = imported_docking_results
    docking_plan["ligand_preparation_attempts"] = ligand_preparation_attempts[:50]
    write_json(docking_plan_path, docking_plan)
    write_json(receptor_package_path, receptor_package)
    validation_metrics = write_stage4_validation_metrics(
        validation_metrics_path,
        benchmark_rows,
        similarity_rows,
        docking_plan,
    )

    files = {
        "real_descriptors": str(descriptor_path),
        "similarity_to_controls": str(similarity_path),
        "diverse_selection": str(diversity_path),
        "decoys": str(decoy_path),
        "benchmark_panel": str(benchmark_path),
        "benchmark_sdf": str(benchmark_sdf_path),
        "receptor_package": str(receptor_package_path),
        "docking_inputs": str(docking_inputs_path),
        "docking_plan": str(docking_plan_path),
        "docking_scores_template": str(docking_scores_path),
        "validation_metrics": str(validation_metrics_path),
        "molecule_2d_dir": str(images_dir),
        "ligands_sdf": str(sdf_path),
        "report": str(report_path),
    }
    valid_count = sum(1 for row in descriptor_rows if str(row.get("valid", "")) == "1")
    write_json(
        assets_path,
        {
            "schema_version": "0.1",
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "stage": 4,
            "round": args.round,
            "target_id": target_id,
            "candidate_count": len(descriptor_rows),
            "valid_count": valid_count,
            "invalid_count": len(descriptor_rows) - valid_count,
            "controls_count": len(controls),
            "decoy_count": len(decoy_rows),
            "benchmark_panel_count": len(benchmark_rows),
            "benchmark_sdf_molecules_written": benchmark_sdf_written,
            "sdf_ligands_written": sdf_written,
            "rendered_2d_images": image_rows,
            "receptor_package": receptor_package,
            "validation_metrics": validation_metrics,
            "docking_results_count": len(docking_results),
            "real_libraries": {
                "rdkit": capabilities.get("rdkit", {}),
                "numpy": capabilities.get("modules", {}).get("numpy", {}) if isinstance(capabilities.get("modules"), dict) else {},
                "pandas": capabilities.get("modules", {}).get("pandas", {}) if isinstance(capabilities.get("modules"), dict) else {},
            },
            "optional_libraries": {
                "openbabel": capabilities.get("modules", {}).get("openbabel", {}) if isinstance(capabilities.get("modules"), dict) else {},
                "meeko": capabilities.get("modules", {}).get("meeko", {}) if isinstance(capabilities.get("modules"), dict) else {},
                "vina": capabilities.get("modules", {}).get("vina", {}) if isinstance(capabilities.get("modules"), dict) else {},
                "posebusters": capabilities.get("modules", {}).get("posebusters", {}) if isinstance(capabilities.get("modules"), dict) else {},
            },
            "docking_backend": capabilities.get("docking_backend", {}),
            "docking_plan": docking_plan,
            "ligand_preparation": capabilities.get("ligand_preparation", {}),
            "files": files,
            "rescore": {
                "requested": bool(getattr(args, "rescore", False)),
                "score_file": str(project / "scores" / f"round_{args.round}_scores.csv") if getattr(args, "rescore", False) else "",
                "ranked_file": str(project / "ranked" / f"round_{args.round}_ranked.csv") if getattr(args, "rescore", False) else "",
                "feedback_file": str(project / "feedback" / f"round_{args.round}_feedback.json") if getattr(args, "rescore", False) else "",
            },
            "sdf_warnings": sdf_warnings,
            "boundary": [
                "RDKit descriptors, fingerprints, scaffolds, and SDF export are real chemistry-library operations.",
                "Similarity and descriptor values are computational screening features, not measured antiviral efficacy.",
                "Docking remains unavailable until Vina/GNINA plus ligand/receptor preparation tools are installed and executed.",
            ],
        },
    )
    write_text(
        report_path,
        render_stage4_report(
            args.round,
            target_id,
            descriptor_rows,
            controls,
            similarity_rows,
            diversity_rows,
            capabilities,
            files,
            sdf_written,
            sdf_warnings,
            receptor_package,
            docking_plan,
            validation_metrics,
            len(image_rows),
        ),
    )

    print(f"Wrote Stage 4 receptor package: {receptor_package_path}")
    print(f"Wrote Stage 4 real descriptors: {descriptor_path}")
    print(f"Wrote Stage 4 control similarity: {similarity_path}")
    print(f"Wrote Stage 4 diverse selection: {diversity_path}")
    print(f"Wrote Stage 4 decoys: {decoy_path}")
    print(f"Wrote Stage 4 benchmark panel: {benchmark_path}")
    print(f"Wrote Stage 4 docking plan: {docking_plan_path}")
    print(f"Wrote Stage 4 validation metrics: {validation_metrics_path}")
    print(f"Wrote Stage 4 ligand SDF: {sdf_path}")
    print(f"Wrote Stage 4 assets: {assets_path}")
    print(f"Wrote Stage 4 report: {report_path}")

    if getattr(args, "rescore", False):
        rank_top = int(getattr(args, "rank_top", 0) or args.top)
        feedback_top = int(getattr(args, "feedback_top", 0) or rank_top)
        rescore_external_scores = getattr(args, "external_scores", None)
        if not rescore_external_scores and any(row.get("docking_score") not in {"", None} for row in docking_results):
            rescore_external_scores = str(docking_scores_path)
        score_candidates(
            argparse.Namespace(
                project=str(project),
                round=args.round,
                external_scores=rescore_external_scores,
                real_descriptors=str(descriptor_path),
            )
        )
        rank_candidates(argparse.Namespace(project=str(project), round=args.round, top=rank_top))
        feedback(argparse.Namespace(project=str(project), round=args.round, top=feedback_top))
        print("Stage 4 rescore completed: scores -> ranked -> feedback")


def stage5_read_json(path: Path) -> Dict[str, object]:
    if not path.exists():
        return {}
    try:
        data = read_json(path)
    except (json.JSONDecodeError, OSError):
        return {}
    return data if isinstance(data, dict) else {}


def stage5_to_int(value: object, default: int = 0) -> int:
    if value in {"", None}:
        return default
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def stage5_to_float(value: object, default: object = "") -> object:
    if value in {"", None}:
        return default
    try:
        return round(float(str(value)), 4)
    except (TypeError, ValueError):
        return default


def stage5_relative_path(path: Path, base: Path) -> str:
    try:
        return os.path.relpath(path.resolve(), base.resolve())
    except OSError:
        return str(path)


def stage5_file_record(path: Path, base: Path) -> Dict[str, object]:
    return {
        "path": str(path),
        "relative_path": stage5_relative_path(path, base),
        "exists": path.exists(),
    }


def stage5_first_stage2_target(stage2_rows: Sequence[Dict[str, str]], stage2_assets: Dict[str, object]) -> Dict[str, object]:
    if stage2_rows:
        return dict(stage2_rows[0])
    targets = stage2_assets.get("targets", [])
    if isinstance(targets, list) and targets and isinstance(targets[0], dict):
        item = targets[0]
        quality = item.get("evidence_quality", {})
        if not isinstance(quality, dict):
            quality = {}
        return {
            "target_id": item.get("target_id", ""),
            "target": item.get("target", ""),
            "evidence_score": quality.get("score", ""),
            "readiness": quality.get("readiness", ""),
        }
    return {}


def stage5_target_summary(
    config: Dict[str, object],
    brief: Dict[str, object],
    stage2_row: Dict[str, object],
    stage4_assets: Dict[str, object],
) -> Dict[str, object]:
    target_config = config.get("target", {})
    if not isinstance(target_config, dict):
        target_config = {}
    site = brief.get("binding_site", {})
    if not isinstance(site, dict):
        site = {}
    target_id = (
        stage4_assets.get("target_id")
        or brief.get("target_catalog_id")
        or target_config.get("target_catalog_id")
        or stage2_row.get("target_id")
        or ""
    )
    return {
        "id": target_id,
        "name": brief.get("target_name") or target_config.get("name") or stage2_row.get("target", ""),
        "disease_context": brief.get("disease_context", ""),
        "evidence_score": stage5_to_float(stage2_row.get("evidence_score", "")),
        "evidence_readiness": stage2_row.get("readiness", ""),
        "primary_pdb": stage2_row.get("primary_pdb") or target_config.get("pdb_id", ""),
        "reference_ligand": site.get("reference_ligand") or stage2_row.get("reference_ligand", ""),
        "key_residues": site.get("key_residues", []),
        "positive_controls": stage2_row.get("positive_controls", ""),
    }


def stage5_count_advanced(ranked_rows: Sequence[Dict[str, str]]) -> int:
    return sum(1 for row in ranked_rows if str(row.get("decision", "")).lower() == "advance")


def stage5_count_filtered(stage3_assets: Dict[str, object], filtered_rows: Sequence[Dict[str, str]], candidates: Sequence[Dict[str, str]]) -> int:
    if "passed_filter" in stage3_assets:
        return stage5_to_int(stage3_assets.get("passed_filter"))
    passed = [row for row in filtered_rows if str(row.get("passed_filter", "")).lower() in {"1", "true", "yes", "pass", "passed"}]
    if passed:
        return len(passed)
    return len(candidates)


def stage5_panel_count(panel_counts: Dict[str, object], benchmark_rows: Sequence[Dict[str, str]], names: Sequence[str]) -> int:
    count = 0
    for name in names:
        count += stage5_to_int(panel_counts.get(name, 0))
    if count:
        return count
    wanted = set(names)
    return sum(1 for row in benchmark_rows if row.get("panel_type") in wanted)


def stage5_table(rows: Sequence[Dict[str, str]], limit: int, preferred_fields: Sequence[str]) -> List[Dict[str, str]]:
    output = []
    for row in rows[:limit]:
        selected = {}
        for field in preferred_fields:
            if field in row:
                selected[field] = row.get(field, "")
        for key, value in row.items():
            if key not in selected and len(selected) < len(preferred_fields):
                selected[key] = value
        output.append(selected)
    return output


def stage5_molecule_images(project: Path, stage5_dir: Path, round_no: int, stage4_assets: Dict[str, object]) -> List[Dict[str, str]]:
    images: List[Dict[str, str]] = []
    seen = set()
    rendered = stage4_assets.get("rendered_2d_images", [])
    if isinstance(rendered, list):
        for item in rendered:
            if not isinstance(item, dict):
                continue
            raw_path = item.get("path", "")
            if not raw_path:
                continue
            path = Path(str(raw_path)).expanduser()
            if not path.is_absolute():
                path = project / path
            if path in seen:
                continue
            seen.add(path)
            images.append(
                {
                    "panel_type": str(item.get("panel_type", "")),
                    "id": str(item.get("id", "")),
                    "path": str(path),
                    "relative_path": stage5_relative_path(path, stage5_dir),
                    "exists": str(path.exists()).lower(),
                }
            )
    image_dir = project / "stage4" / f"round_{round_no}_2d"
    if image_dir.exists():
        for path in sorted(image_dir.glob("*.png")):
            if path in seen:
                continue
            name = path.stem
            panel_type = name.split("_", 1)[0] if "_" in name else "molecule"
            item_id = name.split("_", 1)[1] if "_" in name else name
            images.append(
                {
                    "panel_type": panel_type,
                    "id": item_id,
                    "path": str(path),
                    "relative_path": stage5_relative_path(path, stage5_dir),
                    "exists": "true",
                }
            )
    return images


def build_stage5_dashboard_data(project: Path, round_no: int, title: str) -> Dict[str, object]:
    stage5_dir = project / "stage5"
    config = stage5_read_json(project / "config.json")
    brief = stage5_read_json(project / "briefs" / "target_brief.json")
    stage2_rows = read_csv(project / "evidence" / "stage2_target_sources.csv")
    stage2_assets = stage5_read_json(project / "evidence" / "stage2_closed_loop_assets.json")
    stage2_row = stage5_first_stage2_target(stage2_rows, stage2_assets)
    stage3_assets = stage5_read_json(project / "stage3" / f"round_{round_no}_stage3_assets.json")
    stage4_assets = stage5_read_json(project / "stage4" / f"round_{round_no}_stage4_assets.json")
    candidates = read_csv(project / "candidates" / f"round_{round_no}_candidates.csv")
    filtered_rows = read_csv(project / "filtered" / f"round_{round_no}_filtered.csv")
    scores = read_csv(project / "scores" / f"round_{round_no}_scores.csv")
    ranked = read_csv(project / "ranked" / f"round_{round_no}_ranked.csv")
    seeds = read_csv(project / "seeds" / f"round_{round_no}_seeds.csv")
    feedback_doc = stage5_read_json(project / "feedback" / f"round_{round_no}_feedback.json")
    descriptors = read_csv(project / "stage4" / f"round_{round_no}_real_descriptors.csv")
    similarity = read_csv(project / "stage4" / f"round_{round_no}_similarity_to_controls.csv")
    benchmark = read_csv(project / "stage4" / f"round_{round_no}_benchmark_panel.csv")
    docking_plan = stage5_read_json(project / "stage4" / f"round_{round_no}_docking_plan.json")
    validation = stage5_read_json(project / "stage4" / f"round_{round_no}_validation_metrics.json")

    ranked.sort(key=lambda row: stage5_to_int(row.get("rank"), 999999))
    panel_counts = validation.get("panel_counts", {})
    if not isinstance(panel_counts, dict):
        panel_counts = {}
    validation_readiness = validation.get("readiness", {})
    if not isinstance(validation_readiness, dict):
        validation_readiness = {}
    valid_rdkit = stage5_to_int(stage4_assets.get("valid_count"))
    if not valid_rdkit:
        valid_rdkit = sum(1 for row in descriptors if str(row.get("valid", "")) == "1")
    benchmark_count = len(benchmark)
    advanced = stage5_count_advanced(ranked)
    molecule_images = stage5_molecule_images(project, stage5_dir, round_no, stage4_assets)
    raw_candidates = stage5_to_int(stage3_assets.get("raw_candidates"), len(candidates))
    filtered_candidates = stage5_count_filtered(stage3_assets, filtered_rows, candidates)
    docking_status = str(
        validation_readiness.get("docking")
        or validation.get("docking_status", "")
        or docking_plan.get("status", "")
        or nested_get(stage4_assets, ["docking_plan", "status"])
        or "missing"
    )

    files = {
        "config": stage5_file_record(project / "config.json", stage5_dir),
        "target_brief": stage5_file_record(project / "briefs" / "target_brief.json", stage5_dir),
        "stage2_matrix": stage5_file_record(project / "evidence" / "stage2_target_sources.csv", stage5_dir),
        "stage2_assets": stage5_file_record(project / "evidence" / "stage2_closed_loop_assets.json", stage5_dir),
        "stage3_assets": stage5_file_record(project / "stage3" / f"round_{round_no}_stage3_assets.json", stage5_dir),
        "filtered_candidates": stage5_file_record(project / "filtered" / f"round_{round_no}_filtered.csv", stage5_dir),
        "candidates": stage5_file_record(project / "candidates" / f"round_{round_no}_candidates.csv", stage5_dir),
        "scores": stage5_file_record(project / "scores" / f"round_{round_no}_scores.csv", stage5_dir),
        "ranked": stage5_file_record(project / "ranked" / f"round_{round_no}_ranked.csv", stage5_dir),
        "feedback": stage5_file_record(project / "feedback" / f"round_{round_no}_feedback.json", stage5_dir),
        "seeds": stage5_file_record(project / "seeds" / f"round_{round_no}_seeds.csv", stage5_dir),
        "stage4_assets": stage5_file_record(project / "stage4" / f"round_{round_no}_stage4_assets.json", stage5_dir),
        "real_descriptors": stage5_file_record(project / "stage4" / f"round_{round_no}_real_descriptors.csv", stage5_dir),
        "benchmark_panel": stage5_file_record(project / "stage4" / f"round_{round_no}_benchmark_panel.csv", stage5_dir),
        "docking_plan": stage5_file_record(project / "stage4" / f"round_{round_no}_docking_plan.json", stage5_dir),
        "validation_metrics": stage5_file_record(project / "stage4" / f"round_{round_no}_validation_metrics.json", stage5_dir),
    }

    data = {
        "schema_version": "0.1",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "stage": 5,
        "round": round_no,
        "title": title,
        "project": {
            "name": project.name,
            "path": str(project),
        },
        "target": stage5_target_summary(config, brief, stage2_row, stage4_assets),
        "readiness": {
            "target_evidence": stage2_row.get("readiness", "missing"),
            "candidate_intake": "ready" if candidates else "missing",
            "rdkit_validation": validation_readiness.get("rdkit_descriptors") or ("ready" if valid_rdkit else "missing"),
            "control_panel": validation_readiness.get("control_panel") or ("ready" if stage5_panel_count(panel_counts, benchmark, ["positive_control", "reference_control"]) else "missing"),
            "decoy_panel": validation_readiness.get("decoy_panel") or ("ready" if stage5_panel_count(panel_counts, benchmark, ["decoy"]) else "missing"),
            "docking": docking_status,
            "feedback": "ready" if feedback_doc or seeds else "missing",
        },
        "metrics": {
            "raw_candidates": raw_candidates,
            "filtered_candidates": filtered_candidates,
            "candidate_rows": len(candidates),
            "scored_candidates": len(scores),
            "ranked_candidates": len(ranked),
            "advanced": advanced,
            "valid_rdkit": valid_rdkit,
            "invalid_rdkit": max(len(descriptors) - valid_rdkit, 0),
            "controls": stage5_panel_count(panel_counts, benchmark, ["positive_control", "reference_control"]),
            "decoys": stage5_panel_count(panel_counts, benchmark, ["decoy"]),
            "benchmark_panel": benchmark_count,
            "molecule_images": len(molecule_images),
        },
        "files": files,
        "tables": {
            "stage2_targets": stage5_table(
                stage2_rows,
                10,
                ["target_id", "target", "evidence_score", "readiness", "primary_pdb", "positive_controls", "pubmed_articles"],
            ),
            "ranked_top": stage5_table(
                ranked,
                20,
                ["rank", "id", "smiles", "total_proxy", "docking_proxy", "pose_proxy", "score_source", "decision"],
            ),
            "real_descriptors": stage5_table(
                descriptors,
                20,
                ["id", "canonical_smiles", "valid", "qed", "mw", "logp", "tpsa", "lipinski_violations", "murcko_scaffold"],
            ),
            "benchmark_panel": stage5_table(
                benchmark,
                50,
                ["panel_type", "id", "canonical_smiles", "qed", "mw", "nearest_control", "max_similarity_to_controls"],
            ),
            "control_similarity_top": stage5_table(
                sorted(similarity, key=lambda row: stage5_to_float(row.get("tanimoto"), 0.0), reverse=True),
                20,
                ["candidate_id", "control_drug", "tanimoto", "control_role", "fingerprint"],
            ),
            "seeds": stage5_table(seeds, 20, ["id", "smiles", "note"]),
            "validation_metrics": validation,
            "docking_plan": docking_plan,
            "feedback": feedback_doc,
        },
        "molecule_images": molecule_images,
        "boundary": [
            "Computational screening dashboard only.",
            "RDKit descriptors and fingerprint similarity are real library calculations, not measured antiviral efficacy.",
            "Docking is displayed only when a backend is installed and actually run; skipped or not_available means no docking score was produced.",
            "Known drugs are controls and benchmarks, not newly discovered candidates.",
        ],
    }
    return data


def stage5_embedded_json(data: Dict[str, object]) -> str:
    return json.dumps(data, ensure_ascii=False, indent=2).replace("</", "<\\/")


def stage5_html_escape(value: object) -> str:
    text = str(value if value is not None else "")
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#39;")
    )


def render_stage5_html(data: Dict[str, object]) -> str:
    title = str(data.get("title") or "AI Molecule Closed Loop Dashboard")
    safe_title = stage5_html_escape(title)
    embedded = stage5_embedded_json(data)
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{safe_title}</title>
  <link rel="stylesheet" href="styles.css">
</head>
<body>
  <header class="topbar">
    <div>
      <p class="eyebrow">Stage 5 Product Dashboard</p>
      <h1>{safe_title}</h1>
      <p class="subline" id="project-line"></p>
    </div>
    <div class="status-pill" id="docking-pill"></div>
  </header>

  <main>
    <section class="section" aria-labelledby="overview-title">
      <div class="section-head">
        <h2 id="overview-title">闭环概览</h2>
        <p>Computational screening only. 不代表真实药效、剂量、安全性或临床有效性。</p>
      </div>
      <div class="metric-grid" id="metric-grid"></div>
      <div class="readiness-strip" id="readiness-strip"></div>
    </section>

    <section class="section two-col" aria-labelledby="target-title">
      <div>
        <h2 id="target-title">靶点与证据</h2>
        <dl class="kv" id="target-kv"></dl>
      </div>
      <div>
        <h2>验证通路</h2>
        <div id="pipeline"></div>
      </div>
    </section>

    <section class="section" aria-labelledby="ranked-title">
      <div class="section-head">
        <h2 id="ranked-title">Top 候选分子</h2>
      </div>
      <div class="table-wrap">
        <table id="ranked-table"></table>
      </div>
    </section>

    <section class="section" aria-labelledby="benchmark-title">
      <div class="section-head">
        <h2 id="benchmark-title">候选 / 控药 / Decoy 面板</h2>
      </div>
      <div class="table-wrap">
        <table id="benchmark-table"></table>
      </div>
    </section>

    <section class="section" aria-labelledby="gallery-title">
      <div class="section-head">
        <h2 id="gallery-title">分子图谱</h2>
      </div>
      <div class="gallery" id="gallery"></div>
    </section>

    <section class="section two-col" aria-labelledby="docking-title">
      <div>
        <h2 id="docking-title">Docking 状态</h2>
        <dl class="kv" id="docking-kv"></dl>
      </div>
      <div>
        <h2>反馈种子</h2>
        <div class="table-wrap compact">
          <table id="seeds-table"></table>
        </div>
      </div>
    </section>

    <section class="section" aria-labelledby="files-title">
      <div class="section-head">
        <h2 id="files-title">资产文件</h2>
      </div>
      <div class="file-grid" id="file-grid"></div>
    </section>
  </main>

  <script id="dashboard-data" type="application/json">{embedded}</script>
  <script src="app.js"></script>
  <script>window.__dashboardDataPath = "dashboard_data.json";</script>
</body>
</html>
"""


def render_stage5_css() -> str:
    return """* {
  box-sizing: border-box;
}

:root {
  color-scheme: light;
  --bg: #f7f8fa;
  --text: #1f2933;
  --muted: #657181;
  --line: #d8dee7;
  --panel: #ffffff;
  --green: #1f7a55;
  --amber: #b7791f;
  --red: #b42318;
  --blue: #2563a8;
}

body {
  margin: 0;
  background: var(--bg);
  color: var(--text);
  font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  letter-spacing: 0;
}

.topbar {
  display: flex;
  align-items: end;
  justify-content: space-between;
  gap: 24px;
  padding: 28px 40px 22px;
  border-bottom: 1px solid var(--line);
  background: #ffffff;
}

.eyebrow {
  margin: 0 0 6px;
  color: var(--blue);
  font-size: 12px;
  font-weight: 700;
  text-transform: uppercase;
}

h1, h2, p {
  margin: 0;
}

h1 {
  font-size: 28px;
  line-height: 1.2;
  font-weight: 720;
}

h2 {
  font-size: 17px;
  line-height: 1.25;
  margin-bottom: 14px;
}

.subline {
  color: var(--muted);
  font-size: 13px;
  margin-top: 8px;
}

.status-pill {
  min-width: 142px;
  padding: 8px 12px;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: #f9fafb;
  text-align: center;
  font-size: 13px;
  font-weight: 650;
}

main {
  max-width: 1320px;
  margin: 0 auto;
  padding: 24px 28px 48px;
}

.section {
  padding: 24px 0;
  border-bottom: 1px solid var(--line);
}

.section-head {
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  gap: 18px;
  margin-bottom: 16px;
}

.section-head p {
  color: var(--muted);
  font-size: 13px;
}

.two-col {
  display: grid;
  grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
  gap: 32px;
}

.metric-grid {
  display: grid;
  grid-template-columns: repeat(6, minmax(0, 1fr));
  gap: 10px;
}

.metric {
  min-height: 84px;
  padding: 14px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: var(--panel);
}

.metric strong {
  display: block;
  font-size: 24px;
  line-height: 1.15;
}

.metric span {
  display: block;
  color: var(--muted);
  font-size: 12px;
  margin-top: 8px;
}

.readiness-strip {
  display: flex;
  flex-wrap: wrap;
  gap: 8px;
  margin-top: 14px;
}

.badge {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  min-height: 30px;
  padding: 6px 10px;
  border: 1px solid var(--line);
  border-radius: 999px;
  background: #fff;
  color: var(--muted);
  font-size: 12px;
}

.badge.ready {
  color: var(--green);
  border-color: #badbcc;
  background: #f0fdf7;
}

.badge.warn {
  color: var(--amber);
  border-color: #f0d399;
  background: #fff9eb;
}

.badge.missing {
  color: var(--red);
  border-color: #f3c4bf;
  background: #fff5f5;
}

.kv {
  display: grid;
  grid-template-columns: minmax(120px, 180px) minmax(0, 1fr);
  gap: 10px 14px;
  margin: 0;
}

.kv dt {
  color: var(--muted);
  font-size: 13px;
}

.kv dd {
  margin: 0;
  min-width: 0;
  overflow-wrap: anywhere;
  font-size: 13px;
}

.pipeline {
  display: grid;
  gap: 8px;
}

.step {
  display: grid;
  grid-template-columns: 28px minmax(0, 1fr) auto;
  align-items: center;
  gap: 10px;
  padding: 10px 12px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fff;
  font-size: 13px;
}

.step-index {
  display: inline-grid;
  place-items: center;
  width: 24px;
  height: 24px;
  border-radius: 50%;
  background: #edf2f7;
  color: #394150;
  font-weight: 700;
}

.table-wrap {
  width: 100%;
  overflow-x: auto;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fff;
}

.table-wrap.compact {
  max-height: 320px;
}

table {
  width: 100%;
  border-collapse: collapse;
  font-size: 12px;
}

th, td {
  padding: 10px 12px;
  border-bottom: 1px solid var(--line);
  text-align: left;
  vertical-align: top;
}

th {
  color: #4a5563;
  background: #f4f6f8;
  font-weight: 700;
}

td {
  overflow-wrap: anywhere;
}

.gallery {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
  gap: 12px;
}

.mol {
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fff;
  overflow: hidden;
}

.mol img {
  display: block;
  width: 100%;
  aspect-ratio: 18 / 13;
  object-fit: contain;
  background: #ffffff;
}

.mol footer {
  display: flex;
  justify-content: space-between;
  gap: 8px;
  padding: 9px 10px;
  border-top: 1px solid var(--line);
  font-size: 12px;
}

.file-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(260px, 1fr));
  gap: 10px;
}

.file-link {
  display: block;
  min-height: 54px;
  padding: 11px 12px;
  border: 1px solid var(--line);
  border-radius: 8px;
  background: #fff;
  color: var(--text);
  text-decoration: none;
  overflow-wrap: anywhere;
  font-size: 12px;
}

.file-link strong {
  display: block;
  margin-bottom: 5px;
  font-size: 13px;
}

.empty {
  padding: 18px;
  color: var(--muted);
  border: 1px dashed var(--line);
  border-radius: 8px;
  background: #fff;
  font-size: 13px;
}

@media (max-width: 980px) {
  .topbar {
    display: block;
    padding: 22px 20px;
  }

  .status-pill {
    display: inline-block;
    margin-top: 14px;
  }

  main {
    padding: 18px 18px 36px;
  }

  .metric-grid {
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }

  .two-col {
    grid-template-columns: 1fr;
  }

  .section-head {
    display: block;
  }

  .section-head p {
    margin-top: 8px;
  }
}

@media (max-width: 560px) {
  h1 {
    font-size: 23px;
  }

  .metric-grid {
    grid-template-columns: 1fr;
  }

  .kv {
    grid-template-columns: 1fr;
  }
}
"""


def render_stage5_js() -> str:
    return """const embedded = document.getElementById("dashboard-data");

function parseData() {
  if (embedded && embedded.textContent.trim()) {
    return Promise.resolve(JSON.parse(embedded.textContent));
  }
  return fetch(window.__dashboardDataPath || "dashboard_data.json").then((response) => response.json());
}

function text(value) {
  if (value === null || value === undefined || value === "") return "-";
  if (Array.isArray(value)) return value.join(", ");
  if (typeof value === "object") return JSON.stringify(value);
  return String(value);
}

function html(value) {
  return text(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function statusClass(value) {
  const raw = text(value).toLowerCase();
  if (raw.includes("ready") || raw.includes("available") || raw.includes("completed")) return "ready";
  if (raw.includes("missing") || raw.includes("failed") || raw.includes("blocked")) return "missing";
  return "warn";
}

function setText(id, value) {
  const node = document.getElementById(id);
  if (node) node.textContent = value;
}

function renderMetrics(metrics) {
  const labels = [
    ["raw_candidates", "原始候选"],
    ["filtered_candidates", "通过过滤"],
    ["valid_rdkit", "RDKit 有效"],
    ["controls", "控药"],
    ["decoys", "Decoy"],
    ["advanced", "反馈进入下一轮"]
  ];
  const root = document.getElementById("metric-grid");
  root.innerHTML = labels.map(([key, label]) => `
    <div class="metric">
      <strong>${html(metrics[key])}</strong>
      <span>${html(label)}</span>
    </div>
  `).join("");
}

function renderReadiness(readiness) {
  const root = document.getElementById("readiness-strip");
  root.innerHTML = Object.entries(readiness || {}).map(([key, value]) => {
    const cls = statusClass(value);
    return `<span class="badge ${cls}"><strong>${html(key)}</strong>${html(value)}</span>`;
  }).join("");
}

function renderKv(id, entries) {
  const root = document.getElementById(id);
  root.innerHTML = entries.map(([key, value]) => `<dt>${html(key)}</dt><dd>${html(value)}</dd>`).join("");
}

function renderPipeline(data) {
  const readiness = data.readiness || {};
  const steps = [
    ["靶点证据", readiness.target_evidence],
    ["候选输入", readiness.candidate_intake],
    ["RDKit 校验", readiness.rdkit_validation],
    ["控药/Decoy 面板", `${readiness.control_panel} / ${readiness.decoy_panel}`],
    ["Docking", readiness.docking],
    ["反馈", readiness.feedback]
  ];
  document.getElementById("pipeline").innerHTML = `<div class="pipeline">${steps.map((step, idx) => `
    <div class="step">
      <span class="step-index">${idx + 1}</span>
      <strong>${html(step[0])}</strong>
      <span class="badge ${statusClass(step[1])}">${html(step[1])}</span>
    </div>
  `).join("")}</div>`;
}

function renderTable(id, rows) {
  const root = document.getElementById(id);
  if (!rows || !rows.length) {
    root.outerHTML = `<div class="empty">No rows</div>`;
    return;
  }
  const headers = Object.keys(rows[0]);
  root.innerHTML = `
    <thead><tr>${headers.map((header) => `<th>${html(header)}</th>`).join("")}</tr></thead>
    <tbody>${rows.map((row) => `
      <tr>${headers.map((header) => `<td>${html(row[header])}</td>`).join("")}</tr>
    `).join("")}</tbody>
  `;
}

function renderGallery(images) {
  const root = document.getElementById("gallery");
  if (!images || !images.length) {
    root.innerHTML = `<div class="empty">No molecule images</div>`;
    return;
  }
  root.innerHTML = images.map((item) => `
    <article class="mol">
      <img src="${html(item.relative_path)}" alt="${html(item.id)}">
      <footer><strong>${html(item.id)}</strong><span>${html(item.panel_type)}</span></footer>
    </article>
  `).join("");
}

function renderFiles(files) {
  const root = document.getElementById("file-grid");
  root.innerHTML = Object.entries(files || {}).map(([key, value]) => {
    const cls = value.exists ? "ready" : "missing";
    const href = value.exists ? html(value.relative_path) : "#";
    return `<a class="file-link" href="${href}">
      <strong>${html(key)} <span class="badge ${cls}">${value.exists ? "present" : "missing"}</span></strong>
      ${html(value.relative_path)}
    </a>`;
  }).join("");
}

parseData().then((data) => {
  document.title = text(data.title);
  setText("project-line", `${text(data.project?.name)} · round ${text(data.round)} · generated ${text(data.generated_at)}`);
  setText("docking-pill", `Docking: ${text(data.readiness?.docking)}`);
  renderMetrics(data.metrics || {});
  renderReadiness(data.readiness || {});
  renderKv("target-kv", [
    ["Target ID", data.target?.id],
    ["Target", data.target?.name],
    ["Disease", data.target?.disease_context],
    ["Evidence score", data.target?.evidence_score],
    ["Readiness", data.target?.evidence_readiness],
    ["PDB", data.target?.primary_pdb],
    ["Reference ligand", data.target?.reference_ligand],
    ["Positive controls", data.target?.positive_controls],
    ["Key residues", data.target?.key_residues]
  ]);
  renderPipeline(data);
  renderTable("ranked-table", data.tables?.ranked_top || []);
  renderTable("benchmark-table", data.tables?.benchmark_panel || []);
  renderGallery(data.molecule_images || []);
  renderKv("docking-kv", [
    ["Status", data.tables?.docking_plan?.status || data.readiness?.docking],
    ["Backend", data.tables?.docking_plan?.selected_backend],
    ["Message", data.tables?.docking_plan?.message],
    ["Expected scores", data.tables?.docking_plan?.expected_scores_csv]
  ]);
  renderTable("seeds-table", data.tables?.seeds || []);
  renderFiles(data.files || {});
}).catch((error) => {
  document.body.innerHTML = `<main><div class="empty">Dashboard data could not be loaded: ${html(error)}</div></main>`;
});
"""


def render_stage5_report(data: Dict[str, object], dashboard_path: Path) -> str:
    metrics = data.get("metrics", {})
    readiness = data.get("readiness", {})
    target = data.get("target", {})
    if not isinstance(metrics, dict):
        metrics = {}
    if not isinstance(readiness, dict):
        readiness = {}
    if not isinstance(target, dict):
        target = {}
    return "\n".join(
        [
            f"# Stage 5 Product Dashboard - Round {data.get('round', '')}",
            "",
            "## Summary",
            "",
            f"- Dashboard: `{dashboard_path}`",
            f"- Target: `{target.get('id', '')}` / {target.get('name', '')}",
            f"- Raw candidates: {metrics.get('raw_candidates', 0)}",
            f"- Filtered candidates: {metrics.get('filtered_candidates', 0)}",
            f"- RDKit-valid candidates: {metrics.get('valid_rdkit', 0)}",
            f"- Benchmark panel size: {metrics.get('benchmark_panel', 0)}",
            f"- Advanced to feedback: {metrics.get('advanced', 0)}",
            "",
            "## Readiness",
            "",
        ]
        + [f"- {key}: `{value}`" for key, value in readiness.items()]
        + [
            "",
            "## Interpretation Boundary",
            "",
            "- Computational screening dashboard only.",
            "- RDKit descriptors, similarity, and benchmark panels support prioritization, not efficacy proof.",
            "- Docking status must be `completed` with real external score files before docking values can be discussed.",
            "- Wet-lab validation remains outside this CLI.",
        ]
    ) + "\n"


def stage5_dashboard(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    ensure_project_dirs(project)
    title = getattr(args, "title", "") or f"{project.name} round {args.round} closed-loop dashboard"
    data = build_stage5_dashboard_data(project, args.round, title)
    stage5_dir = project / "stage5"
    data_path = stage5_dir / "dashboard_data.json"
    html_path = stage5_dir / "index.html"
    css_path = stage5_dir / "styles.css"
    js_path = stage5_dir / "app.js"
    report_path = project / "reports" / "stage5_dashboard_report.md"

    write_json(data_path, data)
    write_text(html_path, render_stage5_html(data))
    write_text(css_path, render_stage5_css())
    write_text(js_path, render_stage5_js())
    write_text(report_path, render_stage5_report(data, html_path))

    print(f"Wrote Stage 5 dashboard data: {data_path}")
    print(f"Wrote Stage 5 dashboard HTML: {html_path}")
    print(f"Wrote Stage 5 CSS: {css_path}")
    print(f"Wrote Stage 5 JS: {js_path}")
    print(f"Wrote Stage 5 report: {report_path}")
    return html_path


def stage6_status_from_readiness(value: object, *, docking: bool = False) -> str:
    raw = str(value or "").lower()
    if raw in {"ready", "available", "completed", "ready_for_mvp_closed_loop"} or "ready" in raw:
        return "pass"
    if docking and raw in {"skipped", "not_available", "planned", "attempted_no_scores", "blocked_missing_receptor"}:
        return "warn"
    if raw in {"", "missing", "failed", "blocked"} or "missing" in raw or "failed" in raw:
        return "fail"
    return "warn"


def stage6_quality_gate_rows(data: Dict[str, object]) -> List[Dict[str, object]]:
    readiness = data.get("readiness", {})
    metrics = data.get("metrics", {})
    target = data.get("target", {})
    tables = data.get("tables", {})
    if not isinstance(readiness, dict):
        readiness = {}
    if not isinstance(metrics, dict):
        metrics = {}
    if not isinstance(target, dict):
        target = {}
    if not isinstance(tables, dict):
        tables = {}
    docking_plan = tables.get("docking_plan", {})
    if not isinstance(docking_plan, dict):
        docking_plan = {}

    controls = stage5_to_int(metrics.get("controls", 0))
    decoys = stage5_to_int(metrics.get("decoys", 0))
    valid_rdkit = stage5_to_int(metrics.get("valid_rdkit", 0))
    advanced = stage5_to_int(metrics.get("advanced", 0))
    ranked = stage5_to_int(metrics.get("ranked_candidates", 0))
    evidence_score = stage5_to_float(target.get("evidence_score"), "")

    gates = [
        {
            "gate_id": "target_evidence",
            "gate_name": "Target evidence package",
            "status": stage6_status_from_readiness(readiness.get("target_evidence")),
            "evidence": f"readiness={readiness.get('target_evidence', '')}; evidence_score={evidence_score}",
            "required_next_step": "Keep target source links and positive-control rationale visible in reports.",
        },
        {
            "gate_id": "candidate_intake",
            "gate_name": "Candidate intake and filtering",
            "status": "pass" if stage5_to_int(metrics.get("filtered_candidates", 0)) > 0 else "fail",
            "evidence": f"raw={metrics.get('raw_candidates', 0)}; filtered={metrics.get('filtered_candidates', 0)}",
            "required_next_step": "Preserve raw and filtered candidate files for auditability.",
        },
        {
            "gate_id": "rdkit_validation",
            "gate_name": "RDKit descriptor validation",
            "status": "pass" if valid_rdkit > 0 else "fail",
            "evidence": f"valid_rdkit={valid_rdkit}; invalid_rdkit={metrics.get('invalid_rdkit', 0)}",
            "required_next_step": "Use real RDKit descriptors as the minimum chemistry-validity layer.",
        },
        {
            "gate_id": "control_panel",
            "gate_name": "Known-control benchmark panel",
            "status": "pass" if controls > 0 else "fail",
            "evidence": f"controls={controls}",
            "required_next_step": "Keep known drugs as controls only; never label them as generated discoveries.",
        },
        {
            "gate_id": "decoy_panel",
            "gate_name": "Decoy/background panel",
            "status": "pass" if decoys > 0 else "warn",
            "evidence": f"decoys={decoys}",
            "required_next_step": "Expand decoys before claiming enrichment behavior.",
        },
        {
            "gate_id": "rank_feedback",
            "gate_name": "Ranking and feedback loop",
            "status": "pass" if ranked > 0 and advanced > 0 else "warn",
            "evidence": f"ranked={ranked}; advanced={advanced}",
            "required_next_step": "Feed advanced candidates into the next generator round or validation queue.",
        },
        {
            "gate_id": "real_docking",
            "gate_name": "Real docking / pose validation",
            "status": stage6_status_from_readiness(readiness.get("docking") or docking_plan.get("status", ""), docking=True),
            "evidence": f"docking_status={readiness.get('docking') or docking_plan.get('status', '')}; backend={docking_plan.get('selected_backend', '')}",
            "required_next_step": "Install/configure Vina or GNINA plus ligand/receptor preparation before making structure-score claims.",
        },
        {
            "gate_id": "claim_boundary",
            "gate_name": "Scientific claim boundary",
            "status": "pass",
            "evidence": "Outputs are labeled computational screening and not efficacy proof.",
            "required_next_step": "Keep this boundary visible in frontend, reports, and presentation material.",
        },
    ]
    return gates


def stage6_gate_summary(gates: Sequence[Dict[str, object]]) -> Dict[str, object]:
    counts = Counter(str(row.get("status", "")) for row in gates)
    if counts.get("fail", 0) > 0:
        overall = "blocked"
    elif counts.get("warn", 0) > 0:
        overall = "computational_demo_ready_real_docking_missing"
    else:
        overall = "ready_for_computational_demo"
    return {
        "overall_status": overall,
        "pass": counts.get("pass", 0),
        "warn": counts.get("warn", 0),
        "fail": counts.get("fail", 0),
        "total": len(gates),
    }


def stage6_descriptors_by_id(project: Path, round_no: int) -> Dict[str, Dict[str, str]]:
    rows = read_csv(project / "stage4" / f"round_{round_no}_real_descriptors.csv")
    return {row.get("id", ""): row for row in rows if row.get("id")}


def stage6_best_benchmark_by_id(project: Path, round_no: int) -> Dict[str, Dict[str, str]]:
    rows = read_csv(project / "stage4" / f"round_{round_no}_benchmark_panel.csv")
    return {row.get("id", ""): row for row in rows if row.get("panel_type") == "candidate" and row.get("id")}


def stage6_hit_triage_rows(project: Path, round_no: int, top: int) -> List[Dict[str, object]]:
    ranked = read_csv(project / "ranked" / f"round_{round_no}_ranked.csv")
    ranked.sort(key=lambda row: stage5_to_int(row.get("rank"), 999999))
    descriptors = stage6_descriptors_by_id(project, round_no)
    benchmark = stage6_best_benchmark_by_id(project, round_no)
    selected = [row for row in ranked if row.get("decision") == "advance"][:top]
    if not selected:
        selected = ranked[:top]

    rows: List[Dict[str, object]] = []
    for row in selected:
        item_id = row.get("id", "")
        desc = descriptors.get(item_id, {})
        bench = benchmark.get(item_id, {})
        score = stage5_to_float(row.get("total_proxy"), 0.0)
        qed = stage5_to_float(desc.get("qed") or row.get("qed_proxy"), "")
        lipinski = stage5_to_int(desc.get("lipinski_violations", 0))
        valid = str(desc.get("valid", ""))
        if valid == "1" and lipinski <= 1:
            tier = "tier_1_structure_validation"
        elif valid == "1":
            tier = "tier_2_property_review"
        else:
            tier = "reject_or_repair"
        next_action = (
            "Run real docking and PoseBusters before promotion; keep as computational hit only."
            if tier == "tier_1_structure_validation"
            else "Review chemistry flags before docking."
        )
        rows.append(
            {
                "rank": row.get("rank", ""),
                "id": item_id,
                "smiles": row.get("smiles", ""),
                "total_score": score,
                "qed": qed,
                "mw": desc.get("mw", ""),
                "lipinski_violations": lipinski,
                "rdkit_valid": valid,
                "nearest_control": bench.get("nearest_control", ""),
                "control_similarity": bench.get("max_similarity_to_controls", ""),
                "validation_tier": tier,
                "priority": "high" if stage5_to_float(score, 0.0) >= 0.64 and tier == "tier_1_structure_validation" else "medium",
                "next_action": next_action,
                "claim_allowed": "computational prioritization only",
            }
        )
    return rows


def stage6_assay_queue_rows(
    triage_rows: Sequence[Dict[str, object]],
    data: Dict[str, object],
    project: Path,
    round_no: int,
) -> List[Dict[str, object]]:
    target = data.get("target", {})
    tables = data.get("tables", {})
    if not isinstance(target, dict):
        target = {}
    if not isinstance(tables, dict):
        tables = {}
    docking_plan = tables.get("docking_plan", {})
    if not isinstance(docking_plan, dict):
        docking_plan = {}
    expected_scores = docking_plan.get("expected_scores_csv") or str(project / "stage4" / f"round_{round_no}_docking_scores_template.csv")
    rows: List[Dict[str, object]] = []
    for item in triage_rows:
        item_id = item.get("id", "")
        rows.append(
            {
                "queue_type": "computational_docking",
                "priority": item.get("priority", ""),
                "id": item_id,
                "smiles": item.get("smiles", ""),
                "target_id": target.get("id", ""),
                "input_or_output": expected_scores,
                "acceptance_criterion": "Real Vina/GNINA score exists and pose_pass is true or documented.",
                "owner_note": "Prepare receptor and ligand files before running docking.",
            }
        )
        rows.append(
            {
                "queue_type": "pose_quality",
                "priority": item.get("priority", ""),
                "id": item_id,
                "smiles": item.get("smiles", ""),
                "target_id": target.get("id", ""),
                "input_or_output": "PoseBusters or equivalent pose sanity report",
                "acceptance_criterion": "No severe geometry or interaction artifacts.",
                "owner_note": "Pose QC is required before discussing structure-based ranking.",
            }
        )
    rows.append(
        {
            "queue_type": "wet_lab_assay_planning",
            "priority": "planning_only",
            "id": "top_hits_panel",
            "smiles": "; ".join(str(item.get("smiles", "")) for item in triage_rows if item.get("smiles")),
            "target_id": target.get("id", ""),
            "input_or_output": "assay_plan_metadata_only",
            "acceptance_criterion": "Only plan assay handoff after real docking/pose QC; no efficacy claim before measured data.",
            "owner_note": "Suggested path: target-specific biochemical assay with known positive controls and decoys.",
        }
    )
    return rows


def stage6_risk_rows(data: Dict[str, object]) -> List[Dict[str, object]]:
    readiness = data.get("readiness", {})
    metrics = data.get("metrics", {})
    if not isinstance(readiness, dict):
        readiness = {}
    if not isinstance(metrics, dict):
        metrics = {}
    risks = [
        {
            "risk_id": "proxy_score_overclaim",
            "severity": "high",
            "trigger": "Scores include proxy or RDKit-only fields.",
            "mitigation": "Display proxy/RDKit labels and require external docking CSV before structure-score claims.",
        },
        {
            "risk_id": "docking_missing",
            "severity": "high" if str(readiness.get("docking", "")).lower() in {"skipped", "not_available", "missing"} else "medium",
            "trigger": f"docking={readiness.get('docking', '')}",
            "mitigation": "Install Vina/GNINA and receptor/ligand preparation tools; rerun score with external docking results.",
        },
        {
            "risk_id": "small_decoy_panel",
            "severity": "medium" if stage5_to_int(metrics.get("decoys", 0)) < 20 else "low",
            "trigger": f"decoys={metrics.get('decoys', 0)}",
            "mitigation": "Expand decoys and benchmark controls before enrichment claims.",
        },
        {
            "risk_id": "wet_lab_gap",
            "severity": "high",
            "trigger": "No measured assay results in project files.",
            "mitigation": "Treat wet-lab validation as a future handoff path, not a current result.",
        },
    ]
    return risks


def render_stage6_report(
    round_no: int,
    assets: Dict[str, object],
    gates: Sequence[Dict[str, object]],
    triage_rows: Sequence[Dict[str, object]],
    queue_rows: Sequence[Dict[str, object]],
) -> str:
    lines = [
        f"# Stage 6 Validation Operations - Round {round_no}",
        "",
        "## Summary",
        "",
        f"- Overall status: `{assets.get('overall_status', '')}`",
        f"- Docking status: `{assets.get('docking_status', '')}`",
        f"- Triage hits: {len(triage_rows)}",
        f"- Validation queue rows: {len(queue_rows)}",
        "- Interpretation: validation operations plan, not efficacy proof.",
        "",
        "## Quality Gates",
        "",
        "| gate | status | evidence | next step |",
        "|---|---|---|---|",
    ]
    for row in gates:
        lines.append(
            f"| {row.get('gate_id', '')} | `{row.get('status', '')}` | "
            f"{row.get('evidence', '')} | {row.get('required_next_step', '')} |"
        )
    lines.extend(
        [
            "",
            "## Hit Triage",
            "",
            "| rank | id | total | tier | next action |",
            "|---:|---|---:|---|---|",
        ]
    )
    for row in triage_rows:
        lines.append(
            f"| {row.get('rank', '')} | {row.get('id', '')} | {row.get('total_score', '')} | "
            f"{row.get('validation_tier', '')} | {row.get('next_action', '')} |"
        )
    lines.extend(
        [
            "",
            "## Queue Types",
            "",
        ]
    )
    queue_counts = Counter(str(row.get("queue_type", "")) for row in queue_rows)
    for key, value in queue_counts.items():
        lines.append(f"- {key}: {value}")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- Stage 6 turns screening outputs into validation operations and quality gates.",
            "- It does not establish potency, activity, safety, or therapeutic efficacy.",
            "- `warn` on real docking means the product demo can continue, but structure-score claims must wait.",
        ]
    )
    return "\n".join(lines) + "\n"


def stage6_validate(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    ensure_project_dirs(project)
    round_no = args.round
    stage5_path = project / "stage5" / "dashboard_data.json"
    data = stage5_read_json(stage5_path)
    if not data or stage5_to_int(data.get("round"), -1) != round_no:
        data = build_stage5_dashboard_data(project, round_no, f"{project.name} round {round_no} closed-loop dashboard")
        write_json(stage5_path, data)
    gates = stage6_quality_gate_rows(data)
    summary = stage6_gate_summary(gates)
    triage_rows = stage6_hit_triage_rows(project, round_no, int(getattr(args, "top", 10) or 10))
    queue_rows = stage6_assay_queue_rows(triage_rows, data, project, round_no)
    risk_rows = stage6_risk_rows(data)

    stage6_dir = project / "stage6"
    assets_path = stage6_dir / f"round_{round_no}_validation_assets.json"
    gates_path = stage6_dir / f"round_{round_no}_quality_gates.csv"
    triage_path = stage6_dir / f"round_{round_no}_hit_triage.csv"
    queue_path = stage6_dir / f"round_{round_no}_assay_queue.csv"
    risk_path = stage6_dir / f"round_{round_no}_risk_register.csv"
    runbook_path = stage6_dir / f"round_{round_no}_validation_runbook.md"
    report_path = project / "reports" / f"stage6_round_{round_no}_validation_report.md"

    write_csv(gates_path, gates, ["gate_id", "gate_name", "status", "evidence", "required_next_step"])
    write_csv(
        triage_path,
        triage_rows,
        [
            "rank",
            "id",
            "smiles",
            "total_score",
            "qed",
            "mw",
            "lipinski_violations",
            "rdkit_valid",
            "nearest_control",
            "control_similarity",
            "validation_tier",
            "priority",
            "next_action",
            "claim_allowed",
        ],
    )
    write_csv(
        queue_path,
        queue_rows,
        ["queue_type", "priority", "id", "smiles", "target_id", "input_or_output", "acceptance_criterion", "owner_note"],
    )
    write_csv(risk_path, risk_rows, ["risk_id", "severity", "trigger", "mitigation"])

    assets = {
        "schema_version": "0.1",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "stage": 6,
        "round": round_no,
        "overall_status": summary["overall_status"],
        "docking_status": nested_get(data, ["readiness", "docking"]) or "",
        "quality_gate_summary": summary,
        "triage_count": len(triage_rows),
        "queue_count": len(queue_rows),
        "risk_count": len(risk_rows),
        "files": {
            "stage5_dashboard_data": str(stage5_path),
            "quality_gates": str(gates_path),
            "hit_triage": str(triage_path),
            "assay_queue": str(queue_path),
            "risk_register": str(risk_path),
            "validation_runbook": str(runbook_path),
            "report": str(report_path),
        },
        "next_actions": [
            "Run real docking and pose QC for tier_1_structure_validation hits.",
            "Feed external docking CSV back through score --external-scores.",
            "Expand decoys and controls before claiming enrichment.",
            "Keep wet-lab assay planning separate from current computational evidence.",
        ],
        "boundary": [
            "Stage 6 is validation operations planning.",
            "No therapeutic, clinical, safety, or measured-potency claim is created here.",
        ],
    }
    write_json(assets_path, assets)
    write_text(runbook_path, render_stage6_report(round_no, assets, gates, triage_rows, queue_rows))
    write_text(report_path, render_stage6_report(round_no, assets, gates, triage_rows, queue_rows))

    print(f"Wrote Stage 6 validation assets: {assets_path}")
    print(f"Wrote Stage 6 quality gates: {gates_path}")
    print(f"Wrote Stage 6 hit triage: {triage_path}")
    print(f"Wrote Stage 6 assay queue: {queue_path}")
    print(f"Wrote Stage 6 risk register: {risk_path}")
    print(f"Wrote Stage 6 report: {report_path}")
    return assets_path


def stage7_deliverable_record(path: Path, stage7_dir: Path, description: str) -> Dict[str, object]:
    return {
        "description": description,
        "path": str(path),
        "relative_path": stage5_relative_path(path, stage7_dir),
        "exists": path.exists(),
    }


def render_stage7_executive_summary(
    title: str,
    round_no: int,
    stage5_data: Dict[str, object],
    stage6_assets: Dict[str, object],
) -> str:
    target = stage5_data.get("target", {})
    metrics = stage5_data.get("metrics", {})
    readiness = stage5_data.get("readiness", {})
    if not isinstance(target, dict):
        target = {}
    if not isinstance(metrics, dict):
        metrics = {}
    if not isinstance(readiness, dict):
        readiness = {}
    lines = [
        f"# {title}",
        "",
        "## One-Line Positioning",
        "",
        "This is a productized AI molecule design closed loop for computational screening, target-evidence traceability, real RDKit validation, and validation-operations planning.",
        "",
        "## Current Demo Scope",
        "",
        f"- Round: `{round_no}`",
        f"- Target: `{target.get('id', '')}` / {target.get('name', '')}",
        f"- Evidence readiness: `{readiness.get('target_evidence', '')}`",
        f"- Candidate rows: {metrics.get('candidate_rows', 0)}",
        f"- RDKit-valid candidates: {metrics.get('valid_rdkit', 0)}",
        f"- Advanced hits: {metrics.get('advanced', 0)}",
        f"- Benchmark panel: {metrics.get('benchmark_panel', 0)} molecules",
        f"- Stage 6 status: `{stage6_assets.get('overall_status', '')}`",
        f"- Docking status: `{readiness.get('docking', '')}`",
        "",
        "## Product Claim",
        "",
        "The system demonstrates an auditable computational screening workflow. It does not claim measured potency, antiviral efficacy, clinical benefit, safety, dosing, or synthesis feasibility.",
        "",
        "## Next Product Milestones",
        "",
        "1. Install and configure real docking and pose-quality backends.",
        "2. Import docking/pose CSV into the score and feedback loop.",
        "3. Build the Stage 8 web product around project setup, evidence, candidate funnel, validation operations, and delivery export.",
    ]
    return "\n".join(lines) + "\n"


def render_stage7_reproducibility(project: Path, round_no: int) -> str:
    cli = Path(__file__).resolve()
    lines = [
        f"# Reproducibility Runbook - Round {round_no}",
        "",
        "Run these commands from the repository workspace.",
        "",
        "```bash",
        f"python3 {cli} target-select {project} --disease 甲流 --top 5",
        f"python3 {cli} brief-from-target {project} --disease 甲流 --target influenza_a_h1n1_na --force",
        f"python3 {cli} evidence-stage2 {project} --disease influenza --top 5",
        f"python3 {cli} stage3-screen {project} --round {round_no} --top 8",
        f"python3 {cli} stage4-real {project} --round {round_no} --target influenza_a_h1n1_na --top 8 --decoys 8 --rescore",
        f"python3 {cli} stage5-dashboard {project} --round {round_no}",
        f"python3 {cli} stage6-validate {project} --round {round_no} --top 8",
        f"python3 {cli} stage7-package {project} --round {round_no}",
        "```",
        "",
        "If real docking is available, insert this after Stage 4 writes docking inputs:",
        "",
        "```bash",
        f"python3 {cli} score {project} --round {round_no} --external-scores <real_docking_scores.csv>",
        f"python3 {cli} rank {project} --round {round_no} --top 8",
        f"python3 {cli} feedback {project} --round {round_no} --top 8",
        "```",
        "",
        "Boundary: repeated runs reproduce project files and computational prioritization, not measured efficacy.",
    ]
    return "\n".join(lines) + "\n"


def stage7_demo_checklist_rows() -> List[Dict[str, object]]:
    return [
        {
            "section": "demo_story",
            "item": "Open with the target-evidence problem and why undefined targets break AI molecule design demos.",
            "status": "ready",
            "owner_note": "Use Stage 2 target evidence and Stage 5 dashboard.",
        },
        {
            "section": "demo_story",
            "item": "Show candidate funnel from generation/intake to RDKit validation and feedback.",
            "status": "ready",
            "owner_note": "Use Stage 5 metrics and ranked table.",
        },
        {
            "section": "technical_depth",
            "item": "Explain that RDKit descriptors, fingerprints, scaffolds, SDF, controls, and decoys are real local computations.",
            "status": "ready",
            "owner_note": "Avoid saying docking was completed unless external scores exist.",
        },
        {
            "section": "risk_boundary",
            "item": "State that this is computational screening and validation planning, not efficacy proof.",
            "status": "must_show",
            "owner_note": "Keep in UI and presentation.",
        },
        {
            "section": "next_milestone",
            "item": "Position real docking, pose QC, and assay handoff as the next product milestones.",
            "status": "ready",
            "owner_note": "Use Stage 6 queue and risk register.",
        },
    ]


def render_stage8_frontend_spec(project: Path, round_no: int, title: str) -> str:
    return f"""# Stage 8 Frontend Product Specification

## Product Goal

Build a full web product that turns the current CLI pipeline into an operable AI molecule screening workspace. The frontend should let a user define a target, inspect evidence, intake or generate candidates, run validation stages, review quality gates, and export a delivery package without reading raw CSV/JSON files.

## Primary Users

- Project owner: needs a credible demo and progress dashboard.
- Computational chemistry reviewer: needs evidence traceability, descriptors, controls, decoys, docking status, and risk boundaries.
- Competition or investor audience: needs a clear story without overclaiming drug efficacy.

## Information Architecture

### 1. Project Command Center

First screen after opening a project. It should show project name, target, round, current stage, key readiness badges, candidate funnel, and the most important next action. It should not look like a marketing landing page. Use a dense operational dashboard with clear sections and compact cards.

Core widgets:

- Stage progress rail: Stage 1 target, Stage 2 evidence, Stage 3 intake, Stage 4 real library, Stage 5 dashboard, Stage 6 validation ops, Stage 7 delivery.
- Readiness cards: target evidence, candidate intake, RDKit, control panel, decoy panel, docking, feedback.
- Next-action panel: generated from Stage 6 quality gates and validation queue.
- Open files panel: links to dashboard, reports, manifest, and reproducibility runbook.

Primary data:

- `stage5/dashboard_data.json`
- `stage6/round_{round_no}_validation_assets.json`
- `stage7/round_{round_no}_delivery_manifest.json`

### 2. Target Evidence Workspace

Purpose: make the target decision defensible before molecule generation.

Layout:

- Left panel: disease/query, selected target, PDB, reference ligand, key residues, known controls.
- Main table: target candidates, evidence score, readiness, PDB count, PubMed metadata count, assay path.
- Detail drawer: source URLs, structure notes, positive controls, caveats.

Required interactions:

- Select target from catalog.
- Compare influenza A NA, PA endonuclease, M2, influenza B targets.
- Regenerate target brief and generator prompt.
- Display warning when target lacks positive controls or PDB structure.

### 3. Candidate Funnel

Purpose: show how generated or imported candidates move through filtering, scoring, ranking, and feedback.

Layout:

- Funnel counters: raw, filtered, valid RDKit, ranked, advanced.
- Candidate table with sortable columns: ID, SMILES, total score, QED, MW, Lipinski violations, decision, score source.
- Molecule image strip from Stage 4 PNG files.
- Feedback seeds panel for the next round.

Required interactions:

- Import CSV/JSON/URL candidate source.
- Trigger OpenAI candidate generation with key source selection, without storing API keys.
- Run filtering/scoring/ranking commands.
- Export selected candidates.

### 4. Real Library Validation

Purpose: distinguish real local chemistry calculations from proxy scoring.

Layout:

- RDKit descriptor table.
- Fingerprint similarity to controls.
- Candidate/control/decoy benchmark panel.
- SDF and docking input readiness status.
- Docking backend detector: Vina, GNINA, OpenBabel, Meeko, PoseBusters.

Required interactions:

- Run Stage 4.
- Download or open SDF/CSV artifacts.
- Attach external docking scores CSV.
- Re-score and re-rank with real docking results.

### 5. Validation Operations

Purpose: convert computational output into a responsible validation plan.

Layout:

- Quality gate matrix with pass/warn/fail.
- Hit triage table.
- Computational docking queue.
- Pose-quality queue.
- Wet-lab assay planning row shown as planning-only.
- Risk register.

Required interactions:

- Mark queue items as planned/running/completed manually.
- Attach external score files.
- Regenerate Stage 6 report.
- Block efficacy language while docking or assay data is missing.

### 6. Delivery Room

Purpose: create a polished output package for a competition, team review, or investor demo.

Layout:

- Executive summary preview.
- Delivery manifest with file presence checks.
- Reproducibility command runbook.
- Demo checklist.
- Claims boundary panel.
- Export buttons for static dashboard, report bundle, CSV package, and future PPT.

Required interactions:

- Generate Stage 7 package.
- Open generated files.
- Copy reproducibility commands.
- Surface missing deliverables before export.

## Navigation

Use a left sidebar with stages and a top project bar:

- Project selector
- Target name
- Round selector
- Run status
- Global command button

Avoid hero sections. The product is an operational tool, so the first screen must be the dashboard, not a landing page.

## Visual Style

- Quiet, technical, compact.
- Use neutral background, strong table readability, restrained accent colors.
- Use green/warn/red readiness badges.
- Use molecule images as real visual assets.
- No decorative gradients, floating marketing cards, or oversized hero typography.

## Data Contract

Frontend should treat CLI outputs as the source of truth:

- Stage 1/2: `targets/target_selection.csv`, `evidence/stage2_target_sources.csv`, `evidence/stage2_closed_loop_assets.json`
- Stage 3: `stage3/round_N_stage3_assets.json`, `filtered/round_N_filtered.csv`, `candidates/round_N_candidates.csv`
- Stage 4: `stage4/round_N_stage4_assets.json`, `stage4/round_N_real_descriptors.csv`, `stage4/round_N_benchmark_panel.csv`, `stage4/round_N_validation_metrics.json`
- Stage 5: `stage5/dashboard_data.json`
- Stage 6: `stage6/round_N_validation_assets.json`, `stage6/round_N_quality_gates.csv`, `stage6/round_N_hit_triage.csv`, `stage6/round_N_assay_queue.csv`
- Stage 7: `stage7/round_N_delivery_manifest.json`

## Backend/API Shape

Initial product can wrap the CLI through a local API:

- `GET /projects`
- `GET /projects/:id/dashboard?round=N`
- `POST /projects/:id/stage3-screen`
- `POST /projects/:id/stage4-real`
- `POST /projects/:id/stage5-dashboard`
- `POST /projects/:id/stage6-validate`
- `POST /projects/:id/stage7-package`
- `GET /projects/:id/files/:path`

Long-running jobs should stream logs and write final artifacts to the existing project directories.

## Safety And Claims Boundary

The UI must prevent accidental overclaiming:

- Replace “validated drug” with “computational hit” unless measured assay data exists.
- Mark docking as missing/skipped/not available unless real backend output is imported.
- Display known drugs as controls only.
- Show “not efficacy proof” in dashboard, validation, and delivery views.

## MVP Build Order

1. Project Command Center from Stage 5 JSON.
2. Validation Operations from Stage 6 CSV/JSON.
3. Delivery Room from Stage 7 manifest and reports.
4. Candidate Funnel table and molecule gallery.
5. Target Evidence Workspace.
6. CLI job runner and logs.

## Acceptance Criteria

- User can open one project and understand target, round, stage status, and next action in under 30 seconds.
- User can trace every displayed metric back to a local file.
- User cannot mistake proxy/RDKit-only output for real docking or wet-lab validation.
- Stage 6 quality gates directly drive visible next actions.
- Stage 7 export can be generated and opened from the frontend.

Generated for `{project}` round `{round_no}` as part of `{title}`.
"""


def stage7_manifest(
    project: Path,
    round_no: int,
    stage7_dir: Path,
    stage5_data: Dict[str, object],
    stage6_assets: Dict[str, object],
    stage8_spec_path: Path,
) -> Dict[str, object]:
    deliverables = {
        "stage5_dashboard": stage7_deliverable_record(project / "stage5" / "index.html", stage7_dir, "Static dashboard for product demo."),
        "stage5_data": stage7_deliverable_record(project / "stage5" / "dashboard_data.json", stage7_dir, "Aggregated frontend data package."),
        "stage6_validation_assets": stage7_deliverable_record(project / "stage6" / f"round_{round_no}_validation_assets.json", stage7_dir, "Validation operations assets."),
        "stage6_quality_gates": stage7_deliverable_record(project / "stage6" / f"round_{round_no}_quality_gates.csv", stage7_dir, "Quality gate matrix."),
        "stage6_hit_triage": stage7_deliverable_record(project / "stage6" / f"round_{round_no}_hit_triage.csv", stage7_dir, "Top-hit validation triage."),
        "stage6_assay_queue": stage7_deliverable_record(project / "stage6" / f"round_{round_no}_assay_queue.csv", stage7_dir, "Computational and planning queue."),
        "stage6_risk_register": stage7_deliverable_record(project / "stage6" / f"round_{round_no}_risk_register.csv", stage7_dir, "Validation and communication risk register."),
        "stage6_report": stage7_deliverable_record(project / "reports" / f"stage6_round_{round_no}_validation_report.md", stage7_dir, "Stage 6 validation operations report."),
        "stage4_descriptors": stage7_deliverable_record(project / "stage4" / f"round_{round_no}_real_descriptors.csv", stage7_dir, "RDKit descriptors."),
        "stage4_benchmark_panel": stage7_deliverable_record(project / "stage4" / f"round_{round_no}_benchmark_panel.csv", stage7_dir, "Candidate/control/decoy panel."),
        "stage4_report": stage7_deliverable_record(project / "reports" / f"stage4_round_{round_no}_report.md", stage7_dir, "Stage 4 real-library validation report."),
        "stage3_report": stage7_deliverable_record(project / "reports" / f"stage3_round_{round_no}_report.md", stage7_dir, "Stage 3 candidate intake and screening report."),
        "ranked_candidates": stage7_deliverable_record(project / "ranked" / f"round_{round_no}_ranked.csv", stage7_dir, "Ranked candidates."),
        "executive_summary": stage7_deliverable_record(stage7_dir / f"round_{round_no}_executive_summary.md", stage7_dir, "Product executive summary."),
        "reproducibility_runbook": stage7_deliverable_record(stage7_dir / f"round_{round_no}_reproducibility.md", stage7_dir, "Commands for reproducing the workflow artifacts."),
        "demo_checklist": stage7_deliverable_record(stage7_dir / f"round_{round_no}_investor_demo_checklist.csv", stage7_dir, "Demo and communication checklist."),
        "stage8_frontend_spec": stage7_deliverable_record(stage8_spec_path, stage7_dir, "Frontend product specification for Stage 8."),
    }
    metrics = stage5_data.get("metrics", {})
    readiness = stage5_data.get("readiness", {})
    if not isinstance(metrics, dict):
        metrics = {}
    if not isinstance(readiness, dict):
        readiness = {}
    return {
        "schema_version": "0.1",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "stage": 7,
        "round": round_no,
        "project": str(project),
        "delivery_status": "ready_for_demo_package" if all(item["exists"] for item in deliverables.values()) else "missing_some_deliverables",
        "stage6_status": stage6_assets.get("overall_status", ""),
        "metrics_snapshot": metrics,
        "readiness_snapshot": readiness,
        "deliverables": deliverables,
        "claims_boundary": [
            "Computational screening and validation planning only.",
            "No measured potency, efficacy, clinical benefit, dosing, or safety claim.",
            "Known drugs are controls and benchmarks.",
            "Real docking claims require imported or generated external docking scores.",
        ],
    }


def render_stage7_delivery_report(
    round_no: int,
    manifest: Dict[str, object],
    executive_summary: Path,
    reproducibility: Path,
    stage8_spec: Path,
) -> str:
    deliverables = manifest.get("deliverables", {})
    if not isinstance(deliverables, dict):
        deliverables = {}
    lines = [
        f"# Stage 7 Delivery Package - Round {round_no}",
        "",
        f"- Delivery status: `{manifest.get('delivery_status', '')}`",
        f"- Stage 6 status: `{manifest.get('stage6_status', '')}`",
        f"- Executive summary: `{executive_summary}`",
        f"- Reproducibility runbook: `{reproducibility}`",
        f"- Stage 8 frontend spec: `{stage8_spec}`",
        "",
        "## Deliverables",
        "",
        "| key | exists | path | description |",
        "|---|---|---|---|",
    ]
    for key, item in deliverables.items():
        if not isinstance(item, dict):
            continue
        lines.append(f"| {key} | {item.get('exists', '')} | `{item.get('path', '')}` | {item.get('description', '')} |")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- Stage 7 packages product artifacts and reproducibility material.",
            "- It does not add new scientific validation.",
            "- The Stage 8 specification describes the next frontend product, not an implemented web app.",
        ]
    )
    return "\n".join(lines) + "\n"


def stage7_package(args: argparse.Namespace) -> Path:
    project = project_path(args.project)
    ensure_project_dirs(project)
    round_no = args.round
    title = getattr(args, "title", "") or f"{project.name} round {round_no} delivery package"
    stage5_path = project / "stage5" / "dashboard_data.json"
    stage6_path = project / "stage6" / f"round_{round_no}_validation_assets.json"
    stage5_data = stage5_read_json(stage5_path)
    if not stage5_data or stage5_to_int(stage5_data.get("round"), -1) != round_no:
        data = build_stage5_dashboard_data(project, round_no, f"{project.name} round {round_no} closed-loop dashboard")
        write_json(stage5_path, data)
        stage5_data = data
    stage6_assets = stage5_read_json(stage6_path)
    if not stage6_assets:
        stage6_validate(argparse.Namespace(project=str(project), round=round_no, top=8))
        stage6_assets = stage5_read_json(stage6_path)

    stage7_dir = project / "stage7"
    summary_path = stage7_dir / f"round_{round_no}_executive_summary.md"
    repro_path = stage7_dir / f"round_{round_no}_reproducibility.md"
    checklist_path = stage7_dir / f"round_{round_no}_investor_demo_checklist.csv"
    stage8_spec_path = stage7_dir / "stage8_frontend_product_spec.md"
    manifest_path = stage7_dir / f"round_{round_no}_delivery_manifest.json"
    report_path = project / "reports" / f"stage7_round_{round_no}_delivery_report.md"

    write_text(summary_path, render_stage7_executive_summary(title, round_no, stage5_data, stage6_assets))
    write_text(repro_path, render_stage7_reproducibility(project, round_no))
    write_csv(checklist_path, stage7_demo_checklist_rows(), ["section", "item", "status", "owner_note"])
    write_text(stage8_spec_path, render_stage8_frontend_spec(project, round_no, title))
    manifest = stage7_manifest(project, round_no, stage7_dir, stage5_data, stage6_assets, stage8_spec_path)
    write_json(manifest_path, manifest)
    write_text(report_path, render_stage7_delivery_report(round_no, manifest, summary_path, repro_path, stage8_spec_path))

    print(f"Wrote Stage 7 delivery manifest: {manifest_path}")
    print(f"Wrote Stage 7 executive summary: {summary_path}")
    print(f"Wrote Stage 7 reproducibility runbook: {repro_path}")
    print(f"Wrote Stage 7 demo checklist: {checklist_path}")
    print(f"Wrote Stage 8 frontend spec: {stage8_spec_path}")
    print(f"Wrote Stage 7 report: {report_path}")
    return manifest_path


def run_demo(args: argparse.Namespace) -> None:
    project = project_path(args.project)
    if not (project / "config.json").exists():
        init_args = argparse.Namespace(project=str(project), force=False)
        init_project(init_args)
    for round_no in range(1, args.rounds + 1):
        print(f"\n=== Round {round_no} ===")
        generate_candidates(
            argparse.Namespace(project=str(project), round=round_no, n=args.n, source_csv=None)
        )
        score_candidates(
            argparse.Namespace(project=str(project), round=round_no, external_scores=args.external_scores)
        )
        rank_candidates(argparse.Namespace(project=str(project), round=round_no, top=args.top))
        feedback(argparse.Namespace(project=str(project), round=round_no, top=args.top))
    print(f"\nDemo completed: {project}")


def doctor(args: argparse.Namespace) -> None:
    project = project_path(args.project) if args.project else None
    config = load_config(project) if project and (project / "config.json").exists() else default_config()
    tool_paths = config.get("tool_paths", {})
    if not isinstance(tool_paths, dict):
        tool_paths = {}

    print("Local repository/tool path check:")
    for name, raw_path in tool_paths.items():
        path = Path(str(raw_path)).expanduser()
        status = "present" if path.exists() else "missing"
        print(f"  {name:14s} {status:8s} {path}")

    print("\nExecutable check:")
    for executable in [
        "reinvent",
        "drugex",
        "vina",
        "gnina",
        "obabel",
        "mk_prepare_ligand.py",
        "mk_prepare_receptor.py",
        "bust",
        "openfe",
    ]:
        found = shutil.which(executable)
        print(f"  {executable:14s} {'found' if found else 'not found':8s} {found or ''}")

    capabilities = stage4_capabilities()
    modules = capabilities.get("modules", {})
    print("\nPython real-library check:")
    if isinstance(modules, dict):
        for name in ["rdkit", "numpy", "pandas", "openbabel", "meeko", "vina", "posebusters"]:
            status = modules.get(name, {})
            if not isinstance(status, dict):
                continue
            version = status.get("version", "")
            origin = status.get("origin", "")
            print(f"  {name:14s} {status.get('status', ''):12s} {version!s:10s} {origin}")

    docking = capabilities.get("docking_backend", {})
    if not isinstance(docking, dict):
        docking = {}
    ligand_prep = capabilities.get("ligand_preparation", {})
    if not isinstance(ligand_prep, dict):
        ligand_prep = {}

    print("\nInterpretation:")
    print("  Missing executables are expected unless those packages were installed into the current shell.")
    print("  The CLI can still run proxy closed-loop demos and import external scoring CSV files.")
    print("  Stage 4 real-library path requires RDKit; it exports descriptors, fingerprints, scaffolds, and ligand SDF.")
    print(f"  Docking backend: {docking.get('status', '')} - {docking.get('message', '')}")
    print(f"  Ligand preparation: {ligand_prep.get('message', '')}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ai_mol_loop",
        description="Prototype CLI for AI molecule generation, scoring, ranking, and feedback.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    explain_parser = sub.add_parser("explain", help="Explain the closed-loop workflow.")
    explain_parser.set_defaults(func=lambda args: explain())

    init_parser = sub.add_parser("init", help="Initialize a closed-loop project directory.")
    init_parser.add_argument("project", help="Project directory to create or update.")
    init_parser.add_argument("--force", action="store_true", help="Overwrite config.json if it exists.")
    init_parser.set_defaults(func=init_project)

    target_select_parser = sub.add_parser(
        "target-select",
        help="Rank disease-relevant targets from the local target catalog.",
    )
    target_select_parser.add_argument(
        "project",
        nargs="?",
        help="Optional project directory. If provided, writes target ranking and report files.",
    )
    target_select_parser.add_argument(
        "--disease",
        default="influenza_a",
        help="Disease/query filter, e.g. influenza_a, h1n1, 甲流, influenza_b, 乙流, influenza.",
    )
    target_select_parser.add_argument("--target", help="Optional target id/name substring filter.")
    target_select_parser.add_argument("--top", type=int, default=5)
    target_select_parser.add_argument("--catalog", help="Optional target catalog JSON path.")
    target_select_parser.set_defaults(func=target_select)

    brief_from_target_parser = sub.add_parser(
        "brief-from-target",
        help="Create target brief and generator prompt from a catalog target.",
    )
    brief_from_target_parser.add_argument("project")
    brief_from_target_parser.add_argument(
        "--disease",
        default="influenza_a",
        help="Disease/query filter used when --target is omitted.",
    )
    brief_from_target_parser.add_argument(
        "--target",
        help="Target id from the catalog. If omitted, uses the highest-ranked target for --disease.",
    )
    brief_from_target_parser.add_argument("--catalog", help="Optional target catalog JSON path.")
    brief_from_target_parser.add_argument("--free-text", default="", help="Extra user requirement added to the prompt.")
    brief_from_target_parser.add_argument("--max-heavy-atoms", type=int, default=55)
    brief_from_target_parser.add_argument("--max-molecular-weight", type=float, default=550.0)
    brief_from_target_parser.add_argument("--force", action="store_true", help="Overwrite existing brief and prompt.")
    brief_from_target_parser.set_defaults(func=brief_from_target)

    evidence_parser = sub.add_parser(
        "evidence-refresh",
        help="Build or refresh target evidence packages from source registries, RCSB, and PubMed metadata.",
    )
    evidence_parser.add_argument(
        "--disease",
        default="influenza",
        help="Disease/query filter, e.g. influenza, 甲流, 乙流, h1n1.",
    )
    evidence_parser.add_argument("--target", help="Specific target id to refresh.")
    evidence_parser.add_argument("--catalog", help="Optional target catalog JSON path.")
    evidence_parser.add_argument("--sources", help="Optional evidence source registry JSON path.")
    evidence_parser.add_argument("--output", help="Output evidence directory.")
    evidence_parser.add_argument("--retmax", type=int, default=5, help="PubMed articles per query.")
    evidence_parser.add_argument("--timeout", type=int, default=20, help="Network timeout in seconds.")
    evidence_parser.add_argument("--offline", action="store_true", help="Build evidence packages without network calls.")
    evidence_parser.set_defaults(func=evidence_refresh)

    stage2_parser = sub.add_parser(
        "evidence-stage2",
        help="Build project-level stage 2 target-source matrix and closed-loop evidence assets.",
    )
    stage2_parser.add_argument("project", help="Project directory to receive stage 2 evidence outputs.")
    stage2_parser.add_argument(
        "--disease",
        default="influenza",
        help="Disease/query filter, e.g. influenza, 甲流, 乙流, h1n1.",
    )
    stage2_parser.add_argument("--target", help="Specific target id to package.")
    stage2_parser.add_argument("--top", type=int, default=5, help="Number of ranked targets to include when --target is omitted.")
    stage2_parser.add_argument("--catalog", help="Optional target catalog JSON path.")
    stage2_parser.add_argument("--sources", help="Optional evidence source registry JSON path.")
    stage2_parser.add_argument("--evidence-dir", help="Evidence package directory. Defaults to targets/influenza/evidence.")
    stage2_parser.add_argument("--refresh", action="store_true", help="Refresh evidence packages before writing project-level outputs.")
    stage2_parser.add_argument("--retmax", type=int, default=5, help="PubMed articles per query when --refresh is used.")
    stage2_parser.add_argument("--timeout", type=int, default=20, help="Network timeout in seconds when --refresh is used.")
    stage2_parser.add_argument("--offline", action="store_true", help="Use cached/offline evidence when --refresh is used.")
    stage2_parser.set_defaults(func=evidence_stage2)

    brief_parser = sub.add_parser(
        "brief",
        help="Create a structured target brief and generator prompt from target requirements.",
    )
    brief_parser.add_argument("project")
    brief_parser.add_argument("--target-name", required=True, help="Short target name, e.g. SARS-CoV-2 Mpro.")
    brief_parser.add_argument("--disease", default="", help="Disease or indication context.")
    brief_parser.add_argument("--rationale", default="", help="Why this target is relevant.")
    brief_parser.add_argument("--free-text", default="", help="User-facing natural-language requirement.")
    brief_parser.add_argument("--from-file", help="Read the user-facing natural-language requirement from a text file.")
    brief_parser.add_argument("--protein", default="", help="Protein name.")
    brief_parser.add_argument("--gene", default="", help="Gene symbol.")
    brief_parser.add_argument("--pdb-id", default="", help="PDB identifier if available.")
    brief_parser.add_argument("--protein-pdb", default="", help="Local protein PDB path if available.")
    brief_parser.add_argument("--pocket", default="", help="Binding-pocket description.")
    brief_parser.add_argument("--reference-ligand", default="", help="Reference ligand name or SDF path.")
    brief_parser.add_argument("--key-residues", default="", help="Comma/semicolon separated key residues.")
    brief_parser.add_argument("--center", default="", help="Docking box center as 'x,y,z'.")
    brief_parser.add_argument("--size", default="", help="Docking box size as 'x,y,z'.")
    brief_parser.add_argument(
        "--pocket-source",
        default="manual_or_literature",
        help="Pocket source: co_crystal, p2rank, literature, manual, etc.",
    )
    brief_parser.add_argument(
        "--style",
        default="de_novo_and_analog_design",
        help="Design style: de_novo, scaffold_hopping, fragment_growing, analog_design, etc.",
    )
    brief_parser.add_argument("--primary-goal", default="improve predicted binding while retaining drug-like properties.")
    brief_parser.add_argument("--activity", default="in silico binding enrichment, not measured potency.")
    brief_parser.add_argument("--must-have", default="", help="Semicolon separated desired molecular features.")
    brief_parser.add_argument("--avoid", default="", help="Semicolon separated features to avoid.")
    brief_parser.add_argument("--desired-properties", default="", help="Semicolon separated property targets.")
    brief_parser.add_argument("--selectivity", default="", help="Selectivity or off-target notes.")
    brief_parser.add_argument("--max-heavy-atoms", type=int, default=55)
    brief_parser.add_argument("--max-molecular-weight", type=float, default=500.0)
    brief_parser.add_argument("--preserve-scaffold", action="store_true")
    brief_parser.add_argument("--no-diversity", action="store_false", dest="diversity")
    brief_parser.set_defaults(diversity=True)
    brief_parser.add_argument("--force", action="store_true", help="Overwrite existing brief and prompt.")
    brief_parser.set_defaults(func=brief_command)

    examples_parser = sub.add_parser(
        "prompt-examples",
        help="Print target-requirement prompt examples for molecule generation.",
    )
    examples_parser.add_argument("project", nargs="?", help="Optional project directory.")
    examples_parser.add_argument("--write", action="store_true", help="Also write examples into project/prompts/.")
    examples_parser.set_defaults(func=prompt_examples)

    generate_parser = sub.add_parser("generate", help="Generate or import candidate SMILES.")
    generate_parser.add_argument("project")
    generate_parser.add_argument("--round", type=int, default=1)
    generate_parser.add_argument("--n", type=int, default=50, help="Number of proxy candidates to generate.")
    generate_parser.add_argument("--source-csv", help="Import candidates from a CSV with a smiles column.")
    generate_parser.set_defaults(func=generate_candidates)

    score_parser = sub.add_parser("score", help="Score candidates with proxy metrics and optional external results.")
    score_parser.add_argument("project")
    score_parser.add_argument("--round", type=int, default=1)
    score_parser.add_argument(
        "--external-scores",
        help="CSV with id or smiles plus docking_score/affinity/gnina_score/vina_score and optional pose_pass.",
    )
    score_parser.add_argument(
        "--real-descriptors",
        help="Optional Stage 4 RDKit descriptor CSV. Defaults to stage4/round_N_real_descriptors.csv when present.",
    )
    score_parser.set_defaults(func=score_candidates)

    rank_parser = sub.add_parser("rank", help="Rank scored candidates and mark molecules for advancement.")
    rank_parser.add_argument("project")
    rank_parser.add_argument("--round", type=int, default=1)
    rank_parser.add_argument("--top", type=int, default=10)
    rank_parser.set_defaults(func=rank_candidates)

    feedback_parser = sub.add_parser("feedback", help="Write next-round seed and feedback package.")
    feedback_parser.add_argument("project")
    feedback_parser.add_argument("--round", type=int, default=1)
    feedback_parser.add_argument("--top", type=int, default=10)
    feedback_parser.set_defaults(func=feedback)

    stage3_parser = sub.add_parser(
        "stage3-screen",
        help="Productized stage 3 candidate intake, optional OpenAI generation, filtering, scoring, ranking, and feedback.",
    )
    stage3_parser.add_argument("project")
    stage3_parser.add_argument("--round", type=int, default=3)
    stage3_parser.add_argument("--n", type=int, default=30, help="Number of OpenAI/proxy candidates to request or generate.")
    stage3_parser.add_argument("--top", type=int, default=10, help="Top candidates to advance after ranking.")
    stage3_parser.add_argument("--source-csv", help="Candidate CSV with a smiles column.")
    stage3_parser.add_argument("--source-json", help="Candidate JSON file with candidates/molecules rows.")
    stage3_parser.add_argument("--source-url", action="append", default=[], help="HTTP(S) URL containing CSV/JSON/Markdown/text candidate molecules.")
    stage3_parser.add_argument("--context-url", action="append", default=[], help="HTTP(S) URL used as OpenAI context, not parsed as candidate rows.")
    stage3_parser.add_argument("--use-openai", action="store_true", help="Use OpenAI Responses API to generate additional candidate molecules.")
    stage3_parser.add_argument("--openai-model", default="gpt-5.2", help="OpenAI model for candidate generation.")
    stage3_parser.add_argument("--api-key", help="Temporary OpenAI API key. Prefer OPENAI_API_KEY or --api-key-file; never stored in outputs.")
    stage3_parser.add_argument("--api-key-file", help="File containing an OpenAI API key. The key content is not copied into outputs.")
    stage3_parser.add_argument("--api-key-env", default="OPENAI_API_KEY", help="Environment variable name containing the OpenAI API key.")
    stage3_parser.add_argument("--prompt", default="", help="Extra natural-language requirement for OpenAI candidate generation.")
    stage3_parser.add_argument("--timeout", type=int, default=60, help="URL/OpenAI request timeout in seconds.")
    stage3_parser.add_argument("--max-url-bytes", type=int, default=800000, help="Maximum bytes to read from each source/context URL.")
    stage3_parser.add_argument("--max-heavy-atoms", type=int, default=55)
    stage3_parser.add_argument("--max-molecular-weight", type=float, default=550.0)
    stage3_parser.add_argument("--max-lipinski-violations", type=int, default=1)
    stage3_parser.add_argument("--allow-risk", action="store_true", help="Do not reject simple structural risk flags during early filtering.")
    stage3_parser.add_argument("--external-scores", help="Optional external docking/pose scores CSV passed through to score.")
    stage3_parser.add_argument("--no-score", action="store_true", help="Only intake and filter candidates; skip scoring/ranking/feedback.")
    stage3_parser.set_defaults(func=stage3_screen)

    stage4_parser = sub.add_parser(
        "stage4-real",
        help="Run Stage 4 with real local chemistry libraries: RDKit descriptors, fingerprints, SDF export, and control similarity.",
    )
    stage4_parser.add_argument("project")
    stage4_parser.add_argument("--round", type=int, default=4, help="Candidate round to validate.")
    stage4_parser.add_argument(
        "--target",
        default="",
        help="Target id used to select known controls, e.g. influenza_a_h1n1_na. Defaults to config target_catalog_id.",
    )
    stage4_parser.add_argument("--input-csv", help="Optional candidate CSV. Defaults to project/candidates/round_N_candidates.csv.")
    stage4_parser.add_argument(
        "--controls-csv",
        help="Optional controls CSV with drug,target_id,smiles/status_for_workflow columns. Defaults to influenza known_drugs.csv plus local PubChem SMILES cache.",
    )
    stage4_parser.add_argument("--receptor-pdb", help="Optional local receptor PDB file. Used for receptor readiness and GNINA planning.")
    stage4_parser.add_argument("--pdb-id", default="", help="PDB id for the receptor package. Defaults to stage 2/catalog/config metadata.")
    stage4_parser.add_argument("--pocket-center", default="", help="Docking box center as 'x,y,z'. Stored in project target.pocket.center.")
    stage4_parser.add_argument("--pocket-size", default="", help="Docking box size as 'x,y,z'. Stored in project target.pocket.size.")
    stage4_parser.add_argument("--pocket-source", default="", help="Source label for the docking box coordinates.")
    stage4_parser.add_argument("--fetch-receptor", action="store_true", help="Fetch receptor PDB from RCSB when --pdb-id is available.")
    stage4_parser.add_argument("--top", type=int, default=10, help="Number of diverse RDKit-ready molecules to select.")
    stage4_parser.add_argument("--decoys", type=int, default=8, help="Number of local drug-like decoys to include in the benchmark panel.")
    stage4_parser.add_argument("--max-conformers", type=int, default=1, help="Embed one 3D conformer per ligand when greater than zero.")
    stage4_parser.add_argument("--seed", type=int, default=61453, help="RDKit conformer random seed.")
    stage4_parser.add_argument("--no-sdf", action="store_true", help="Skip ligand SDF generation.")
    stage4_parser.add_argument("--no-render-2d", action="store_false", dest="render_2d", help="Skip PNG rendering for molecule cards.")
    stage4_parser.set_defaults(render_2d=True)
    stage4_parser.add_argument("--docking-backend", choices=["auto", "vina", "gnina"], default="auto", help="Docking backend to plan or run.")
    stage4_parser.add_argument("--run-docking", action="store_true", help="Attempt real docking only when backend and receptor/ligand inputs are available.")
    stage4_parser.add_argument("--docking-timeout", type=int, default=600, help="Timeout in seconds for each docking subprocess.")
    stage4_parser.add_argument(
        "--rescore",
        action="store_true",
        help="After Stage 4 assets are written, rerun score/rank/feedback with the RDKit descriptors.",
    )
    stage4_parser.add_argument("--rank-top", type=int, default=0, help="Top candidates for rank when --rescore is used. Defaults to --top.")
    stage4_parser.add_argument("--feedback-top", type=int, default=0, help="Top candidates for feedback when --rescore is used. Defaults to --rank-top.")
    stage4_parser.add_argument(
        "--external-scores",
        help="Optional real docking/pose CSV passed to score during --rescore.",
    )
    stage4_parser.set_defaults(func=stage4_real)

    stage45_parser = sub.add_parser(
        "stage4-validate-controls",
        help="Run Stage 4.5 control calibration by docking candidates, known controls, and decoys in one pool.",
    )
    stage45_parser.add_argument("project")
    stage45_parser.add_argument("--round", type=int, default=1, help="Candidate round to calibrate.")
    stage45_parser.add_argument(
        "--target",
        default="",
        help="Target id used to select known controls, e.g. influenza_a_h1n1_na. Defaults to Stage 4 assets target_id.",
    )
    stage45_parser.add_argument("--top-candidates", type=int, default=12, help="Top ranked candidates to include in the calibration pool.")
    stage45_parser.add_argument("--decoys", type=int, default=8, help="Number of local drug-like decoys to include.")
    stage45_parser.add_argument("--controls-csv", help="Optional controls CSV. Defaults to influenza known_drugs.csv plus local SMILES cache.")
    stage45_parser.add_argument("--docking-backend", choices=["auto", "vina", "gnina"], default="auto", help="Docking backend to run.")
    stage45_parser.add_argument("--docking-timeout", type=int, default=600, help="Timeout in seconds for each docking subprocess.")
    stage45_parser.add_argument("--seed", type=int, default=61453, help="RDKit conformer random seed.")
    stage45_parser.add_argument("--no-docking", action="store_true", help="Write the control calibration plan without executing docking.")
    stage45_parser.set_defaults(func=stage45_validate_controls)

    stage46_parser = sub.add_parser(
        "stage4-retrospective-benchmark",
        help="Run Stage 4.6 retrospective benchmark from Stage 4.5 control/decoy docking scores.",
    )
    stage46_parser.add_argument("project")
    stage46_parser.add_argument("--round", type=int, default=1, help="Candidate round to benchmark.")
    stage46_parser.add_argument(
        "--positive-types",
        default=",".join(STAGE46_DEFAULT_POSITIVE_TYPES),
        help="Comma-separated panel_type values treated as known positives.",
    )
    stage46_parser.add_argument(
        "--negative-types",
        default=",".join(STAGE46_DEFAULT_NEGATIVE_TYPES),
        help="Comma-separated panel_type values treated as negatives.",
    )
    stage46_parser.add_argument("--top-k", default="1,3,5,10", help="Comma-separated Top-K cutoffs for control recovery.")
    stage46_parser.set_defaults(func=stage46_retrospective_benchmark)

    stage5_parser = sub.add_parser(
        "stage5-dashboard",
        help="Build a productized local static dashboard from Stage 1-4 project outputs.",
    )
    stage5_parser.add_argument("project")
    stage5_parser.add_argument("--round", type=int, default=4, help="Round to summarize in the dashboard.")
    stage5_parser.add_argument(
        "--title",
        default="",
        help="Optional dashboard title. Defaults to '<project> round <N> closed-loop dashboard'.",
    )
    stage5_parser.set_defaults(func=stage5_dashboard)

    stage6_parser = sub.add_parser(
        "stage6-validate",
        help="Build validation operations assets: quality gates, hit triage, assay queue, and risk register.",
    )
    stage6_parser.add_argument("project")
    stage6_parser.add_argument("--round", type=int, default=4, help="Round to validate operationally.")
    stage6_parser.add_argument("--top", type=int, default=8, help="Top advanced or ranked hits to triage.")
    stage6_parser.set_defaults(func=stage6_validate)

    stage7_parser = sub.add_parser(
        "stage7-package",
        help="Build delivery package, reproducibility runbook, demo checklist, and Stage 8 frontend spec.",
    )
    stage7_parser.add_argument("project")
    stage7_parser.add_argument("--round", type=int, default=4, help="Round to package for delivery.")
    stage7_parser.add_argument(
        "--title",
        default="",
        help="Optional delivery title. Defaults to '<project> round <N> delivery package'.",
    )
    stage7_parser.set_defaults(func=stage7_package)

    demo_parser = sub.add_parser("run-demo", help="Run init/generate/score/rank/feedback in one command.")
    demo_parser.add_argument("project")
    demo_parser.add_argument("--rounds", type=int, default=2)
    demo_parser.add_argument("--n", type=int, default=30)
    demo_parser.add_argument("--top", type=int, default=8)
    demo_parser.add_argument("--external-scores", help="Optional external scores CSV used in every round.")
    demo_parser.set_defaults(func=run_demo)

    doctor_parser = sub.add_parser("doctor", help="Check local repository paths and optional executables.")
    doctor_parser.add_argument("project", nargs="?")
    doctor_parser.set_defaults(func=doctor)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
