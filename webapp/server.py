#!/usr/bin/env python3
"""AI 分子设计闭环 Web 可视化平台 — FastAPI 后端

启动: python3 webapp/server.py
"""

from __future__ import annotations
import argparse, csv, json, queue, re, shutil, subprocess, sys, threading, traceback, uuid, zipfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ai_mol_loop"))

from ai_mol_loop import (
    DEFAULT_ROOT, DEFAULT_DESIGN_REPOS, DEFAULT_EVAL_TOOLS,
    DEFAULT_TARGET_CATALOG, DEFAULT_PDB_STRUCTURES, DEFAULT_KNOWN_DRUGS,
    DEFAULT_WEIGHTS, DEFAULT_SEEDS, SUBSTITUENTS, TEMPLATE_MOLECULES,
    ensure_project_dirs, read_csv, write_csv, read_json, write_json, write_text,
    load_config, plausible_smiles, estimate_scores, weighted_total,
    load_external_scores, float_field, render_round_summary, project_path,
    default_config, seed_file_for_generation, load_seed_smiles,
    aromatic_variant, chain_variant, generate_candidates, score_candidates,
    rank_candidates, feedback,
    load_target_catalog, sorted_targets, get_target_by_id, first_structure_asset,
    target_to_brief, extract_stage4_cocrystal_pocket, stage4_box_from_receptor_package,
    find_stage4_project_receptor, target_select, brief_from_target, evidence_stage2, stage3_screen, stage4_real, stage45_validate_controls,
    stage46_retrospective_benchmark, stage4_capabilities, stage5_dashboard,
    build_stage5_dashboard_data, stage6_validate, stage7_package,
    target_brief_from_args, render_generator_prompt,
    split_items, parse_triplet, PROMPT_EXAMPLES, stage45_export_reference_ligand_pdb,
)

from fastapi import FastAPI, HTTPException, Query, Response
from fastapi.responses import FileResponse, HTMLResponse, StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
import uvicorn

app = FastAPI(title="AI 分子设计闭环可视化平台", version="1.2.0")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

STATIC = ROOT / "webapp" / "static"
PROJECTS_ROOT = ROOT / "webapp" / "projects"
PROJECTS_ROOT.mkdir(parents=True, exist_ok=True)
app.mount("/static", StaticFiles(directory=str(STATIC)), name="static")

# ── 运行引擎 ──────────────────────────────────────────────────────────────
run_queues: Dict[str, queue.Queue] = {}
run_status: Dict[str, Dict[str, Any]] = {}
run_lock = threading.Lock()

job_status: Dict[str, Dict[str, Any]] = {}
job_lock = threading.Lock()

def emit(run_id: str, event: str, data: Any = None):
    with run_lock:
        q = run_queues.get(run_id)
        if q: q.put({"event": event, "data": data})

def start_run(run_id: str):
    with run_lock:
        run_queues[run_id] = queue.Queue()
        run_status[run_id] = {"status": "running", "started_at": datetime.now().isoformat(), "rounds": []}

def finish_run(run_id: str, error: Optional[str] = None):
    with run_lock:
        if run_id in run_status:
            run_status[run_id]["status"] = "error" if error else "completed"
            run_status[run_id]["finished_at"] = datetime.now().isoformat()
            if error: run_status[run_id]["error"] = error
        q = run_queues.get(run_id)
        if q: q.put({"event": "done", "data": {"error": error}})

# ── Models ────────────────────────────────────────────────────────────────
class InitProjectRequest(BaseModel):
    name: str; force: bool = False

class ConfigUpdateRequest(BaseModel):
    config: Dict[str, Any]

class RunRequest(BaseModel):
    rounds: int = 2; n: int = 24; top: int = 6
    external_scores: Optional[str] = None

class Stage1TargetSelectRequest(BaseModel):
    disease: str = "influenza_a"
    target: str = ""
    top: int = 5
    catalog: Optional[str] = None

class Stage1BriefFromTargetRequest(BaseModel):
    disease: str = "influenza_a"
    target: str = ""
    catalog: Optional[str] = None
    free_text: str = ""
    max_heavy_atoms: int = 55
    max_molecular_weight: float = 550.0
    force: bool = False

class TargetIntakeRequest(BaseModel):
    source_kind: str = "text"
    source: str = ""
    disease: str = "influenza_a"
    target_hint: str = ""
    target: str = ""
    round: int = 1

class TargetPackRequest(BaseModel):
    round: int = 1
    target: str = ""
    pdb_id: str = ""
    reference_ligand: str = ""
    target_hint: str = ""

class Stage4PocketPackRequest(BaseModel):
    round: int = 1
    target: str = ""
    pdb_id: str = ""
    reference_ligand: str = ""

class ExportPackageRequest(BaseModel):
    round: int = 1

class Stage8FullExportRequest(BaseModel):
    round: int = 1
    title: str = ""

class ProductJobRequest(BaseModel):
    task: str
    project: str = ""
    round: int = 1
    payload: Dict[str, Any] = {}

class Stage2EvidenceRequest(BaseModel):
    disease: str = "influenza"
    target: str = ""
    top: int = 5
    catalog: Optional[str] = None
    sources: Optional[str] = None
    evidence_dir: Optional[str] = None
    refresh: bool = False
    retmax: int = 5
    timeout: int = 20
    offline: bool = True

class Stage4RealRequest(BaseModel):
    round: int = 1
    target: str = ""
    input_csv: Optional[str] = None
    controls_csv: Optional[str] = None
    receptor_pdb: Optional[str] = None
    pdb_id: str = ""
    pocket_center: str = ""
    pocket_size: str = ""
    pocket_source: str = "stage4_manual"
    fetch_receptor: bool = False
    top: int = 10
    decoys: int = 8
    max_conformers: int = 1
    seed: int = 61453
    no_sdf: bool = False
    render_2d: bool = True
    docking_backend: str = "auto"
    run_docking: bool = False
    docking_timeout: int = 600
    rescore: bool = False
    rank_top: int = 0
    feedback_top: int = 0
    external_scores: Optional[str] = None

class Stage4SmokeTestRequest(BaseModel):
    timeout: int = 15

class Stage45ValidateRequest(BaseModel):
    round: int = 1
    target: str = ""
    top_candidates: int = 12
    decoys: int = 8
    controls_csv: Optional[str] = None
    docking_backend: str = "auto"
    docking_timeout: int = 600
    seed: int = 61453
    no_docking: bool = False

class Stage46BenchmarkRequest(BaseModel):
    round: int = 1
    positive_types: str = "positive_control,reference_control,control"
    negative_types: str = "decoy"
    top_k: str = "1,3,5,10"

class Stage5DashboardRequest(BaseModel):
    round: int = 1
    title: str = ""

class Stage6ValidateRequest(BaseModel):
    round: int = 1
    top: int = 8

class Stage7PackageRequest(BaseModel):
    round: int = 1
    title: str = ""

class Stage8DemoPackageRequest(BaseModel):
    round: int = 1
    top: int = 8
    title: str = ""

class Stage8DemoRunnerRequest(BaseModel):
    round: int = 1
    disease: str = "甲流"
    target: str = "influenza_a_h1n1_na"
    source_kind: str = "url"
    source: str = "https://www.rcsb.org/structure/3TI6"
    target_hint: str = "neuraminidase oseltamivir pocket"
    candidates: int = 24
    top: int = 6
    decoys: int = 8
    run_docking: bool = False
    docking_backend: str = "auto"
    docking_timeout: int = 600
    title: str = ""

class Stage8ReportRequest(BaseModel):
    round: int = 1
    title: str = ""

class Stage8EvidencePackRequest(BaseModel):
    round: int = 1
    title: str = ""

class Stage8RepairRequest(BaseModel):
    round: int = 1
    disease: str = "甲流"
    target: str = "influenza_a_h1n1_na"
    source_kind: str = "url"
    source: str = "https://www.rcsb.org/structure/3TI6"
    target_hint: str = "neuraminidase oseltamivir pocket"
    candidates: int = 24
    top: int = 6
    decoys: int = 8
    run_docking: bool = False
    docking_backend: str = "auto"
    docking_timeout: int = 600

class Stage8AcceptanceDemoRequest(BaseModel):
    project_name: str = "flu_na_acceptance_demo"
    force: bool = False
    round: int = 1
    candidates: int = 24
    top: int = 6
    decoys: int = 8
    run_docking: bool = False

class TargetCatalogCustomRequest(BaseModel):
    target_id: str
    display_name: str
    disease: str = ""
    pdb_id: str = ""
    reference_ligand: str = ""
    known_drugs: str = ""
    mechanism: str = ""
    recommendation: str = "custom_review"

class Stage3CandidateRequest(BaseModel):
    round: int = 1
    source_mode: str = "proxy"
    source_text: str = ""
    source_csv: Optional[str] = None
    source_url: Optional[List[str]] = None
    context_url: Optional[List[str]] = None
    use_openai: bool = False
    openai_model: str = "gpt-5.2"
    api_key: Optional[str] = None
    api_key_file: Optional[str] = None
    api_key_env: str = "OPENAI_API_KEY"
    prompt: str = ""
    n: int = 24
    top: int = 6
    timeout: int = 60
    max_url_bytes: int = 800000
    max_heavy_atoms: int = 55
    max_molecular_weight: float = 550.0
    max_lipinski_violations: int = 1
    allow_risk: bool = False
    external_scores: Optional[str] = None
    no_score: bool = False

class Stage4RepairRequest(BaseModel):
    round: int = 1
    persist: bool = True

class BriefRequest(BaseModel):
    target_name: str; disease: str = ""; rationale: str = ""; free_text: str = ""
    protein: str = ""; gene: str = ""; pdb_id: str = ""; protein_pdb: str = ""
    pocket: str = ""; reference_ligand: str = ""; key_residues: str = ""
    center: str = ""; size: str = ""; pocket_source: str = "manual_or_literature"
    style: str = "de_novo_and_analog_design"
    primary_goal: str = "improve predicted binding while retaining drug-like properties."
    activity: str = "in silico binding enrichment, not measured potency."
    must_have: str = ""; avoid: str = ""; desired_properties: str = ""; selectivity: str = ""
    max_heavy_atoms: int = 55; max_molecular_weight: float = 500.0
    preserve_scaffold: bool = False; diversity: bool = True; force: bool = False

# ── 辅助 ──────────────────────────────────────────────────────────────────
def get_project_dir(name: str) -> Path:
    p = PROJECTS_ROOT / name
    if not p.exists(): raise HTTPException(404, f"项目 {name} 不存在")
    return p

def project_latest_round(project_dir: Path) -> int:
    rounds = []
    candidates_dir = project_dir / "candidates"
    if candidates_dir.exists():
        for f in candidates_dir.glob("round_*_candidates.csv"):
            try:
                rounds.append(int(f.stem.split("_")[1]))
            except (IndexError, ValueError):
                continue
    return max(rounds) if rounds else 1

def project_asset_status(project_dir: Path, round_no: int) -> Dict[str, Dict[str, Any]]:
    candidates_path = project_dir / "candidates" / f"round_{round_no}_candidates.csv"
    candidate_rows = csv_to_dicts(candidates_path)
    candidate_count = len(candidate_rows)
    candidate_status = "ready" if candidate_count else "missing"

    stage4_dir = project_dir / "stage4"
    receptor_package_path = stage4_dir / f"round_{round_no}_receptor_package.json"
    receptor_package = read_json_if_exists(receptor_package_path)
    receptor_pdb = stage4_existing_path(receptor_package.get("local_receptor_pdb", ""))
    receptor_pdbqt = stage4_existing_path(receptor_package.get("local_receptor_pdbqt", ""))
    if not receptor_pdb:
        receptor_pdb = stage4_existing_path(receptor_package.get("prepared_receptor_pdb", ""))
    if not receptor_pdbqt:
        receptor_pdbqt = stage4_existing_path(receptor_package.get("prepared_receptor_pdbqt", ""))
    receptor_status = "ready" if receptor_pdb or receptor_pdbqt else ("partial" if receptor_package else "missing")

    docking_plan_path = stage4_dir / f"round_{round_no}_docking_plan.json"
    docking_plan = read_json_if_exists(docking_plan_path)
    docking_scores_path = stage4_dir / f"round_{round_no}_docking_scores_template.csv"
    docking_rows = csv_to_dicts(docking_scores_path)
    scored_rows = [row for row in docking_rows if str(row.get("docking_score") or row.get("score") or row.get("affinity") or "").strip()]
    plan_status = str(docking_plan.get("status") or "").strip() or ("completed" if scored_rows else "")
    if scored_rows:
        docking_status = "completed" if plan_status in {"completed", "imported_external_scores"} else "scored"
    elif plan_status:
        docking_status = plan_status
    else:
        docking_status = "missing"

    return {
        "candidate": {
            "status": candidate_status,
            "count": candidate_count,
            "round": round_no,
            "file": str(candidates_path) if candidates_path.exists() else "",
        },
        "receptor": {
            "status": receptor_status,
            "round": round_no,
            "pdb_id": str(receptor_package.get("pdb_id") or receptor_package.get("structure_id") or ""),
            "target_id": str(receptor_package.get("target_id") or ""),
            "receptor_pdb": receptor_pdb,
            "receptor_pdbqt": receptor_pdbqt,
            "file": str(receptor_package_path) if receptor_package_path.exists() else "",
        },
        "docking": {
            "status": docking_status,
            "round": round_no,
            "plan_status": plan_status or "missing",
            "backend": str(docking_plan.get("selected_backend") or docking_plan.get("backend") or (scored_rows[0].get("backend") if scored_rows else "") or ""),
            "score_count": len(scored_rows),
            "file": str(docking_scores_path) if docking_scores_path.exists() else "",
        },
    }

def list_projects() -> List[Dict]:
    projects = []
    for d in sorted(PROJECTS_ROOT.iterdir()):
        if d.is_dir():
            config_path = d / "config.json"
            cfg = {}
            if config_path.exists():
                try: cfg = json.loads(config_path.read_text())
                except Exception: pass
            has_brief = (d / "briefs" / "target_brief.json").exists()
            latest_round = project_latest_round(d)
            projects.append({
                "name": d.name, "has_config": config_path.exists(),
                "has_brief": has_brief,
                "target": cfg.get("target", {}).get("name", "") if isinstance(cfg.get("target"), dict) else "",
                "latest_round": latest_round,
                "asset_status": project_asset_status(d, latest_round),
                "created": datetime.fromtimestamp(d.stat().st_ctime).isoformat(),
            })
    return projects

def get_round_files(project_dir: Path) -> List[Dict]:
    rounds = []
    candidates_dir = project_dir / "candidates"
    if candidates_dir.exists():
        for f in sorted(candidates_dir.glob("round_*_candidates.csv")):
            try: rn = int(f.stem.split("_")[1])
            except (IndexError, ValueError): continue
            rounds.append({
                "round": rn, "candidates": f.name,
                "scores": (project_dir / "scores" / f"round_{rn}_scores.csv").exists(),
                "ranked": (project_dir / "ranked" / f"round_{rn}_ranked.csv").exists(),
                "feedback": (project_dir / "feedback" / f"round_{rn}_feedback.json").exists(),
                "report": (project_dir / "reports" / f"round_{rn}_summary.md").exists(),
            })
    return rounds

def csv_to_dicts(path: Path) -> List[Dict]:
    if not path.exists(): return []
    return read_csv(path)

def read_text_if_exists(path: Path) -> str:
    return path.read_text(encoding="utf-8") if path.exists() else ""

def read_json_if_exists(path: Path) -> Dict[str, Any]:
    if not path.exists(): return {}
    try:
        data = read_json(path)
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}

def http_from_cli_exit(exc: SystemExit):
    message = str(exc) or "CLI command failed"
    status = 409 if "already exists" in message else 400
    raise HTTPException(status, message)

def safe_project_artifact_path(project_dir: Path, raw_path: str) -> Path:
    if not str(raw_path or "").strip():
        raise HTTPException(400, "下载路径不能为空")
    project_root = project_dir.resolve()
    requested = Path(str(raw_path))
    candidate = requested if requested.is_absolute() else project_dir / requested
    try:
        resolved = candidate.resolve()
        resolved.relative_to(project_root)
    except ValueError:
        raise HTTPException(400, "下载路径必须位于项目目录内")
    except OSError:
        raise HTTPException(400, "下载路径无效")
    if not resolved.exists() or not resolved.is_file():
        raise HTTPException(404, "项目产物不存在")
    return resolved

def project_artifact_url(name: str, project_dir: Path, path: Path) -> str:
    try:
        rel_path = path.resolve().relative_to(project_dir.resolve()).as_posix()
    except ValueError:
        return ""
    return f"/api/projects/{quote(name, safe='')}/artifact?path={quote(rel_path, safe='')}"

def stage8_download_links(name: str, project_dir: Path, round_no: int) -> Dict[str, str]:
    files = {
        "stage5_dashboard": project_dir / "stage5" / "index.html",
        "stage5_dashboard_data": project_dir / "stage5" / "dashboard_data.json",
        "stage5_report": project_dir / "reports" / "stage5_dashboard_report.md",
        "stage6_assets": project_dir / "stage6" / f"round_{round_no}_validation_assets.json",
        "stage6_quality_gates": project_dir / "stage6" / f"round_{round_no}_quality_gates.csv",
        "stage6_hit_triage": project_dir / "stage6" / f"round_{round_no}_hit_triage.csv",
        "stage6_assay_queue": project_dir / "stage6" / f"round_{round_no}_assay_queue.csv",
        "stage7_manifest": project_dir / "stage7" / f"round_{round_no}_delivery_manifest.json",
        "executive_summary": project_dir / "stage7" / f"round_{round_no}_executive_summary.md",
        "reproducibility_runbook": project_dir / "stage7" / f"round_{round_no}_reproducibility.md",
        "demo_checklist": project_dir / "stage7" / f"round_{round_no}_investor_demo_checklist.csv",
        "stage8_frontend_spec": project_dir / "stage7" / "stage8_frontend_product_spec.md",
        "stage7_report": project_dir / "reports" / f"stage7_round_{round_no}_delivery_report.md",
        "ranked_candidates": project_dir / "ranked" / f"round_{round_no}_ranked.csv",
        "candidate_scores": project_dir / "scores" / f"round_{round_no}_scores.csv",
        "raw_candidates": project_dir / "candidates" / f"round_{round_no}_candidates.csv",
    }
    links: Dict[str, str] = {}
    for key, path in files.items():
        if path.exists() and path.is_file():
            url = project_artifact_url(name, project_dir, path)
            if url:
                links[key] = url
    return links

def safe_name(raw: str) -> str:
    return re.sub(r"[^a-zA-Z0-9_\-.一-龥]", "_", raw)

def stage4_image_records(project_name: str, project_dir: Path, round_no: int, assets: Dict[str, Any]) -> List[Dict[str, str]]:
    records: List[Dict[str, str]] = []
    image_dir = project_dir / "stage4" / f"round_{round_no}_2d"
    rendered = assets.get("rendered_2d_images", []) if isinstance(assets, dict) else []
    if isinstance(rendered, list):
        for item in rendered:
            if not isinstance(item, dict): continue
            raw_path = item.get("path", "")
            filename = Path(str(raw_path)).name
            if not filename: continue
            local = image_dir / filename
            if not local.exists(): continue
            records.append({
                "panel_type": str(item.get("panel_type", "")),
                "id": str(item.get("id", "")),
                "filename": filename,
                "url": f"/api/projects/{project_name}/stage4/images/{round_no}/{filename}",
            })
    if not records and image_dir.exists():
        for img in sorted(image_dir.glob("*.png")):
            records.append({
                "panel_type": "",
                "id": img.stem,
                "filename": img.name,
                "url": f"/api/projects/{project_name}/stage4/images/{round_no}/{img.name}",
            })
    return records

def stage4_target_presets() -> List[Dict[str, Any]]:
    catalog = read_json_if_exists(DEFAULT_TARGET_CATALOG)
    targets = catalog.get("targets", []) if isinstance(catalog, dict) else []
    structures = csv_to_dicts(DEFAULT_PDB_STRUCTURES)
    drugs = csv_to_dicts(DEFAULT_KNOWN_DRUGS)
    by_target_structures: Dict[str, List[Dict[str, str]]] = {}
    for row in structures:
        by_target_structures.setdefault(str(row.get("target_id", "")), []).append(row)
    by_target_drugs: Dict[str, List[str]] = {}
    for row in drugs:
        if row.get("status_for_workflow") in {"positive_control", "reference_control"}:
            by_target_drugs.setdefault(str(row.get("target_id", "")), []).append(str(row.get("drug", "")))
    presets: List[Dict[str, Any]] = []
    for target in targets:
        if not isinstance(target, dict):
            continue
        tid = str(target.get("id", ""))
        structures_for_target = by_target_structures.get(tid, [])
        recommended = ""
        for item in structures_for_target:
            if "recommended" in str(item.get("use", "")):
                recommended = str(item.get("pdb_id", ""))
                break
        if not recommended and structures_for_target:
            recommended = str(structures_for_target[0].get("pdb_id", ""))
        site = target.get("binding_site", {})
        site = site if isinstance(site, dict) else {}
        presets.append(
            {
                "target_id": tid,
                "label": str(target.get("display_name") or tid),
                "short_name": str(target.get("short_name") or tid),
                "recommendation": str(target.get("recommendation", "")),
                "recommended_pdb": recommended,
                "positive_controls": by_target_drugs.get(tid, []),
                "reference_ligand": str(site.get("reference_ligand", "")),
                "key_residues": site.get("key_residues", []),
                "pocket": {
                    "status": "needs_curated_coordinates",
                    "center": [],
                    "size": [],
                    "source": str(site.get("source", "")),
                    "strategy": str(site.get("box_strategy", "")),
                    "description": str(site.get("description", "")),
                },
                "structures": structures_for_target,
                "boundary": "该预设提供靶点、控药和结构建议；真实 docking box 需来自共晶配体质心、文献或人工结构准备。",
            }
        )
    smoke_receptor = DEFAULT_EVAL_TOOLS / "AutoDock-Vina" / "example" / "basic_docking" / "solution" / "1iep_receptor.pdbqt"
    presets.append(
        {
            "target_id": "vina_official_1iep_smoke",
            "label": "Vina 官方 1IEP 工具链烟测",
            "short_name": "Vina smoke",
            "recommendation": "toolchain_smoke_test",
            "recommended_pdb": "1IEP",
            "positive_controls": [],
            "reference_ligand": "1iep_ligand",
            "key_residues": [],
            "local_receptor_pdbqt": str(smoke_receptor),
            "pocket": {
                "status": "ready",
                "center": [15.19, 53.903, 16.917],
                "size": [20.0, 20.0, 20.0],
                "source": "autodock_vina_official_example",
                "strategy": "Use only to verify local Vina execution; not a disease-target claim.",
                "description": "AutoDock Vina bundled 1IEP example search box.",
            },
            "structures": [],
            "boundary": "该预设只用于验证本地 docking 工具链，不代表甲流或其他疾病靶点。",
        }
    )
    return presets

def stage4_score_band(score: Any) -> str:
    try:
        value = float(score)
    except Exception:
        return "missing"
    if value <= -7.5:
        return "strong"
    if value <= -6.0:
        return "moderate"
    return "weak"

def stage4_note_value(notes: str, key: str) -> str:
    match = re.search(r"(?:^|;)" + re.escape(key) + r"=([^;]+)", notes or "")
    return match.group(1) if match else ""

def stage4_existing_path(path: str) -> str:
    if not path:
        return ""
    candidate = Path(path)
    return str(candidate) if candidate.exists() else ""

def stage4_path_exists(path: Any) -> bool:
    raw = str(path or "")
    return bool(raw and Path(raw).exists())

def stage4_find_receptor_files(project_dir: Path, round_no: int, rec: Dict[str, Any], stage45_inputs: List[Dict[str, Any]], stage45_scores: List[Dict[str, Any]]) -> Dict[str, str]:
    receptor_dir = project_dir / "stage4" / "receptors"
    pdb_candidates: List[str] = []
    pdbqt_candidates: List[str] = []
    for key in ["local_receptor_pdb", "prepared_receptor_pdb", "receptor_pdb"]:
        value = rec.get(key)
        if value:
            pdb_candidates.append(str(value))
    for key in ["local_receptor_pdbqt", "prepared_receptor_pdbqt", "receptor_pdbqt"]:
        value = rec.get(key)
        if value:
            pdbqt_candidates.append(str(value))
    for row in stage45_inputs:
        if not isinstance(row, dict):
            continue
        if row.get("receptor_pdb"):
            pdb_candidates.append(str(row.get("receptor_pdb")))
        if row.get("receptor_pdbqt"):
            pdbqt_candidates.append(str(row.get("receptor_pdbqt")))
    for row in stage45_scores:
        if not isinstance(row, dict):
            continue
        receptor = str(row.get("receptor") or "")
        if receptor.endswith(".pdbqt"):
            pdbqt_candidates.append(receptor)
        elif receptor.endswith(".pdb"):
            pdb_candidates.append(receptor)
    pdb_id = str(rec.get("pdb_id") or "").strip()
    if pdb_id and receptor_dir.exists():
        pdb_candidates.extend(str(path) for path in sorted(receptor_dir.glob(f"{pdb_id}*.pdb")))
        pdbqt_candidates.extend(str(path) for path in sorted(receptor_dir.glob(f"{pdb_id}*.pdbqt")))
    if receptor_dir.exists():
        pdb_candidates.extend(str(path) for path in sorted(receptor_dir.glob("*.pdb")))
        pdbqt_candidates.extend(str(path) for path in sorted(receptor_dir.glob("*.pdbqt")))

    def first_existing(paths: List[str]) -> str:
        for raw in paths:
            path = Path(str(raw))
            if path.exists() and path.is_file():
                return str(path)
        return ""

    return {"pdb": first_existing(pdb_candidates), "pdbqt": first_existing(pdbqt_candidates)}

def stage4_result_record(
    row: Dict[str, Any],
    run_dir: Path,
    evidence_source: str,
) -> Optional[Dict[str, Any]]:
    score = row.get("docking_score", "") or row.get("score", "") or row.get("affinity", "") or row.get("vina_score", "")
    if score in {"", None}:
        return None
    item_id = str(row.get("id", ""))
    notes = str(row.get("notes", ""))
    log_file = stage4_existing_path(stage4_note_value(notes, "log")) or stage4_existing_path(str(run_dir / f"{item_id}.log"))
    posebusters_report = (
        stage4_existing_path(stage4_note_value(notes, "posebusters_report"))
        or stage4_existing_path(str(run_dir / f"{item_id}_posebusters.csv"))
    )
    pose_pdbqt = stage4_existing_path(str(run_dir / f"{item_id}_pose.pdbqt"))
    pose_sdf = stage4_existing_path(str(run_dir / f"{item_id}_pose.sdf"))
    pose_pass = str(row.get("pose_pass", "")).lower() in {"true", "pass", "passed", "yes", "1"}
    band = str(row.get("score_band") or stage4_score_band(score))
    return {
        **row,
        "score": score,
        "docking_score": score,
        "score_band": band,
        "pose_pass_bool": pose_pass,
        "log_file": log_file,
        "pose_file": pose_pdbqt or pose_sdf,
        "pose_sdf": pose_sdf,
        "posebusters_report": posebusters_report,
        "evidence_source": evidence_source,
        "interpretation": {
            "score_band": band,
            "pose_quality": "passed" if pose_pass else "not_passed_or_not_checked",
            "summary": (
                "Docking 分数较强，且 PoseBusters 通过；可进入人工复核。"
                if band == "strong" and pose_pass
                else "Docking 有可读分数；仍需结合 pose、控药和化学可行性复核。"
                if band in {"strong", "moderate"}
                else "Docking 信号偏弱或缺失，优先级较低。"
            ),
        },
    }

def stage4_docking_results(project_dir: Path, round_no: int) -> List[Dict[str, Any]]:
    stage4_dir = project_dir / "stage4"
    scores_path = stage4_dir / f"round_{round_no}_docking_scores_template.csv"
    rows = csv_to_dicts(scores_path)
    run_dir = stage4_dir / "docking_runs"
    results: List[Dict[str, Any]] = []
    for row in rows:
        record = stage4_result_record(row, run_dir, "stage4")
        if record:
            results.append(record)
    return sorted(results, key=lambda item: float(item.get("docking_score", 999) or 999))

def stage45_docking_results(project_dir: Path, round_no: int) -> List[Dict[str, Any]]:
    stage45_dir = project_dir / "stage4_5"
    rows = csv_to_dicts(stage45_dir / f"round_{round_no}_control_docking_scores.csv")
    run_dir = stage45_dir / "control_docking_runs"
    results: List[Dict[str, Any]] = []
    for row in rows:
        record = stage4_result_record(row, run_dir, "stage4_5")
        if record:
            results.append(record)
    return sorted(results, key=lambda item: float(item.get("docking_score", 999) or 999))

def stage4_recoverable_docking_results(results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    recoverable = []
    for row in results:
        pose_file = str(row.get("pose_file") or row.get("pose_sdf") or "")
        notes = str(row.get("notes") or "")
        if pose_file and Path(pose_file).exists() and "missing_pose" not in notes and "missing_pose_output" not in notes:
            recoverable.append(row)
    return recoverable

def stage4_invalid_backend_reason(project_dir: Path, round_no: int) -> str:
    plan = read_json_if_exists(project_dir / "stage4" / f"round_{round_no}_docking_plan.json")
    if not isinstance(plan, dict):
        return ""
    tool_status = plan.get("tool_status", {}) if isinstance(plan.get("tool_status"), dict) else {}
    backend_candidates = [
        str(plan.get("selected_backend") or ""),
        str(plan.get("requested_backend") or ""),
        "vina",
        "gnina",
    ]
    seen = set()
    for backend in backend_candidates:
        if not backend or backend in seen:
            continue
        seen.add(backend)
        tool = tool_status.get(backend, {}) if isinstance(tool_status.get(backend, {}), dict) else {}
        if tool.get("status") == "invalid":
            validation = str(tool.get("validation") or "invalid_backend")
            return f"{backend}_invalid:{validation}"
    return ""

def stage4_apply_recovered_docking(payload: Dict[str, Any], project_dir: Path, round_no: int) -> Dict[str, Any]:
    rec = payload.get("receptor_package", {}) if isinstance(payload.get("receptor_package"), dict) else {}
    plan = payload.get("docking_plan", {}) if isinstance(payload.get("docking_plan"), dict) else {}
    validation = payload.get("validation_metrics", {}) if isinstance(payload.get("validation_metrics"), dict) else {}
    changed_metadata = False
    tool_status = plan.get("tool_status", {}) if isinstance(plan.get("tool_status"), dict) else {}
    selected_backend = str(plan.get("selected_backend") or "vina")
    selected_tool = tool_status.get(selected_backend, {}) if isinstance(tool_status.get(selected_backend, {}), dict) else {}
    recovery_block_reason = ""
    if selected_tool.get("status") == "invalid":
        recovery_block_reason = f"{selected_backend}_invalid:{selected_tool.get('validation', '')}"
    stage45_inputs = csv_to_dicts(project_dir / "stage4_5" / f"round_{round_no}_control_docking_inputs.csv")
    stage45_scores = csv_to_dicts(project_dir / "stage4_5" / f"round_{round_no}_control_docking_scores.csv")
    stage45_validation = read_json_if_exists(project_dir / "stage4_5" / f"round_{round_no}_control_validation.json")
    stage45_docking = stage45_validation.get("docking", {}) if isinstance(stage45_validation.get("docking"), dict) else {}
    stage45_results = stage45_docking_results(project_dir, round_no)
    recovered_results = [] if recovery_block_reason else stage4_recoverable_docking_results(stage45_results)
    receptor_files = stage4_find_receptor_files(project_dir, round_no, rec, stage45_inputs, stage45_scores)

    if receptor_files.get("pdb") and not stage4_path_exists(rec.get("local_receptor_pdb")):
        rec["local_receptor_pdb"] = receptor_files["pdb"]
        changed_metadata = True
    if receptor_files.get("pdbqt") and not stage4_path_exists(rec.get("local_receptor_pdbqt")):
        rec["local_receptor_pdbqt"] = receptor_files["pdbqt"]
        changed_metadata = True
    if receptor_files.get("pdb") or receptor_files.get("pdbqt"):
        ready_receptor_statuses = {
            "prepared_local_receptor",
            "receptor_pdbqt_available",
            "prepared_receptor",
            "ready",
        }
        current_status = str(rec.get("preparation_status") or "")
        if current_status not in ready_receptor_statuses:
            changed_metadata = True
            rec["preparation_status"] = "prepared_local_receptor"
        if not rec.get("recovered_from") and (stage45_inputs or stage45_scores) and current_status not in ready_receptor_statuses:
            rec["recovered_from"] = "stage4_5"

    box = plan.get("docking_box", {}) if isinstance(plan.get("docking_box"), dict) else {}
    stage45_box = stage45_docking.get("docking_box", {}) if isinstance(stage45_docking.get("docking_box"), dict) else {}
    if stage45_box.get("status") == "ready" and box.get("status") != "ready":
        plan["docking_box"] = stage45_box
        changed_metadata = True

    if recovered_results and str(plan.get("status", "")).lower() != "completed":
        plan["status"] = "completed"
        plan["selected_backend"] = str(stage45_docking.get("backend") or plan.get("selected_backend") or "vina")
        plan["run_docking_requested"] = True
        plan["recovered_from"] = "stage4_5"
        plan["source_scores_csv"] = str(project_dir / "stage4_5" / f"round_{round_no}_control_docking_scores.csv")
        plan.setdefault("expected_scores_csv", str(project_dir / "stage4" / f"round_{round_no}_docking_scores_template.csv"))
        plan["message"] = "Recovered completed docking evidence from Stage 4.5 control-calibration artifacts."
        changed_metadata = True

    readiness = validation.get("readiness", {}) if isinstance(validation.get("readiness"), dict) else {}
    if recovered_results:
        if validation.get("docking_status") != "completed" or readiness.get("docking") != "completed":
            changed_metadata = True
        validation["docking_status"] = "completed"
        validation["docking_backend"] = str(stage45_docking.get("backend") or plan.get("selected_backend") or "vina")
        readiness["docking"] = "completed"
        validation["readiness"] = readiness

    if recovered_results and not payload.get("docking_results"):
        payload["docking_results"] = recovered_results
    payload["receptor_package"] = rec
    payload["docking_plan"] = plan
    payload["validation_metrics"] = validation
    recovered_status = "none"
    if recovered_results:
        recovered_status = "recovered" if changed_metadata else "consistent"
    payload["recovered_docking_evidence"] = {
        "status": recovered_status,
        "source": "stage4_5" if recovered_results else "",
        "results": len(recovered_results),
        "available_stage45_scores": len(stage45_results),
        "ignored_stage45_scores": max(0, len(stage45_results) - len(recovered_results)),
        "receptor_pdb": rec.get("local_receptor_pdb", ""),
        "receptor_pdbqt": rec.get("local_receptor_pdbqt", ""),
        "docking_box": plan.get("docking_box", {}),
        "metadata_changed": changed_metadata,
        "note": (
            f"Stage 4.5 scores were ignored because current docking backend is invalid ({recovery_block_reason})."
            if recovery_block_reason and stage45_results
            else "" if recovered_results or not stage45_results else "Stage 4.5 scores were ignored because no pose files were available."
        ),
    }
    return payload

def stage4_visualization_file(project_dir: Path, name: str, path_value: Any, preferred_format: str = "") -> Dict[str, Any]:
    raw = str(path_value or "")
    if not raw:
        return {}
    path = Path(raw)
    if not path.exists() or not path.is_file():
        return {}
    url = project_artifact_url(name, project_dir, path)
    if not url:
        return {}
    suffix = path.suffix.lower().lstrip(".")
    return {
        "path": str(path),
        "url": url,
        "format": preferred_format or suffix,
        "filename": path.name,
    }

def stage4_reference_ligand_visualization(name: str, project_dir: Path, rec: Dict[str, Any]) -> Dict[str, Any]:
    export_rec = dict(rec)
    site = dict(export_rec.get("binding_site", {}) if isinstance(export_rec.get("binding_site"), dict) else {})
    if not isinstance(site.get("detected_ligand"), dict):
        receptor_pdb = Path(str(export_rec.get("local_receptor_pdb", "")))
        extracted = extract_stage4_cocrystal_pocket(receptor_pdb, str(site.get("reference_ligand", ""))) if receptor_pdb.exists() else {}
        detected = extracted.get("detected_ligand", {}) if isinstance(extracted, dict) else {}
        if isinstance(detected, dict) and detected:
            site["detected_ligand"] = detected
            site.setdefault("reference_ligand", extracted.get("reference_ligand", ""))
            export_rec["binding_site"] = site
    export = stage45_export_reference_ligand_pdb(export_rec, project_dir / "stage4" / "reference_ligands")
    if export.get("status") != "reference_pose_exported":
        return {}
    ligand = stage4_visualization_file(project_dir, name, export.get("path"), "pdb")
    if not ligand:
        return {}
    ligand.update(
        {
            "role": "co_crystal_reference_ligand",
            "label": "共晶参考配体",
            "resname": export.get("resname", ""),
            "chain": export.get("chain", ""),
            "resseq": export.get("resseq", ""),
            "atom_count": export.get("atom_count", ""),
            "boundary": "This is the co-crystallized reference ligand from the receptor structure, not a candidate docking pose.",
        }
    )
    return ligand

def stage4_visualization_assets(name: str, project_dir: Path, payload: Dict[str, Any]) -> Dict[str, Any]:
    rec = payload.get("receptor_package", {}) if isinstance(payload.get("receptor_package"), dict) else {}
    plan = payload.get("docking_plan", {}) if isinstance(payload.get("docking_plan"), dict) else {}
    box = plan.get("docking_box", {}) if isinstance(plan.get("docking_box"), dict) else {}
    results = payload.get("docking_results", []) if isinstance(payload.get("docking_results"), list) else []

    receptor = (
        stage4_visualization_file(project_dir, name, rec.get("local_receptor_pdb"), "pdb")
        or stage4_visualization_file(project_dir, name, rec.get("prepared_receptor_pdb"), "pdb")
        or stage4_visualization_file(project_dir, name, rec.get("local_receptor_pdbqt"), "pdbqt")
    )
    poses: List[Dict[str, Any]] = []
    for row in results:
        if not isinstance(row, dict):
            continue
        pose = (
            stage4_visualization_file(project_dir, name, row.get("pose_sdf"), "sdf")
            or stage4_visualization_file(project_dir, name, row.get("pose_file"), Path(str(row.get("pose_file", ""))).suffix.lower().lstrip("."))
        )
        if not pose:
            continue
        pose.update(
            {
                "id": str(row.get("id", "")),
                "smiles": str(row.get("smiles", "")),
                "docking_score": row.get("docking_score") or row.get("score") or "",
                "score_band": row.get("score_band") or (row.get("interpretation", {}) or {}).get("score_band", ""),
                "pose_pass": row.get("pose_pass", ""),
                "backend": row.get("backend", ""),
            }
        )
        poses.append(pose)

    status = "ready" if receptor and poses else ("missing_pose" if receptor else "missing_receptor")
    return {
        "viewer": "3Dmol.js",
        "status": status,
        "receptor": receptor,
        "reference_ligand": stage4_reference_ligand_visualization(name, project_dir, rec) if receptor else {},
        "top_pose": poses[0] if poses else {},
        "poses": poses[:25],
        "docking_box": {
            "status": box.get("status", ""),
            "center": box.get("center", []),
            "size": box.get("size", []),
            "source": box.get("source", ""),
        },
        "boundary": "3D viewer displays computational receptor/pose artifacts only; it is not experimental validation.",
    }

def stage4_readiness_guide(payload: Dict[str, Any]) -> Dict[str, Any]:
    rec = payload.get("receptor_package", {}) if isinstance(payload.get("receptor_package"), dict) else {}
    plan = payload.get("docking_plan", {}) if isinstance(payload.get("docking_plan"), dict) else {}
    box = plan.get("docking_box", {}) if isinstance(plan.get("docking_box"), dict) else {}
    files = payload.get("files", {}) if isinstance(payload.get("files"), dict) else {}
    docking = str(plan.get("status", "missing") or "missing")
    actions: List[str] = []
    if not payload.get("has_descriptors"):
        actions.append("先运行 Stage 4 生成 RDKit 描述符和候选 SDF。")
    candidate_file = str(files.get("candidates") or "")
    if candidate_file and not Path(candidate_file).exists():
        actions.append("先生成 candidates/round_N_candidates.csv，再运行 Stage 4。")
    if not rec.get("local_receptor_pdbqt") and not rec.get("local_receptor_pdb"):
        actions.append("提供 receptor PDB/PDBQT，或选择可拉取的 PDB ID。")
    if box.get("status") != "ready":
        actions.append("补入口袋中心和盒子尺寸，来源应来自共晶配体、文献或人工结构准备。")
    if docking in {"not_available", "blocked_missing_receptor", "blocked_missing_box", "attempted_no_scores", "missing", ""}:
        actions.append("确认 Vina/OpenBabel/Meeko/PoseBusters 可用后再勾选真实 docking。")
    if docking == "completed":
        actions.append("查看 docking 结果解释、pose 文件和日志，再决定是否进入下一轮反馈。")
    return {
        "descriptors": "ready" if payload.get("has_descriptors") else "missing",
        "receptor": "ready" if rec.get("local_receptor_pdbqt") or rec.get("local_receptor_pdb") else "missing",
        "box": str(box.get("status", "missing") or "missing"),
        "docking": docking,
        "next_actions": actions,
    }

def stage4_preflight(payload: Dict[str, Any]) -> List[Dict[str, Any]]:
    rec = payload.get("receptor_package", {}) if isinstance(payload.get("receptor_package"), dict) else {}
    plan = payload.get("docking_plan", {}) if isinstance(payload.get("docking_plan"), dict) else {}
    box = plan.get("docking_box", {}) if isinstance(plan.get("docking_box"), dict) else {}
    results = payload.get("docking_results", [])
    files = payload.get("files", {}) if isinstance(payload.get("files"), dict) else {}
    candidate_file = str(files.get("candidates") or "")
    candidate_exists = bool(candidate_file and Path(candidate_file).exists())
    candidate_rows = read_csv(Path(candidate_file)) if candidate_exists else []
    inputs = payload.get("docking_inputs", [])
    ready_ligands = [row for row in inputs if str(row.get("status", "")) == "ready_for_docking"]
    ligand_ready_count = len(ready_ligands)
    candidate_count = len([row for row in inputs if str(row.get("panel_type", "")) == "candidate"])
    plan_status = str(plan.get("status", "") or "")
    docking_completed = plan_status in {"completed", "imported_external_scores"} and bool(results)
    selected_backend = str(plan.get("selected_backend", "") or plan.get("requested_backend", "") or "")
    receptor_file = str(rec.get("local_receptor_pdbqt") or rec.get("local_receptor_pdb") or "")
    receptor_exists = bool(receptor_file and Path(receptor_file).exists())

    def item(
        step_id: str,
        label: str,
        status: str,
        evidence: str,
        next_action: str = "",
        auto_repair: bool = False,
        repair_hint: str = "",
    ) -> Dict[str, Any]:
        return {
            "step_id": step_id,
            "label": label,
            "status": status,
            "evidence": evidence,
            "next_action": next_action,
            "auto_repair": bool(auto_repair),
            "repair_hint": repair_hint,
        }

    tool_status = plan.get("tool_status", {}) if isinstance(plan.get("tool_status"), dict) else {}
    backend_tool = tool_status.get(selected_backend, {}) if isinstance(tool_status.get(selected_backend, {}), dict) else {}
    backend_ready = (
        plan_status in {"planned", "completed", "imported_external_scores"}
        or bool(selected_backend and selected_backend not in {"auto", "not_available"})
    ) and backend_tool.get("status") != "invalid"
    score_file = str(plan.get("expected_scores_csv") or files.get("docking_scores_template") or "")
    scores_ready = bool(results)
    return [
        item(
            "candidate_input",
            "候选输入文件",
            "ready" if candidate_exists and candidate_rows else ("empty" if candidate_exists else "missing"),
            candidate_file or "no candidate csv",
            "点击“修复Stage 4资产”（Stage 4修复）自动生成候选；或先到 Stage 3 生成 candidates/round_N_candidates.csv。"
            if not candidate_exists or not candidate_rows
            else "",
            True if not candidate_exists or not candidate_rows else False,
            "stage3_candidates",
        ),
        item(
            "candidate_assets",
            "候选与RDKit资产",
            "ready" if payload.get("has_descriptors") else "missing",
            f"descriptors={payload.get('has_descriptors')}; candidates={len(payload.get('descriptors', []))}; source_rows={len(candidate_rows)}",
            (
                "点击“修复Stage 4资产”可先补候选；随后运行 Stage 4 生成 RDKit 描述符。"
                if not candidate_exists or not candidate_rows
                else "候选文件已就绪；运行 Stage 4 生成 RDKit 描述符、SDF、控药/Decoy 和 docking plan。"
            )
            if not payload.get("has_descriptors")
            else "",
            False,
            "stage4_real_assets",
        ),
        item(
            "receptor",
            "受体结构",
            "ready" if receptor_exists else ("configured_missing_file" if receptor_file else "missing"),
            receptor_file or "no local receptor path",
            "补充真实存在的 receptor PDB/PDBQT 文件，或输入 PDB ID 并勾选拉取PDB受体。" if not receptor_exists else "",
            False,
            "receptor",
        ),
        item(
            "docking_box",
            "Docking box",
            "ready" if box.get("status") == "ready" else str(box.get("status", "missing") or "missing"),
            f"center={box.get('center', [])}; size={box.get('size', [])}; source={box.get('source', '')}",
            "补入口袋中心和盒子尺寸，并记录来源。" if box.get("status") != "ready" else "",
            False,
            "docking_box",
        ),
        item(
            "ligand_preparation",
            "配体PDBQT准备",
            "ready" if ligand_ready_count or docking_completed else ("pending" if candidate_count else "missing"),
            f"ready={ligand_ready_count}; candidates={candidate_count}; completed_docking_results={len(results) if docking_completed else 0}",
            "安装 Meeko/OpenBabel 并勾选真实 docking，生成每个候选的 PDBQT。" if not ligand_ready_count and not docking_completed else "",
            False,
            "ligand_preparation",
        ),
        item(
            "docking_backend",
            "对接后端",
            "ready" if backend_ready else (plan_status or "missing"),
            f"backend={selected_backend or '-'}; plan_status={plan_status or '-'}",
            "安装 Vina/GNINA，或导入已有真实 docking CSV。" if not backend_ready else "",
            False,
            "docking_backend",
        ),
        item(
            "docking_scores",
            "真实分数文件",
            "ready" if scores_ready else "missing",
            score_file or "no score csv",
            "运行 docking 或在“外部真实分数CSV”中填入文件后重跑 Stage 4。" if not scores_ready else "",
            False,
            "docking_scores",
        ),
        item(
            "interpretation",
            "结果解释",
            "ready" if scores_ready else "pending",
            f"docking_results={len(results)}",
            "有真实分数后再查看 score band、pose 和日志。" if not scores_ready else "",
            False,
            "interpretation",
        ),
    ]

def stage4_payload(name: str, project_dir: Path, round_no: int, command: str = "") -> Dict[str, Any]:
    stage4_dir = project_dir / "stage4"
    descriptor_path = stage4_dir / f"round_{round_no}_real_descriptors.csv"
    similarity_path = stage4_dir / f"round_{round_no}_similarity_to_controls.csv"
    diversity_path = stage4_dir / f"round_{round_no}_diverse_selection.csv"
    decoy_path = stage4_dir / f"round_{round_no}_decoys.csv"
    benchmark_path = stage4_dir / f"round_{round_no}_benchmark_panel.csv"
    receptor_path = stage4_dir / f"round_{round_no}_receptor_package.json"
    docking_inputs_path = stage4_dir / f"round_{round_no}_docking_inputs.csv"
    docking_plan_path = stage4_dir / f"round_{round_no}_docking_plan.json"
    validation_path = stage4_dir / f"round_{round_no}_validation_metrics.json"
    assets_path = stage4_dir / f"round_{round_no}_stage4_assets.json"
    report_path = project_dir / "reports" / f"stage4_round_{round_no}_report.md"
    assets = read_json_if_exists(assets_path)
    docking_results = stage4_docking_results(project_dir, round_no)
    payload = {
        "stage": 4,
        "project": name,
        "round": round_no,
        "has_assets": assets_path.exists(),
        "has_descriptors": descriptor_path.exists(),
        "has_similarity": similarity_path.exists(),
        "has_diverse_selection": diversity_path.exists(),
        "has_decoys": decoy_path.exists(),
        "has_benchmark_panel": benchmark_path.exists(),
        "has_receptor_package": receptor_path.exists(),
        "has_docking_plan": docking_plan_path.exists(),
        "has_validation_metrics": validation_path.exists(),
        "has_report": report_path.exists(),
        "assets": assets,
        "descriptors": csv_to_dicts(descriptor_path),
        "similarity": csv_to_dicts(similarity_path),
        "diverse_selection": csv_to_dicts(diversity_path),
        "decoys": csv_to_dicts(decoy_path),
        "benchmark_panel": csv_to_dicts(benchmark_path),
        "docking_inputs": csv_to_dicts(docking_inputs_path),
        "receptor_package": read_json_if_exists(receptor_path),
        "docking_plan": read_json_if_exists(docking_plan_path),
        "validation_metrics": read_json_if_exists(validation_path),
        "docking_results": docking_results,
        "molecule_images": stage4_image_records(name, project_dir, round_no, assets),
        "report": read_text_if_exists(report_path),
        "files": {
            "candidates": str(project_dir / "candidates" / f"round_{round_no}_candidates.csv"),
            "assets": str(assets_path),
            "real_descriptors": str(descriptor_path),
            "similarity_to_controls": str(similarity_path),
            "diverse_selection": str(diversity_path),
            "decoys": str(decoy_path),
            "benchmark_panel": str(benchmark_path),
            "receptor_package": str(receptor_path),
            "docking_inputs": str(docking_inputs_path),
            "docking_plan": str(docking_plan_path),
            "validation_metrics": str(validation_path),
            "report": str(report_path),
        },
    }
    payload = stage4_apply_recovered_docking(payload, project_dir, round_no)
    if command:
        payload["command"] = command
    payload["visualization_assets"] = stage4_visualization_assets(name, project_dir, payload)
    payload["stage4_readiness_guide"] = stage4_readiness_guide(payload)
    payload["stage4_preflight"] = stage4_preflight(payload)
    return payload

def stage4_repair_generate_candidates(project_dir: Path, round_no: int, target_id: str = "influenza_a_h1n1_na") -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    candidate_path = project_dir / "candidates" / f"round_{safe_round}_candidates.csv"
    candidate_rows = csv_to_dicts(candidate_path)
    if candidate_rows:
        ensure_stage3_assets_from_candidates(project_dir, safe_round, "stage4_repair_existing_candidates")
        return {
            "repair_id": "stage3_candidates",
            "status": "skipped",
            "reason": "candidate csv already exists",
            "file": str(candidate_path),
            "rows": len(candidate_rows),
        }
    try:
        stage3_screen(
            argparse.Namespace(
                project=str(project_dir),
                round=safe_round,
                source_csv=None,
                source_json=None,
                source_url=[],
                context_url=[],
                use_openai=False,
                openai_model="gpt-5.2",
                api_key=None,
                api_key_file=None,
                api_key_env="OPENAI_API_KEY",
                prompt="",
                n=24,
                top=6,
                timeout=60,
                max_url_bytes=800000,
                max_heavy_atoms=55,
                max_molecular_weight=550.0,
                max_lipinski_violations=1,
                allow_risk=False,
                external_scores=None,
                no_score=False,
            )
        )
    except SystemExit as exc:
        return {
            "repair_id": "stage3_candidates",
            "status": "failed",
            "reason": str(exc) or "stage3_screen failed",
            "file": str(candidate_path),
            "rows": 0,
        }
    candidate_rows = csv_to_dicts(candidate_path)
    return {
        "repair_id": "stage3_candidates",
        "status": "completed" if candidate_rows else "failed",
        "reason": "" if candidate_rows else "Stage 3 ran but no candidates passed filters.",
        "file": str(candidate_path),
        "rows": len(candidate_rows),
        "target_id": target_id,
    }

def stage4_repair_project(name: str, project_dir: Path, round_no: int, persist: bool = True) -> Dict[str, Any]:
    payload = stage4_payload(name, project_dir, round_no)
    preflight = payload.get("stage4_preflight", []) if isinstance(payload.get("stage4_preflight"), list) else []
    repairs: List[Dict[str, Any]] = []
    candidate_gate = next((item for item in preflight if item.get("step_id") == "candidate_input"), {})
    if str(candidate_gate.get("status", "")).lower() in {"missing", "empty"}:
        target_id = ""
        assets = payload.get("assets", {}) if isinstance(payload.get("assets"), dict) else {}
        rec = payload.get("receptor_package", {}) if isinstance(payload.get("receptor_package"), dict) else {}
        target_id = str(assets.get("target_id") or rec.get("target_id") or "")
        repairs.append(stage4_repair_generate_candidates(project_dir, round_no, target_id or "influenza_a_h1n1_na"))
        payload = stage4_payload(name, project_dir, round_no)

    recovered = payload.get("recovered_docking_evidence", {})
    metadata_repaired = bool(isinstance(recovered, dict) and recovered.get("status") == "recovered")
    if metadata_repaired:
        repairs.append(
            {
                "repair_id": "recovered_docking_metadata",
                "status": "completed",
                "reason": "Recovered receptor, docking plan, and validation state from Stage 4.5 artifacts.",
                "rows": recovered.get("results", 0),
            }
        )
    repaired = metadata_repaired or any(item.get("status") == "completed" for item in repairs)
    if persist and repaired:
        stage4_dir = project_dir / "stage4"
        rec_path = stage4_dir / f"round_{round_no}_receptor_package.json"
        plan_path = stage4_dir / f"round_{round_no}_docking_plan.json"
        validation_path = stage4_dir / f"round_{round_no}_validation_metrics.json"
        assets_path = stage4_dir / f"round_{round_no}_stage4_assets.json"
        assets = read_json_if_exists(assets_path)
        if assets:
            assets["receptor_package"] = payload.get("receptor_package", {})
            assets["docking_plan"] = payload.get("docking_plan", {})
            assets["validation_metrics"] = payload.get("validation_metrics", {})
            assets["docking_results_count"] = len(payload.get("docking_results", []) or [])
            assets["recovered_docking_evidence"] = recovered
            write_json(assets_path, assets)
        write_json(rec_path, payload.get("receptor_package", {}))
        write_json(plan_path, payload.get("docking_plan", {}))
        write_json(validation_path, payload.get("validation_metrics", {}))
        payload = stage4_payload(name, project_dir, round_no)
    payload["command"] = "stage4-repair"
    payload["repaired"] = repaired
    payload["repairs"] = repairs
    payload["next_actions"] = [
        row.get("next_action", "")
        for row in payload.get("stage4_preflight", [])
        if isinstance(row, dict) and row.get("status") not in {"ready", "ok", "pass"} and row.get("next_action")
    ][:8]
    return payload

def stage4_project_doctor(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    payload = stage4_payload(name, project_dir, round_no)
    stage45 = stage45_payload(name, project_dir, round_no)
    stage46 = stage46_payload(name, project_dir, round_no)
    preflight = payload.get("stage4_preflight", []) if isinstance(payload.get("stage4_preflight"), list) else []
    issues: List[Dict[str, Any]] = []
    for row in preflight:
        if not isinstance(row, dict):
            continue
        status = str(row.get("status", "")).lower()
        if status not in {"ready", "ok", "pass"}:
            issues.append(
                {
                    "area": "stage4_preflight",
                    "severity": "warn" if status in {"pending", "skipped", "empty"} else "fail",
                    "step_id": row.get("step_id", ""),
                    "status": row.get("status", ""),
                    "evidence": row.get("evidence", ""),
                    "next_action": row.get("next_action", ""),
                    "auto_repair": bool(row.get("auto_repair")),
                    "repair_hint": row.get("repair_hint", ""),
                    "repair_endpoint": f"/api/projects/{name}/stage4/repair" if row.get("auto_repair") else "",
                }
            )
    recovered = payload.get("recovered_docking_evidence", {})
    plan = payload.get("docking_plan", {}) if isinstance(payload.get("docking_plan"), dict) else {}
    if isinstance(recovered, dict) and recovered.get("status") == "recovered":
        issues.append(
            {
                "area": "metadata_consistency",
                "severity": "info",
                "step_id": "recovered_docking",
                "status": "recovered",
                "evidence": f"Stage 4 payload recovered {recovered.get('results', 0)} docking rows from Stage 4.5.",
                "next_action": "Use Stage 4 repair to persist recovered receptor, box, and docking status into Stage 4 JSON.",
            }
        )
    return {
        "project": name,
        "round": round_no,
        "status": "needs_repair" if any(item.get("step_id") == "recovered_docking" for item in issues) else ("has_issues" if issues else "ok"),
        "issues": issues,
        "stage4_summary": {
            "receptor": payload.get("stage4_readiness_guide", {}).get("receptor", ""),
            "box": payload.get("stage4_readiness_guide", {}).get("box", ""),
            "docking": payload.get("stage4_readiness_guide", {}).get("docking", ""),
            "docking_results": len(payload.get("docking_results", []) or []),
            "plan_status": plan.get("status", ""),
        },
        "stage45_summary": {
            "has_validation": bool(stage45.get("has_validation")),
            "scores": len(stage45.get("scores", []) or []),
            "docking": (stage45.get("validation", {}) or {}).get("docking", {}) if isinstance(stage45.get("validation"), dict) else {},
        },
        "stage46_summary": {
            "has_benchmark": bool(stage46.get("has_benchmark")),
            "ranking": len(stage46.get("ranking", []) or []),
        },
        "recovered_docking_evidence": recovered,
        "boundary": "Project doctor checks computational assets and metadata consistency only; it does not validate drug efficacy.",
    }

def stage45_payload(name: str, project_dir: Path, round_no: int, command: str = "") -> Dict[str, Any]:
    stage45_dir = project_dir / "stage4_5"
    panel_path = stage45_dir / f"round_{round_no}_control_panel.csv"
    docking_inputs_path = stage45_dir / f"round_{round_no}_control_docking_inputs.csv"
    docking_plan_path = stage45_dir / f"round_{round_no}_control_docking_plan.json"
    scores_path = stage45_dir / f"round_{round_no}_control_docking_scores.csv"
    validation_path = stage45_dir / f"round_{round_no}_control_validation.json"
    report_path = project_dir / "reports" / f"stage4_5_round_{round_no}_control_validation.md"
    raw_scores = csv_to_dicts(scores_path)
    validation = read_json_if_exists(validation_path)
    invalid_backend_reason = stage4_invalid_backend_reason(project_dir, round_no)
    scores = raw_scores
    score_filter = {
        "status": "accepted",
        "raw_scores": len(raw_scores),
        "shown_scores": len(raw_scores),
        "reason": "",
    }
    if invalid_backend_reason and raw_scores:
        scores = []
        score_filter = {
            "status": "ignored_invalid_backend",
            "raw_scores": len(raw_scores),
            "shown_scores": 0,
            "reason": f"Current Stage 4 docking backend is invalid ({invalid_backend_reason}); old Stage 4.5 scores are not shown as usable docking evidence.",
        }
        if isinstance(validation, dict):
            validation = dict(validation)
            docking = dict(validation.get("docking", {}) if isinstance(validation.get("docking"), dict) else {})
            docking["status"] = "not_available"
            docking["backend_status"] = invalid_backend_reason
            docking["ignored_scores"] = len(raw_scores)
            validation["docking"] = docking
    payload = {
        "stage": 4.5,
        "project": name,
        "round": round_no,
        "has_panel": panel_path.exists(),
        "has_docking_inputs": docking_inputs_path.exists(),
        "has_docking_plan": docking_plan_path.exists(),
        "has_scores": scores_path.exists(),
        "has_validation": validation_path.exists(),
        "has_report": report_path.exists(),
        "panel": csv_to_dicts(panel_path),
        "docking_inputs": csv_to_dicts(docking_inputs_path),
        "docking_plan": read_json_if_exists(docking_plan_path),
        "scores": scores,
        "raw_score_count": len(raw_scores),
        "score_filter": score_filter,
        "validation": validation,
        "report": read_text_if_exists(report_path),
        "files": {
            "control_panel": str(panel_path),
            "docking_inputs": str(docking_inputs_path),
            "docking_plan": str(docking_plan_path),
            "docking_scores": str(scores_path),
            "validation": str(validation_path),
            "report": str(report_path),
        },
        "boundary": [
            "Stage 4.5 是计算对照校准，不是实验药效证明。",
            "控药、候选和 decoy 在同一 docking 设置下比较，用于判断分数是否有解释力。",
        ],
    }
    if command:
        payload["command"] = command
    return payload

def stage46_payload(name: str, project_dir: Path, round_no: int, command: str = "") -> Dict[str, Any]:
    stage46_dir = project_dir / "stage4_6"
    benchmark_path = stage46_dir / f"round_{round_no}_retrospective_benchmark.json"
    ranking_path = stage46_dir / f"round_{round_no}_retrospective_ranking.csv"
    report_path = project_dir / "reports" / f"stage4_6_round_{round_no}_retrospective_benchmark.md"
    source_scores = project_dir / "stage4_5" / f"round_{round_no}_control_docking_scores.csv"
    raw_ranking = csv_to_dicts(ranking_path)
    benchmark = read_json_if_exists(benchmark_path)
    invalid_backend_reason = stage4_invalid_backend_reason(project_dir, round_no)
    ranking = raw_ranking
    benchmark_filter = {
        "status": "accepted",
        "raw_ranking": len(raw_ranking),
        "shown_ranking": len(raw_ranking),
        "reason": "",
    }
    if invalid_backend_reason and (raw_ranking or benchmark):
        ranking = []
        benchmark = {}
        benchmark_filter = {
            "status": "ignored_invalid_backend",
            "raw_ranking": len(raw_ranking),
            "shown_ranking": 0,
            "reason": f"Current Stage 4 docking backend is invalid ({invalid_backend_reason}); old Stage 4.6 benchmark is not shown as usable retrospective evidence.",
        }
    payload = {
        "stage": 4.6,
        "project": name,
        "round": round_no,
        "has_benchmark": benchmark_path.exists(),
        "has_ranking": ranking_path.exists(),
        "has_report": report_path.exists(),
        "has_source_scores": source_scores.exists(),
        "benchmark": benchmark,
        "ranking": ranking,
        "raw_ranking_count": len(raw_ranking),
        "benchmark_filter": benchmark_filter,
        "report": read_text_if_exists(report_path),
        "files": {
            "benchmark": str(benchmark_path),
            "ranking": str(ranking_path),
            "source_scores": str(source_scores),
            "report": str(report_path),
        },
        "boundary": [
            "Stage 4.6 是回顾性计算 benchmark，不是湿实验验证。",
            "AUC 与 Top-K 只说明当前 docking 设置能否把已知控药排在 decoy 前面。",
            "候选分子排名只用于优先级排序，不代表真实药效、毒性、安全性或临床价值。",
        ],
    }
    if command:
        payload["command"] = command
    return payload

def stage6_payload(name: str, project_dir: Path, round_no: int, command: str = "") -> Dict[str, Any]:
    stage6_dir = project_dir / "stage6"
    assets_path = stage6_dir / f"round_{round_no}_validation_assets.json"
    gates_path = stage6_dir / f"round_{round_no}_quality_gates.csv"
    triage_path = stage6_dir / f"round_{round_no}_hit_triage.csv"
    queue_path = stage6_dir / f"round_{round_no}_assay_queue.csv"
    risk_path = stage6_dir / f"round_{round_no}_risk_register.csv"
    runbook_path = stage6_dir / f"round_{round_no}_validation_runbook.md"
    report_path = project_dir / "reports" / f"stage6_round_{round_no}_validation_report.md"
    assets = read_json_if_exists(assets_path)
    payload = {
        "stage": 6,
        "project": name,
        "round": round_no,
        "has_assets": assets_path.exists(),
        "has_quality_gates": gates_path.exists(),
        "has_hit_triage": triage_path.exists(),
        "has_assay_queue": queue_path.exists(),
        "has_risk_register": risk_path.exists(),
        "has_runbook": runbook_path.exists(),
        "has_report": report_path.exists(),
        "assets": assets,
        "summary": assets.get("quality_gate_summary", {}) if isinstance(assets, dict) else {},
        "quality_gates": csv_to_dicts(gates_path),
        "hit_triage": csv_to_dicts(triage_path),
        "assay_queue": csv_to_dicts(queue_path),
        "risk_register": csv_to_dicts(risk_path),
        "runbook": read_text_if_exists(runbook_path),
        "report": read_text_if_exists(report_path),
        "files": {
            "assets": str(assets_path),
            "quality_gates": str(gates_path),
            "hit_triage": str(triage_path),
            "assay_queue": str(queue_path),
            "risk_register": str(risk_path),
            "validation_runbook": str(runbook_path),
            "report": str(report_path),
        },
        "boundary": [
            "Stage 6 是验证运营和质量门管理，不是药效证明。",
            "候选分子只能被称为 computational hit 或优先级线索。",
            "真实 docking、Pose QC 和湿实验结果缺失时，界面和报告都不能声称活性、疗效、安全性或临床价值。",
        ],
    }
    if command:
        payload["command"] = command
    return payload

def stage7_payload(name: str, project_dir: Path, round_no: int, command: str = "") -> Dict[str, Any]:
    stage7_dir = project_dir / "stage7"
    manifest_path = stage7_dir / f"round_{round_no}_delivery_manifest.json"
    summary_path = stage7_dir / f"round_{round_no}_executive_summary.md"
    repro_path = stage7_dir / f"round_{round_no}_reproducibility.md"
    checklist_path = stage7_dir / f"round_{round_no}_investor_demo_checklist.csv"
    stage8_spec_path = stage7_dir / "stage8_frontend_product_spec.md"
    report_path = project_dir / "reports" / f"stage7_round_{round_no}_delivery_report.md"
    manifest = read_json_if_exists(manifest_path)
    deliverables = manifest.get("deliverables", {}) if isinstance(manifest, dict) else {}
    boundary = manifest.get("claims_boundary", []) if isinstance(manifest, dict) else []
    if not boundary:
        boundary = [
            "Stage 7 只打包产品交付材料和复现说明，不新增科学验证。",
            "所有输出仍属于计算筛选与验证计划，不构成药效、毒性、安全性或临床价值证明。",
            "已知药物在本项目中仅作为控药和 benchmark。",
        ]
    payload = {
        "stage": 7,
        "project": name,
        "round": round_no,
        "has_manifest": manifest_path.exists(),
        "has_executive_summary": summary_path.exists(),
        "has_reproducibility": repro_path.exists(),
        "has_checklist": checklist_path.exists(),
        "has_stage8_spec": stage8_spec_path.exists(),
        "has_report": report_path.exists(),
        "manifest": manifest,
        "deliverables": deliverables,
        "executive_summary": read_text_if_exists(summary_path),
        "reproducibility": read_text_if_exists(repro_path),
        "checklist": csv_to_dicts(checklist_path),
        "stage8_spec": read_text_if_exists(stage8_spec_path),
        "report": read_text_if_exists(report_path),
        "files": {
            "manifest": str(manifest_path),
            "executive_summary": str(summary_path),
            "reproducibility": str(repro_path),
            "demo_checklist": str(checklist_path),
            "stage8_frontend_spec": str(stage8_spec_path),
            "report": str(report_path),
        },
        "boundary": boundary,
    }
    if command:
        payload["command"] = command
    return payload

def stage8_status_from_readiness(value: Any) -> str:
    raw = str(value or "").strip().lower()
    if not raw:
        return "missing"
    if "missing" in raw or "failed" in raw or "blocked" in raw or "not_available" in raw:
        return "missing"
    if "skipped" in raw or "pending" in raw or "warn" in raw or "needs" in raw:
        return "warn"
    if "ready" in raw or "completed" in raw or "available" in raw or "pass" in raw:
        return "ready"
    return "warn"

def stage8_gate_summary(gates: List[Dict[str, Any]]) -> Dict[str, Any]:
    summary = {"pass": 0, "warn": 0, "fail": 0, "total": len(gates), "overall_status": "missing"}
    for row in gates:
        status = str(row.get("status", "")).lower()
        if status in summary:
            summary[status] += 1
    if summary["fail"]:
        summary["overall_status"] = "blocked_by_quality_gates"
    elif summary["warn"]:
        summary["overall_status"] = "ready_with_warnings"
    elif summary["total"]:
        summary["overall_status"] = "ready_for_computational_demo"
    return summary

def stage8_deliverable_rows(deliverables: Dict[str, Any]) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for key in sorted(deliverables.keys()):
        item = deliverables.get(key, {})
        if not isinstance(item, dict):
            item = {}
        rows.append(
            {
                "key": key,
                "exists": bool(item.get("exists")),
                "path": str(item.get("path") or item.get("relative_path") or ""),
                "relative_path": str(item.get("relative_path") or ""),
                "description": str(item.get("description") or ""),
            }
        )
    return rows

def stage8_stage_rail(
    dashboard: Dict[str, Any],
    stage6: Dict[str, Any],
    stage7: Dict[str, Any],
    has_stage5: bool,
) -> List[Dict[str, Any]]:
    readiness = dashboard.get("readiness", {}) if isinstance(dashboard.get("readiness"), dict) else {}
    metrics = dashboard.get("metrics", {}) if isinstance(dashboard.get("metrics"), dict) else {}
    target = dashboard.get("target", {}) if isinstance(dashboard.get("target"), dict) else {}
    stage6_summary = stage6.get("summary", {}) if isinstance(stage6.get("summary"), dict) else {}
    stage7_manifest = stage7.get("manifest", {}) if isinstance(stage7.get("manifest"), dict) else {}

    stage4_statuses = [
        stage8_status_from_readiness(readiness.get("rdkit_validation")),
        stage8_status_from_readiness(readiness.get("control_panel")),
        stage8_status_from_readiness(readiness.get("decoy_panel")),
    ]
    docking_status = stage8_status_from_readiness(readiness.get("docking"))
    if "missing" in stage4_statuses:
        stage4_status = "missing"
    elif docking_status == "ready":
        stage4_status = "ready"
    else:
        stage4_status = "warn"

    if stage6.get("has_assets"):
        stage6_status = "warn" if int(stage6_summary.get("fail", 0) or 0) else "ready"
    else:
        stage6_status = "missing"

    if stage7.get("has_manifest"):
        delivery_status = str(stage7_manifest.get("delivery_status", ""))
        stage7_status = "ready" if delivery_status.startswith("ready") else "warn"
    else:
        stage7_status = "missing"

    return [
        {
            "stage": 1,
            "label": "靶点需求",
            "status": "ready" if target.get("id") or target.get("name") else "missing",
            "evidence": str(target.get("name") or target.get("id") or "missing target brief"),
            "route": "brief",
        },
        {
            "stage": 2,
            "label": "靶点证据",
            "status": stage8_status_from_readiness(readiness.get("target_evidence")),
            "evidence": f"target_evidence={readiness.get('target_evidence', 'missing')}",
            "route": "stage2",
        },
        {
            "stage": 3,
            "label": "候选输入",
            "status": stage8_status_from_readiness(readiness.get("candidate_intake")),
            "evidence": f"raw={metrics.get('raw_candidates', 0)}; ranked={metrics.get('ranked_candidates', 0)}",
            "route": "results",
        },
        {
            "stage": 4,
            "label": "真实库校验",
            "status": stage4_status,
            "evidence": (
                f"rdkit={readiness.get('rdkit_validation', 'missing')}; "
                f"controls={readiness.get('control_panel', 'missing')}; "
                f"decoys={readiness.get('decoy_panel', 'missing')}; "
                f"docking={readiness.get('docking', 'missing')}"
            ),
            "route": "stage4",
        },
        {
            "stage": 5,
            "label": "产品看板",
            "status": "ready" if has_stage5 else "warn",
            "evidence": "dashboard_data.json exists" if has_stage5 else "live aggregation only; dashboard file not generated",
            "route": "dashboard",
        },
        {
            "stage": 6,
            "label": "验证运营",
            "status": stage6_status,
            "evidence": f"quality_gates={stage6_summary.get('total', 0)}; status={stage6_summary.get('overall_status', 'missing')}",
            "route": "stage6",
        },
        {
            "stage": 7,
            "label": "交付包",
            "status": stage7_status,
            "evidence": str(stage7_manifest.get("delivery_status") or "missing delivery package"),
            "route": "stage7",
        },
    ]

def stage8_next_actions(stage6: Dict[str, Any], stage7: Dict[str, Any], gates: List[Dict[str, Any]]) -> List[str]:
    actions: List[str] = []
    assets = stage6.get("assets", {}) if isinstance(stage6.get("assets"), dict) else {}
    raw_actions = assets.get("next_actions", []) if isinstance(assets, dict) else []
    if isinstance(raw_actions, list):
        actions.extend(str(item) for item in raw_actions if str(item).strip())
    if not stage6.get("has_assets"):
        actions.append("Run Stage 6 validation operations to create quality gates, hit triage, assay queue, and risk register.")
    if not stage7.get("has_manifest"):
        actions.append("Generate Stage 7 delivery package after Stage 6 is available.")
    for gate in gates:
        if str(gate.get("status", "")).lower() in {"warn", "fail"}:
            step = str(gate.get("required_next_step", "")).strip()
            if step:
                actions.append(step)
    if not actions:
        actions.append("Review the stage rail and keep the computational-screening claim boundary visible before external sharing.")
    deduped: List[str] = []
    seen = set()
    for action in actions:
        if action not in seen:
            deduped.append(action)
            seen.add(action)
    return deduped[:8]

def stage8_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    stage5_data_path = project_dir / "stage5" / "dashboard_data.json"
    stage5_html_path = project_dir / "stage5" / "index.html"
    stage5_report_path = project_dir / "reports" / "stage5_dashboard_report.md"
    existing_dashboard = read_json_if_exists(stage5_data_path)
    title = str(existing_dashboard.get("title", "") or f"{name} Stage 8 Command Center")
    dashboard = build_stage5_dashboard_data(project_dir, round_no, title)
    stage4 = stage4_payload(name, project_dir, round_no)
    stage4_guide = stage4.get("stage4_readiness_guide", {}) if isinstance(stage4.get("stage4_readiness_guide"), dict) else {}
    readiness = dashboard.get("readiness", {}) if isinstance(dashboard.get("readiness"), dict) else {}
    if stage4_guide.get("docking") == "completed":
        readiness["docking"] = "completed"
        dashboard["readiness"] = readiness
    stage6 = stage6_payload(name, project_dir, round_no)
    stage7 = stage7_payload(name, project_dir, round_no)
    gates = stage6.get("quality_gates", []) if isinstance(stage6.get("quality_gates"), list) else []
    stage6_summary = stage6.get("summary", {}) if isinstance(stage6.get("summary"), dict) else {}
    gate_summary = stage6_summary or stage8_gate_summary(gates)
    deliverables = stage8_deliverable_rows(stage7.get("deliverables", {}) if isinstance(stage7.get("deliverables"), dict) else {})
    metrics = dashboard.get("metrics", {}) if isinstance(dashboard.get("metrics"), dict) else {}
    assets = stage6.get("assets", {}) if isinstance(stage6.get("assets"), dict) else {}
    stage7_manifest = stage7.get("manifest", {}) if isinstance(stage7.get("manifest"), dict) else {}

    return {
        "stage": 8,
        "project": name,
        "round": round_no,
        "generated_at": datetime.now().isoformat(),
        "has_stage5": stage5_data_path.exists(),
        "has_stage5_html": stage5_html_path.exists(),
        "has_stage6": bool(stage6.get("has_assets")),
        "has_stage7": bool(stage7.get("has_manifest")),
        "target": dashboard.get("target", {}),
        "readiness": dashboard.get("readiness", {}),
        "stage4_consistency": stage4.get("recovered_docking_evidence", {}),
        "metrics": metrics,
        "stage_rail": stage8_stage_rail(dashboard, stage6, stage7, stage5_data_path.exists()),
        "funnel": [
            {"label": "raw_candidates", "value": metrics.get("raw_candidates", 0)},
            {"label": "filtered_candidates", "value": metrics.get("filtered_candidates", 0)},
            {"label": "scored_candidates", "value": metrics.get("scored_candidates", 0)},
            {"label": "ranked_candidates", "value": metrics.get("ranked_candidates", 0)},
            {"label": "advanced", "value": metrics.get("advanced", 0)},
            {"label": "validation_queue", "value": assets.get("queue_count", len(stage6.get("assay_queue", []) or []))},
        ],
        "quality_gate_summary": gate_summary,
        "quality_gates": gates[:12],
        "hit_triage": (stage6.get("hit_triage", []) or [])[:5],
        "assay_queue": (stage6.get("assay_queue", []) or [])[:8],
        "deliverables": deliverables,
        "delivery_status": stage7_manifest.get("delivery_status", ""),
        "stage6_status": assets.get("overall_status", ""),
        "next_actions": stage8_next_actions(stage6, stage7, gates),
        "files": {
            "stage5_dashboard_data": str(stage5_data_path),
            "stage5_dashboard_html": str(stage5_html_path),
            "stage5_report": str(stage5_report_path),
            "stage6_assets": str(project_dir / "stage6" / f"round_{round_no}_validation_assets.json"),
            "stage6_quality_gates": str(project_dir / "stage6" / f"round_{round_no}_quality_gates.csv"),
            "stage7_manifest": str(project_dir / "stage7" / f"round_{round_no}_delivery_manifest.json"),
            "stage8_frontend_spec": str(project_dir / "stage7" / "stage8_frontend_product_spec.md"),
        },
        "download_links": stage8_download_links(name, project_dir, round_no),
        "boundary": [
            "Computational screening and validation-planning command center only.",
            "No efficacy, potency, toxicity, clinical benefit, dosing, or safety claim is created by this workflow.",
            "Docking, RDKit, controls, decoys, and benchmark panels are decision-support evidence, not wet-lab validation.",
            "Known drugs are shown as controls and references, not as newly discovered project outputs.",
        ],
    }

def stage8_review_mode_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    command_center = stage8_payload(name, project_dir, safe_round)
    guide = stage8_demo_guide_payload(name, project_dir, safe_round)
    preflight = stage8_preflight_payload(name, project_dir, safe_round)
    target = command_center.get("target", {}) if isinstance(command_center.get("target"), dict) else {}
    metrics = command_center.get("metrics", {}) if isinstance(command_center.get("metrics"), dict) else {}
    gate_summary = command_center.get("quality_gate_summary", {}) if isinstance(command_center.get("quality_gate_summary"), dict) else {}
    stage4_consistency = command_center.get("stage4_consistency", {}) if isinstance(command_center.get("stage4_consistency"), dict) else {}
    completion = guide.get("completion", {}) if isinstance(guide.get("completion"), dict) else {}
    ready_steps = int(completion.get("ready", 0) or 0)
    warn_steps = int(completion.get("warn", 0) or 0)
    missing_steps = int(completion.get("missing", 0) or 0)
    if missing_steps:
        review_status = "needs_repair"
    elif warn_steps:
        review_status = "ready_with_warnings"
    else:
        review_status = "ready_for_review"

    target_name = str(target.get("name") or target.get("id") or "selected target")
    disease = str(target.get("disease_context") or "selected disease context")
    primary_pdb = str(target.get("primary_pdb") or "-")
    reference_ligand = str(target.get("reference_ligand") or "-")
    readiness = command_center.get("readiness", {}) if isinstance(command_center.get("readiness"), dict) else {}
    docking_ready = stage8_status_from_readiness(readiness.get("docking")) == "ready"
    docking_note = "真实 docking / 外部分数已接入" if docking_ready else "已生成 docking-ready 资产或验证计划"

    storyline = [
        {
            "section_id": "project_goal",
            "title": "项目目标",
            "status": "ready" if target.get("id") or target.get("name") else "missing",
            "headline": f"围绕 {disease} 的 {target_name} 建立候选分子计算筛选演示。",
            "evidence": f"Project={name}; round={safe_round}; ranked={metrics.get('ranked_candidates', 0)}.",
            "reviewer_message": "本项目展示从靶点需求到候选排序、验证运营和证据交付的产品化链路。",
        },
        {
            "section_id": "target_evidence",
            "title": "靶点证据",
            "status": stage8_status_from_readiness(readiness.get("target_evidence")),
            "headline": f"靶点使用公开结构和控药信息约束，PDB={primary_pdb}，参考配体={reference_ligand}。",
            "evidence": f"evidence_score={target.get('evidence_score', '-')}; readiness={target.get('evidence_readiness', '-')}.",
            "reviewer_message": "评审重点是靶点、口袋、控药和公开证据是否足以支撑虚拟筛选，不把它解释为实验验证。",
        },
        {
            "section_id": "candidate_screening",
            "title": "候选筛选",
            "status": stage8_status_from_readiness(readiness.get("candidate_intake")),
            "headline": f"候选漏斗已完成：原始 {metrics.get('raw_candidates', 0)}，排序 {metrics.get('ranked_candidates', 0)}，晋级 {metrics.get('advanced', 0)}。",
            "evidence": f"filtered={metrics.get('filtered_candidates', 0)}; scored={metrics.get('scored_candidates', 0)}.",
            "reviewer_message": "候选由生成或导入入口进入同一评分表，后续反馈只改变下一轮种子，不声称产生临床有效药物。",
        },
        {
            "section_id": "computational_validation",
            "title": "计算验证",
            "status": "ready" if stage8_status_from_readiness(readiness.get("rdkit_validation")) == "ready" else "warn",
            "headline": f"Stage 4 聚合 RDKit、控药、Decoy 和 docking 证据；当前 docking 状态：{docking_note}。",
            "evidence": f"quality_gates={gate_summary.get('total', 0)}; recovered_docking={stage4_consistency.get('has_recovered_scores', False)}.",
            "reviewer_message": "真实 docking 分数只在后端实际运行或外部分数导入后展示，Pose 和对照面板用于可信度解释。",
        },
        {
            "section_id": "delivery_evidence",
            "title": "交付证据",
            "status": "ready" if command_center.get("has_stage7") else "warn",
            "headline": "Stage 5/6/7/8 将 Dashboard、质量门、复现说明和证据包合并为可交付材料。",
            "evidence": f"stage6={command_center.get('has_stage6')}; stage7={command_center.get('has_stage7')}; preflight={preflight.get('status')}.",
            "reviewer_message": "评审可以直接查看 Dashboard、闭环报告和证据包 ZIP 来复现项目材料来源。",
        },
    ]

    download_links = dict(command_center.get("download_links", {}) if isinstance(command_center.get("download_links"), dict) else {})
    extra_downloads = {
        "stage8_report": project_dir / "reports" / f"stage8_round_{safe_round}_closed_loop_report.md",
        "evidence_pack_zip": project_dir / "exports" / f"round_{safe_round}_evidence_pack.zip",
        "evidence_manifest": project_dir / "exports" / f"round_{safe_round}_evidence_manifest.json",
    }
    for key, path in extra_downloads.items():
        if path.exists() and path.is_file():
            url = project_artifact_url(name, project_dir, path)
            if url:
                download_links[key] = url

    primary_actions = {
        "refresh_review": f"GET /api/projects/{name}/stage8/review-mode?round={safe_round}",
        "stage8_report": f"POST /api/projects/{name}/stage8/report",
        "evidence_pack": f"POST /api/projects/{name}/stage8/evidence-pack",
        "demo_doctor": f"GET /api/projects/{name}/demo-doctor?round={safe_round}",
    }
    if review_status == "needs_repair":
        primary_actions["repair"] = f"POST /api/projects/{name}/stage8/repair"

    return {
        "stage": 8,
        "command": "stage8-review-mode",
        "product_mode": "review_mode",
        "project": name,
        "round": safe_round,
        "generated_at": datetime.now().isoformat(),
        "review_status": review_status,
        "positioning": (
            "Computational screening product demo for target-driven candidate generation, "
            "scoring, validation planning, and evidence handoff."
        ),
        "target": {
            "id": target.get("id") or target.get("name") or "",
            "name": target_name,
            "disease_context": disease,
            "primary_pdb": primary_pdb,
            "reference_ligand": reference_ligand,
            "positive_controls": target.get("positive_controls", ""),
            "evidence_score": target.get("evidence_score", ""),
            "evidence_readiness": target.get("evidence_readiness", ""),
        },
        "storyline": storyline,
        "evidence_strength": {
            "ready_steps": ready_steps,
            "warn_steps": warn_steps,
            "missing_steps": missing_steps,
            "completion_percent": completion.get("percent", 0),
            "quality_gate_status": gate_summary.get("overall_status", "missing"),
            "downloads": len(download_links),
        },
        "talk_track": [
            f"本项目先确定 {disease} / {target_name}，再用公开结构、参考配体和已知控药约束候选设计。",
            f"候选分子进入统一入口后完成过滤、代理评分、排序和反馈种子生成，当前排序候选数为 {metrics.get('ranked_candidates', 0)}。",
            "Stage 4 只解释计算证据：RDKit、控药/Decoy、docking 资产或真实 docking 分数；实验活性需要后续湿实验验证。",
            "Stage 5 到 Stage 8 把结果收束成 Dashboard、质量门、交付清单、闭环报告和证据包，便于评审复核。",
        ],
        "primary_actions": primary_actions,
        "download_links": download_links,
        "next_primary_action": guide.get("next_primary_action", {}),
        "command_center": command_center,
        "claim_boundary": [
            "Computational screening and validation planning only.",
            "No efficacy, potency, toxicity, dosing, safety, clinical benefit, or therapeutic claim is created.",
            "Known drugs are controls and references, not newly discovered outputs.",
            "Docking and pose evidence are decision-support signals and require experimental validation.",
        ],
    }

def demo_doctor_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    stage8 = stage8_payload(name, project_dir, round_no)
    project_doctor = stage4_project_doctor(name, project_dir, round_no)
    capabilities = stage4_capabilities()

    stage_rail = stage8.get("stage_rail", []) if isinstance(stage8.get("stage_rail"), list) else []
    stage_counts = {"ready": 0, "warn": 0, "missing": 0}
    pipeline_issues: List[Dict[str, Any]] = []
    for row in stage_rail:
        if not isinstance(row, dict):
            continue
        status = str(row.get("status") or "warn").lower()
        if status not in stage_counts:
            status = "warn"
        stage_counts[status] += 1
        if status != "ready":
            pipeline_issues.append(
                {
                    "check": "stage_pipeline",
                    "severity": "warn",
                    "stage": row.get("stage"),
                    "evidence": f"{row.get('label', 'stage')} is {status}: {row.get('evidence', '')}",
                    "next_action": f"Open {row.get('label', 'the stage')} and complete its missing assets.",
                }
            )
    pipeline_status = "ready" if stage_counts["missing"] == 0 and stage_counts["warn"] == 0 else "warn"

    static_path = STATIC / "vendor" / "3Dmol-min.js"
    static_size = static_path.stat().st_size if static_path.exists() else 0
    static_3dmol = {
        "status": "ready" if static_path.exists() and static_size > 100_000 else "missing",
        "path": str(static_path),
        "bytes": static_size,
        "url": "/static/vendor/3Dmol-min.js",
    }

    docking_backend = capabilities.get("docking_backend", {}) if isinstance(capabilities.get("docking_backend"), dict) else {}
    executables = capabilities.get("executables", {}) if isinstance(capabilities.get("executables"), dict) else {}
    available_backends = docking_backend.get("available_backends", []) if isinstance(docking_backend, dict) else []
    found_tools = [
        key for key, value in executables.items()
        if isinstance(value, dict) and value.get("status") == "found"
    ]
    system_tools = {
        "status": "ready" if available_backends or found_tools else "warn",
        "docking_backend": docking_backend,
        "found_tools": found_tools,
        "proxy_only": not bool(available_backends or found_tools),
    }

    stage8_check = {
        "status": "ready" if stage8.get("has_stage6") and stage8.get("has_stage7") else "warn",
        "has_stage5": bool(stage8.get("has_stage5")),
        "has_stage6": bool(stage8.get("has_stage6")),
        "has_stage7": bool(stage8.get("has_stage7")),
        "downloads": len(stage8.get("download_links", {}) or {}),
        "delivery_status": stage8.get("delivery_status", ""),
    }

    checks = {
        "stage_pipeline": {
            "status": pipeline_status,
            "counts": stage_counts,
            "stages": stage_rail,
        },
        "static_3dmol": static_3dmol,
        "system_tools": system_tools,
        "stage8": stage8_check,
    }

    issues = list(pipeline_issues)
    if static_3dmol["status"] != "ready":
        issues.append(
            {
                "check": "static_3dmol",
                "severity": "warn",
                "evidence": "Local 3Dmol.js vendor file is missing or unexpectedly small.",
                "next_action": "Restore webapp/static/vendor/3Dmol-min.js so receptor/pose viewing works offline.",
            }
        )
    if system_tools["status"] != "ready":
        issues.append(
            {
                "check": "system_tools",
                "severity": "warn",
                "evidence": "No local docking backend or recognized computational executable was detected.",
                "next_action": "Install/activate Vina, GNINA, OpenBabel, Meeko, PoseBusters, or import external real-score CSVs.",
            }
        )
    if stage8_check["status"] != "ready":
        issues.append(
            {
                "check": "stage8",
                "severity": "warn",
                "evidence": "Stage 6 validation operations or Stage 7 delivery manifest is missing.",
                "next_action": "Run Stage 8 demo package, or run Stage 6 then Stage 7 individually.",
            }
        )
    for item in project_doctor.get("issues", []) if isinstance(project_doctor.get("issues"), list) else []:
        issues.append(
            {
                "check": "stage4_project_doctor",
                "severity": item.get("severity", "warn"),
                "evidence": item.get("evidence", ""),
                "next_action": item.get("next_action", ""),
            }
        )

    next_actions = list(stage8.get("next_actions", []) or [])
    for issue in issues:
        action = str(issue.get("next_action", "")).strip()
        if action and action not in next_actions:
            next_actions.append(action)
    if not next_actions:
        next_actions.append("Use the Stage 8 command center for demo presentation; keep the computational-only boundary visible.")

    ready = (
        checks["stage_pipeline"]["status"] == "ready"
        and checks["static_3dmol"]["status"] == "ready"
        and checks["stage8"]["status"] == "ready"
        and not project_doctor.get("issues")
    )
    return {
        "stage": 8,
        "project": name,
        "round": round_no,
        "generated_at": datetime.now().isoformat(),
        "overall_status": "ready_for_demo" if ready else "ready_with_warnings",
        "checks": checks,
        "issues": issues[:20],
        "next_actions": next_actions[:10],
        "stage8": stage8,
        "project_doctor": project_doctor,
        "boundary": [
            "Demo Doctor checks demo readiness, static assets, computational tool availability, and project artifact consistency.",
            "It does not validate biological activity, efficacy, toxicity, safety, dosing, or clinical usefulness.",
        ],
    }

def stage8_step(step_id: str, label: str, status: str, message: str = "", files: Optional[Dict[str, str]] = None, details: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    return {
        "step_id": step_id,
        "label": label,
        "status": status,
        "message": message,
        "files": files or {},
        "details": details or {},
    }

def stage8_run_cli_step(
    steps: List[Dict[str, Any]],
    step_id: str,
    label: str,
    fn,
    args: argparse.Namespace,
    files: Optional[Dict[str, str]] = None,
    details: Optional[Dict[str, Any]] = None,
    warning_ok: bool = False,
) -> bool:
    try:
        fn(args)
    except SystemExit as exc:
        message = str(exc) or "CLI command failed"
        if warning_ok:
            steps.append(stage8_step(step_id, label, "completed_with_warnings", message, files, details))
            return False
        steps.append(stage8_step(step_id, label, "failed", message, files, details))
        return False
    except Exception as exc:
        if warning_ok:
            steps.append(stage8_step(step_id, label, "completed_with_warnings", str(exc), files, details))
            return False
        steps.append(stage8_step(step_id, label, "failed", str(exc), files, details))
        return False
    steps.append(stage8_step(step_id, label, "completed", "completed", files, details))
    return True

def stage8_triplet_text(value: Any) -> str:
    if isinstance(value, str):
        values = parse_triplet(value)
    elif isinstance(value, list):
        try:
            values = [float(item) for item in value[:3]]
        except Exception:
            values = []
    else:
        values = []
    if len(values) != 3:
        return ""
    return ",".join(str(round(float(item), 6)).rstrip("0").rstrip(".") for item in values)

def stage8_seed_receptor_from_local_cache(project_dir: Path, pdb_id: str) -> Dict[str, str]:
    if not pdb_id:
        return {}
    existing_pdb, existing_pdbqt = find_stage4_project_receptor(project_dir, pdb_id)
    if existing_pdb or existing_pdbqt:
        return {"status": "already_present", "pdb": existing_pdb, "pdbqt": existing_pdbqt}
    pdb_key = pdb_id.upper()
    receptor_dir = project_dir / "stage4" / "receptors"
    search_roots = [
        PROJECTS_ROOT,
        ROOT / "snapshots",
        ROOT / "deliverables",
        ROOT / "ai_mol_loop" / "demo_project",
    ]
    pdb_candidates: List[Path] = []
    pdbqt_candidates: List[Path] = []
    for root in search_roots:
        if not root.exists():
            continue
        for candidate in root.rglob(f"{pdb_key}.pdb"):
            if project_dir in candidate.parents:
                continue
            pdb_candidates.append(candidate)
        for pattern in [f"{pdb_key}.pdbqt", f"{pdb_key}_protein_only_obabel.pdbqt", f"{pdb_key}_protein_only.pdbqt"]:
            for candidate in root.rglob(pattern):
                if project_dir in candidate.parents:
                    continue
                pdbqt_candidates.append(candidate)
    copied: Dict[str, str] = {"status": "missing_local_cache"}
    receptor_dir.mkdir(parents=True, exist_ok=True)
    if pdb_candidates:
        dst = receptor_dir / f"{pdb_key}.pdb"
        shutil.copy2(pdb_candidates[0], dst)
        copied.update({"status": "copied", "pdb": str(dst), "pdb_source": str(pdb_candidates[0])})
    if pdbqt_candidates:
        suffix_name = pdbqt_candidates[0].name
        dst = receptor_dir / suffix_name
        shutil.copy2(pdbqt_candidates[0], dst)
        copied.update({"status": "copied", "pdbqt": str(dst), "pdbqt_source": str(pdbqt_candidates[0])})
    return copied

def stage8_stage4_real_args(
    project_dir: Path,
    safe_round: int,
    req: Stage8DemoRunnerRequest,
    target_id: str,
    target_pack: Dict[str, Any],
    first_structure: Dict[str, Any],
    backend: str,
) -> argparse.Namespace:
    pocket = target_pack.get("pocket", {}) if isinstance(target_pack.get("pocket"), dict) else {}
    pocket_center = stage8_triplet_text(pocket.get("center"))
    pocket_size = stage8_triplet_text(pocket.get("size"))
    pdb_id = str(target_pack.get("pdb_id") or first_structure.get("pdb_id", "") or "")
    return argparse.Namespace(
        project=str(project_dir),
        round=safe_round,
        target=target_id,
        input_csv=None,
        controls_csv=None,
        receptor_pdb=None,
        pdb_id=pdb_id,
        pocket_center=pocket_center,
        pocket_size=pocket_size,
        pocket_source=str(pocket.get("source", "") or "stage8_demo_runner"),
        fetch_receptor=bool(req.run_docking),
        top=max(1, min(int(req.top), 50)),
        decoys=max(0, min(int(req.decoys), 100)),
        max_conformers=1,
        seed=61453,
        no_sdf=False,
        render_2d=True,
        docking_backend=backend,
        run_docking=bool(req.run_docking),
        docking_timeout=max(1, min(int(req.docking_timeout), 7200)),
        rescore=False,
        rank_top=0,
        feedback_top=0,
        external_scores=None,
    )

def stage8_demo_runner(name: str, project_dir: Path, req: Stage8DemoRunnerRequest) -> Dict[str, Any]:
    ensure_project_dirs(project_dir)
    safe_round = max(1, int(req.round))
    target_id = str(req.target or "influenza_a_h1n1_na").strip()
    steps: List[Dict[str, Any]] = []

    catalog = load_target_catalog(None)
    target, matched_by = select_catalog_target(
        catalog,
        disease=req.disease or "甲流",
        target_id=target_id,
        target_hint=req.target_hint or "",
        source_text=req.source or "",
        pdb_id=extract_pdb_id_from_text(req.source or ""),
    )
    if not target:
        raise HTTPException(404, "Demo Runner 未能解析出可用靶点")
    target_id = str(target.get("id", target_id))
    structures = target.get("representative_structures", [])
    first_structure = structures[0] if isinstance(structures, list) and structures and isinstance(structures[0], dict) else {}
    site = target.get("binding_site", {}) if isinstance(target.get("binding_site"), dict) else {}

    stage8_run_cli_step(
        steps,
        "target_selection",
        "Stage 1 target selection",
        target_select,
        argparse.Namespace(
            project=str(project_dir),
            disease=req.disease or "甲流",
            target=None,
            top=5,
            catalog=None,
        ),
        {
            "target_selection": str(project_dir / "targets" / "target_selection.csv"),
            "target_report": str(project_dir / "reports" / "target_selection.md"),
        },
        {"selected_target_id": target_id, "matched_by": matched_by},
    )

    intake = stage1_persist_catalog_brief(
        project_dir,
        target,
        req.disease or "甲流",
        req.source_kind or "url",
        req.source or "",
        req.target_hint or "",
        matched_by,
        force=True,
    )
    steps.append(stage8_step(
        "target_intake",
        "Stage 1 target intake",
        "completed",
        f"matched_by={matched_by}",
        {
            "target_brief": str(intake["brief_path"]),
            "generator_prompt": str(intake["prompt_path"]),
            "config": str(intake["config_path"]),
        },
        {"target_id": target_id, "target_name": target.get("display_name", "")},
    ))

    target_pack = stage1_target_pack_data(
        project_dir,
        safe_round,
        target_id_override=target_id,
        pdb_id_override=str(first_structure.get("pdb_id", "")),
        reference_ligand_override=str(site.get("reference_ligand", "")),
    )
    steps.append(stage8_step(
        "target_pack",
        "Stage 1 target pack",
        "completed" if target_pack.get("target_id") else "completed_with_warnings",
        "target pack generated",
        target_pack.get("files", {}),
        {"pdb_id": target_pack.get("pdb_id", ""), "reference_ligand": target_pack.get("reference_ligand", "")},
    ))

    stage8_run_cli_step(
        steps,
        "stage2_evidence",
        "Stage 2 evidence matrix",
        evidence_stage2,
        argparse.Namespace(
            project=str(project_dir),
            disease=req.disease or "甲流",
            target=target_id,
            top=5,
            catalog=None,
            sources=None,
            evidence_dir=None,
            refresh=False,
            retmax=5,
            timeout=20,
            offline=True,
        ),
        {
            "matrix": str(project_dir / "evidence" / "stage2_target_sources.csv"),
            "assets": str(project_dir / "evidence" / "stage2_closed_loop_assets.json"),
            "report": str(project_dir / "reports" / "stage2_evidence_report.md"),
        },
    )

    stage8_run_cli_step(
        steps,
        "candidate_generation",
        "Candidate generation",
        generate_candidates,
        argparse.Namespace(project=str(project_dir), round=safe_round, n=max(1, min(int(req.candidates), 200)), source_csv=None),
        {"candidates": str(project_dir / "candidates" / f"round_{safe_round}_candidates.csv")},
    )
    stage8_run_cli_step(
        steps,
        "proxy_scoring",
        "Proxy scoring",
        score_candidates,
        argparse.Namespace(project=str(project_dir), round=safe_round, external_scores=None, real_descriptors=None),
        {"scores": str(project_dir / "scores" / f"round_{safe_round}_scores.csv")},
    )
    stage8_run_cli_step(
        steps,
        "ranking",
        "Candidate ranking",
        rank_candidates,
        argparse.Namespace(project=str(project_dir), round=safe_round, top=max(1, min(int(req.top), 50))),
        {"ranked": str(project_dir / "ranked" / f"round_{safe_round}_ranked.csv")},
    )
    stage8_run_cli_step(
        steps,
        "feedback",
        "Feedback seed package",
        feedback,
        argparse.Namespace(project=str(project_dir), round=safe_round, top=max(1, min(int(req.top), 50))),
        {
            "feedback": str(project_dir / "feedback" / f"round_{safe_round}_feedback.json"),
            "seeds": str(project_dir / "seeds" / f"round_{safe_round}_seeds.csv"),
        },
    )

    stage4_run_docking = bool(req.run_docking)
    backend = req.docking_backend if req.docking_backend in {"auto", "vina", "gnina"} else "auto"
    receptor_cache = {}
    if stage4_run_docking:
        receptor_cache = stage8_seed_receptor_from_local_cache(project_dir, str(target_pack.get("pdb_id") or first_structure.get("pdb_id", "")))
    stage8_run_cli_step(
        steps,
        "stage4_assets",
        "Stage 4 real-library assets",
        stage4_real,
        stage8_stage4_real_args(
            project_dir,
            safe_round,
            req,
            target_id,
            target_pack,
            first_structure,
            backend,
        ),
        {
            "assets": str(project_dir / "stage4" / f"round_{safe_round}_stage4_assets.json"),
            "docking_plan": str(project_dir / "stage4" / f"round_{safe_round}_docking_plan.json"),
            "receptor_package": str(project_dir / "stage4" / f"round_{safe_round}_receptor_package.json"),
        },
        {"receptor_cache": receptor_cache} if receptor_cache else None,
        warning_ok=True,
    )
    stage8_run_cli_step(
        steps,
        "stage45_controls",
        "Stage 4.5 control calibration",
        stage45_validate_controls,
        argparse.Namespace(
            project=str(project_dir),
            round=safe_round,
            target=target_id,
            top_candidates=max(1, min(int(req.top), 50)),
            decoys=max(0, min(int(req.decoys), 100)),
            controls_csv=None,
            docking_backend=backend,
            docking_timeout=max(1, min(int(req.docking_timeout), 7200)),
            seed=61453,
            no_docking=not stage4_run_docking,
        ),
        {
            "validation": str(project_dir / "stage4_5" / f"round_{safe_round}_control_validation.json"),
            "scores": str(project_dir / "stage4_5" / f"round_{safe_round}_control_docking_scores.csv"),
        },
        warning_ok=True,
    )
    stage8_run_cli_step(
        steps,
        "stage46_benchmark",
        "Stage 4.6 retrospective benchmark",
        stage46_retrospective_benchmark,
        argparse.Namespace(
            project=str(project_dir),
            round=safe_round,
            positive_types="positive_control,reference_control,control",
            negative_types="decoy",
            top_k="1,3,5,10",
        ),
        {
            "benchmark": str(project_dir / "stage4_6" / f"round_{safe_round}_retrospective_benchmark.json"),
            "ranking": str(project_dir / "stage4_6" / f"round_{safe_round}_retrospective_ranking.csv"),
        },
        warning_ok=True,
    )

    title = req.title or f"{name} Stage 8 Demo Runner"
    for step_id, label, fn, args, files in [
        (
            "stage5_dashboard",
            "Stage 5 dashboard",
            stage5_dashboard,
            argparse.Namespace(project=str(project_dir), round=safe_round, title=title),
            {
                "dashboard_data": str(project_dir / "stage5" / "dashboard_data.json"),
                "dashboard_html": str(project_dir / "stage5" / "index.html"),
            },
        ),
        (
            "stage6_validation",
            "Stage 6 validation operations",
            stage6_validate,
            argparse.Namespace(project=str(project_dir), round=safe_round, top=max(1, min(int(req.top), 50))),
            {"assets": str(project_dir / "stage6" / f"round_{safe_round}_validation_assets.json")},
        ),
        (
            "stage7_delivery",
            "Stage 7 delivery package",
            stage7_package,
            argparse.Namespace(project=str(project_dir), round=safe_round, title=title),
            {"manifest": str(project_dir / "stage7" / f"round_{safe_round}_delivery_manifest.json")},
        ),
    ]:
        stage8_run_cli_step(steps, step_id, label, fn, args, files)

    command_center = stage8_payload(name, project_dir, safe_round)
    steps.append(stage8_step(
        "stage8_command_center",
        "Stage 8 command center",
        "completed",
        "aggregated current closed-loop state",
        {"stage5_dashboard_data": str(project_dir / "stage5" / "dashboard_data.json")},
        {"has_stage6": command_center.get("has_stage6"), "has_stage7": command_center.get("has_stage7")},
    ))

    ranked_rows = csv_to_dicts(project_dir / "ranked" / f"round_{safe_round}_ranked.csv")
    candidate_rows = csv_to_dicts(project_dir / "candidates" / f"round_{safe_round}_candidates.csv")
    failed = [step for step in steps if step.get("status") == "failed"]
    artifacts = {
        "target_selection": str(project_dir / "targets" / "target_selection.csv"),
        "target_report": str(project_dir / "reports" / "target_selection.md"),
        "target_pack_json": str(project_dir / "target_pack.json"),
        "stage2_matrix": str(project_dir / "evidence" / "stage2_target_sources.csv"),
        "raw_candidates": str(project_dir / "candidates" / f"round_{safe_round}_candidates.csv"),
        "candidate_scores": str(project_dir / "scores" / f"round_{safe_round}_scores.csv"),
        "ranked_candidates": str(project_dir / "ranked" / f"round_{safe_round}_ranked.csv"),
        "feedback": str(project_dir / "feedback" / f"round_{safe_round}_feedback.json"),
        "stage4_assets": str(project_dir / "stage4" / f"round_{safe_round}_stage4_assets.json"),
        "stage5_dashboard": str(project_dir / "stage5" / "index.html"),
        "stage6_assets": str(project_dir / "stage6" / f"round_{safe_round}_validation_assets.json"),
        "stage7_manifest": str(project_dir / "stage7" / f"round_{safe_round}_delivery_manifest.json"),
    }
    return {
        "stage": 8,
        "command": "stage8-demo-runner",
        "project": name,
        "round": safe_round,
        "status": "failed" if failed else "completed",
        "generated_at": datetime.now().isoformat(),
        "steps": steps,
        "summary": {
            "candidate_count": len(candidate_rows),
            "ranked_count": len(ranked_rows),
            "advanced_count": len([row for row in ranked_rows if row.get("decision") == "advance"]),
            "warning_count": len([step for step in steps if step.get("status") == "completed_with_warnings"]),
            "failed_count": len(failed),
        },
        "artifacts": artifacts,
        "stage8": command_center,
        "boundary": [
            "Demo Runner creates a computational closed-loop demo package only.",
            "It does not prove potency, efficacy, safety, dosing, toxicity, or clinical usefulness.",
            "Docking is executed only when explicitly enabled and local tools are available; otherwise docking-ready assets and validation plans are generated.",
        ],
    }

def stage8_file_check(step_id: str, label: str, path: Path, auto_repair: bool, next_action: str) -> Dict[str, Any]:
    exists = path.exists()
    return {
        "step_id": step_id,
        "label": label,
        "status": "ready" if exists else "missing",
        "evidence": str(path) if exists else f"missing: {path}",
        "auto_repair": bool(auto_repair),
        "next_action": "" if exists else next_action,
    }

def stage8_preflight_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    checks = [
        stage8_file_check("target_selection", "Stage 1 target selection", project_dir / "targets" / "target_selection.csv", True, "Run Stage 1 target selection."),
        stage8_file_check("target_pack", "Target Pack", project_dir / "target_pack.json", True, "Run Target Intake / Target Pack or use Demo Runner repair."),
        stage8_file_check("stage2_evidence", "Stage 2 evidence", project_dir / "evidence" / "stage2_target_sources.csv", True, "Generate Stage 2 evidence matrix."),
        stage8_file_check("stage3_assets", "Stage 3 candidate intake", project_dir / "stage3" / f"round_{safe_round}_stage3_assets.json", True, "Run unified candidate intake."),
        stage8_file_check("candidates", "Candidate CSV", project_dir / "candidates" / f"round_{safe_round}_candidates.csv", True, "Generate or import candidates."),
        stage8_file_check("scores", "Candidate scores", project_dir / "scores" / f"round_{safe_round}_scores.csv", True, "Run proxy scoring or import external scores."),
        stage8_file_check("ranked", "Ranked candidates", project_dir / "ranked" / f"round_{safe_round}_ranked.csv", True, "Run ranking and feedback."),
        stage8_file_check("feedback", "Feedback package", project_dir / "feedback" / f"round_{safe_round}_feedback.json", True, "Generate feedback seeds."),
        stage8_file_check("stage4_assets", "Stage 4 assets", project_dir / "stage4" / f"round_{safe_round}_stage4_assets.json", True, "Generate Stage 4 real-library assets."),
        stage8_file_check("stage5_dashboard", "Stage 5 dashboard", project_dir / "stage5" / "dashboard_data.json", True, "Generate Stage 5 dashboard."),
        stage8_file_check("stage6_validation", "Stage 6 validation", project_dir / "stage6" / f"round_{safe_round}_validation_assets.json", True, "Generate Stage 6 validation operations."),
        stage8_file_check("stage7_delivery", "Stage 7 delivery", project_dir / "stage7" / f"round_{safe_round}_delivery_manifest.json", True, "Generate Stage 7 delivery package."),
    ]
    missing = [item for item in checks if item["status"] == "missing"]
    return {
        "stage": 8,
        "command": "stage8-preflight",
        "project": name,
        "round": safe_round,
        "status": "ready" if not missing else "needs_repair",
        "checks": checks,
        "missing_count": len(missing),
        "can_auto_repair": any(item["auto_repair"] for item in missing),
        "next_actions": [item["next_action"] for item in missing if item.get("next_action")][:8] or ["Project is ready for Stage 8 review."],
        "boundary": "Preflight checks computational workflow artifacts only; it does not validate biological activity.",
    }

def stage8_demo_guide_step(
    step_id: str,
    stage: str,
    label: str,
    status: str,
    evidence: str,
    next_action: str,
    route: str,
    api: str,
) -> Dict[str, Any]:
    return {
        "step_id": step_id,
        "stage": stage,
        "label": label,
        "status": status,
        "evidence": evidence,
        "next_action": "" if status == "ready" else next_action,
        "route": route,
        "api": api,
    }

def stage8_paths_for_round(project_dir: Path, round_no: int) -> Dict[str, Path]:
    return {
        "config": project_dir / "config.json",
        "target_pack": project_dir / "target_pack.json",
        "stage2_matrix": project_dir / "evidence" / "stage2_target_sources.csv",
        "stage3_assets": project_dir / "stage3" / f"round_{round_no}_stage3_assets.json",
        "candidates": project_dir / "candidates" / f"round_{round_no}_candidates.csv",
        "scores": project_dir / "scores" / f"round_{round_no}_scores.csv",
        "ranked": project_dir / "ranked" / f"round_{round_no}_ranked.csv",
        "feedback": project_dir / "feedback" / f"round_{round_no}_feedback.json",
        "stage4_assets": project_dir / "stage4" / f"round_{round_no}_stage4_assets.json",
        "stage5_data": project_dir / "stage5" / "dashboard_data.json",
        "stage5_html": project_dir / "stage5" / "index.html",
        "stage6_assets": project_dir / "stage6" / f"round_{round_no}_validation_assets.json",
        "stage7_manifest": project_dir / "stage7" / f"round_{round_no}_delivery_manifest.json",
        "stage8_report": project_dir / "reports" / f"stage8_round_{round_no}_closed_loop_report.md",
        "evidence_pack_zip": project_dir / "exports" / f"round_{round_no}_evidence_pack.zip",
    }

def stage8_demo_guide_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    paths = stage8_paths_for_round(project_dir, safe_round)
    stage3_core = [paths["candidates"], paths["scores"], paths["ranked"], paths["feedback"]]
    stage3_ready = all(path.exists() for path in stage3_core)
    stage3_status = "ready" if stage3_ready and paths["stage3_assets"].exists() else ("warn" if stage3_ready else "missing")
    stage5_status = "ready" if paths["stage5_data"].exists() and paths["stage5_html"].exists() else ("warn" if paths["stage5_data"].exists() else "missing")

    steps = [
        stage8_demo_guide_step(
            "project_setup",
            "0",
            "选择或创建项目",
            "ready" if paths["config"].exists() else "missing",
            str(paths["config"]) if paths["config"].exists() else "missing config.json",
            "先在项目登记页创建项目，或选择已有 Demo 项目。",
            "projects",
            "POST /api/projects",
        ),
        stage8_demo_guide_step(
            "target_pack",
            "1",
            "生成 Target Pack",
            "ready" if paths["target_pack"].exists() else "missing",
            str(paths["target_pack"]) if paths["target_pack"].exists() else "missing target_pack.json",
            "在靶点需求页生成 Target Pack，或运行 Stage 8 自动修复/Demo Runner。",
            "brief",
            f"POST /api/projects/{name}/target-pack",
        ),
        stage8_demo_guide_step(
            "target_evidence",
            "2",
            "确认靶点证据矩阵",
            "ready" if paths["stage2_matrix"].exists() else "missing",
            str(paths["stage2_matrix"]) if paths["stage2_matrix"].exists() else "missing stage2_target_sources.csv",
            "运行 Stage 2 证据矩阵，确认靶点、PDB、控药和公开证据。",
            "stage2",
            f"POST /api/projects/{name}/stage2/evidence",
        ),
        stage8_demo_guide_step(
            "candidate_intake",
            "3",
            "候选生成与入口统一",
            stage3_status,
            "candidate/scores/ranked/feedback present" if stage3_ready else "missing candidate scoring assets",
            "使用统一候选入口生成或导入 SMILES，并完成评分、排序和反馈种子。",
            "stage8",
            f"POST /api/projects/{name}/stage3/candidates",
        ),
        stage8_demo_guide_step(
            "stage4_validation",
            "4",
            "真实库校验与对照面板",
            "ready" if paths["stage4_assets"].exists() else "missing",
            str(paths["stage4_assets"]) if paths["stage4_assets"].exists() else "missing Stage 4 assets",
            "运行 Stage 4，生成 RDKit 描述符、控药/decoy 面板、docking plan 和可视化资产。",
            "stage4",
            f"POST /api/projects/{name}/stage4/real",
        ),
        stage8_demo_guide_step(
            "stage5_dashboard",
            "5",
            "Dashboard 看板",
            stage5_status,
            "dashboard data/html present" if stage5_status == "ready" else "missing Stage 5 dashboard",
            "生成 Stage 5 Dashboard，让候选漏斗、质量和靶点摘要进入可演示状态。",
            "dashboard",
            f"POST /api/projects/{name}/stage5/dashboard",
        ),
        stage8_demo_guide_step(
            "stage6_validation",
            "6",
            "验证运营与质量门",
            "ready" if paths["stage6_assets"].exists() else "missing",
            str(paths["stage6_assets"]) if paths["stage6_assets"].exists() else "missing Stage 6 validation assets",
            "运行 Stage 6，生成质量门、命中分层、验证队列和风险登记。",
            "stage6",
            f"POST /api/projects/{name}/stage6/validate",
        ),
        stage8_demo_guide_step(
            "stage7_delivery",
            "7",
            "交付包与复现说明",
            "ready" if paths["stage7_manifest"].exists() else "missing",
            str(paths["stage7_manifest"]) if paths["stage7_manifest"].exists() else "missing Stage 7 delivery manifest",
            "运行 Stage 7，生成交付清单、执行摘要、复现说明和前端规格。",
            "stage7",
            f"POST /api/projects/{name}/stage7/package",
        ),
        stage8_demo_guide_step(
            "stage8_report",
            "8",
            "生成闭环报告",
            "ready" if paths["stage8_report"].exists() else "missing",
            str(paths["stage8_report"]) if paths["stage8_report"].exists() else "missing Stage 8 closed-loop report",
            "生成 Stage 8 闭环报告，汇总靶点、候选、Stage 4、质量门和边界声明。",
            "stage8",
            f"POST /api/projects/{name}/stage8/report",
        ),
        stage8_demo_guide_step(
            "evidence_pack",
            "8",
            "导出一键证据包",
            "ready" if paths["evidence_pack_zip"].exists() else "missing",
            str(paths["evidence_pack_zip"]) if paths["evidence_pack_zip"].exists() else "missing evidence pack zip",
            "导出证据包 ZIP，用于答辩、前端交接或材料归档。",
            "stage8",
            f"POST /api/projects/{name}/stage8/evidence-pack",
        ),
    ]
    counts = {
        "ready": len([step for step in steps if step["status"] == "ready"]),
        "warn": len([step for step in steps if step["status"] == "warn"]),
        "missing": len([step for step in steps if step["status"] == "missing"]),
    }
    total = len(steps)
    next_step = next((step for step in steps if step["status"] != "ready"), None)
    if not next_step:
        next_step = {
            "step_id": "present_demo",
            "stage": "8",
            "label": "演示或下载证据包",
            "status": "ready",
            "evidence": "all guided demo steps are ready",
            "next_action": "打开 Stage 8 指挥台，按向导顺序展示并下载证据包。",
            "route": "stage8",
            "api": f"POST /api/projects/{name}/stage8/evidence-pack",
        }
    return {
        "stage": 8,
        "command": "stage8-demo-guide",
        "product_mode": "guided_demo",
        "project": name,
        "round": safe_round,
        "generated_at": datetime.now().isoformat(),
        "steps": steps,
        "completion": {
            **counts,
            "total": total,
            "percent": round((counts["ready"] / total) * 100, 1) if total else 0,
        },
        "next_primary_action": next_step,
        "stage8": stage8_payload(name, project_dir, safe_round),
        "boundary": [
            "Guided demo steps organize computational screening artifacts for presentation and handoff.",
            "They do not prove biological activity, efficacy, toxicity, safety, dosing, or clinical usefulness.",
        ],
    }

def render_stage8_closed_loop_report(name: str, project_dir: Path, round_no: int, title: str = "") -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    stage8 = stage8_payload(name, project_dir, safe_round)
    target_pack = read_json_if_exists(project_dir / "target_pack.json")
    stage3_assets = read_json_if_exists(project_dir / "stage3" / f"round_{safe_round}_stage3_assets.json")
    stage4_assets = read_json_if_exists(project_dir / "stage4" / f"round_{safe_round}_stage4_assets.json")
    ranked = csv_to_dicts(project_dir / "ranked" / f"round_{safe_round}_ranked.csv")
    gates = stage8.get("quality_gates", []) if isinstance(stage8.get("quality_gates"), list) else []
    target = stage8.get("target", {}) if isinstance(stage8.get("target"), dict) else {}
    metrics = stage8.get("metrics", {}) if isinstance(stage8.get("metrics"), dict) else {}
    title = title or f"{name} Stage 8 Closed Loop Report"
    sections = {
        "target_rationale": (
            f"Target `{target_pack.get('target_id') or target.get('id', '')}` was selected for `{target.get('disease_context', '')}` "
            f"with PDB `{target_pack.get('pdb_id') or target.get('primary_pdb', '')}` and reference ligand `{target_pack.get('reference_ligand') or target.get('reference_ligand', '')}`."
        ),
        "candidate_generation": f"Raw candidates={metrics.get('raw_candidates', 0)}, filtered={metrics.get('filtered_candidates', 0)}, ranked={metrics.get('ranked_candidates', 0)}.",
        "scoring_validation": f"Stage 3 sources={', '.join(stage3_assets.get('sources', []) or [])}; Stage 4 status={stage4_assets.get('docking_plan', {}).get('status', stage4_assets.get('docking_status', 'not_available'))}.",
        "quality_gates": f"Quality gates={len(gates)}; delivery={stage8.get('delivery_status', '') or 'not_generated'}.",
        "claim_boundary": "This is a computational screening and validation-planning report, not biological efficacy proof.",
    }
    lines = [
        f"# {title}",
        "",
        f"- Project: `{name}`",
        f"- Round: `{safe_round}`",
        f"- Generated: `{datetime.now().isoformat()}`",
        "",
        "## Target Rationale",
        "",
        sections["target_rationale"],
        "",
        "## Candidate Generation And Funnel",
        "",
        sections["candidate_generation"],
        "",
        "## Top Ranked Candidates",
        "",
        "| rank | id | smiles | total | decision |",
        "|---:|---|---|---:|---|",
    ]
    for row in ranked[:10]:
        lines.append(f"| {row.get('rank', '')} | {row.get('id', '')} | `{row.get('smiles', '')}` | {row.get('total_proxy', '')} | {row.get('decision', '')} |")
    lines.extend(
        [
            "",
            "## Scoring And Validation Evidence",
            "",
            sections["scoring_validation"],
            "",
            "## Quality Gates",
            "",
            sections["quality_gates"],
            "",
            "## Claim Boundary",
            "",
            "- Computational screening and validation planning only.",
            "- No potency, efficacy, toxicity, dosing, safety, clinical benefit, or therapeutic claim is created.",
            "- Known drugs are controls and references, not newly discovered outputs.",
        ]
    )
    report_path = project_dir / "reports" / f"stage8_round_{safe_round}_closed_loop_report.md"
    write_text(report_path, "\n".join(lines) + "\n")
    return {
        "stage": 8,
        "command": "stage8-report",
        "project": name,
        "round": safe_round,
        "title": title,
        "sections": sections,
        "content": "\n".join(lines) + "\n",
        "files": {"report": str(report_path)},
        "boundary": [
            "Computational screening report only.",
            "No efficacy, potency, toxicity, safety, dosing, clinical benefit, or therapeutic claim.",
        ],
    }

def catalog_target_summary(row: Dict[str, Any]) -> Dict[str, Any]:
    structures = row.get("representative_structures", [])
    first = structures[0] if isinstance(structures, list) and structures and isinstance(structures[0], dict) else {}
    site = row.get("binding_site", {}) if isinstance(row.get("binding_site"), dict) else {}
    scores = row.get("scores", {}) if isinstance(row.get("scores"), dict) else {}
    return {
        "id": row.get("id", ""),
        "display_name": row.get("display_name", ""),
        "short_name": row.get("short_name", ""),
        "recommendation": row.get("recommendation", ""),
        "target_family": row.get("target_family", ""),
        "pathogen_scope": "; ".join(str(x) for x in row.get("pathogen_scope", []) if str(x).strip()) if isinstance(row.get("pathogen_scope"), list) else "",
        "known_drugs": "; ".join(str(x) for x in row.get("known_drugs", []) if str(x).strip()) if isinstance(row.get("known_drugs"), list) else "",
        "primary_pdb": first.get("pdb_id", ""),
        "reference_ligand": site.get("reference_ligand", first.get("ligand", "")),
        "evidence_score": scores.get("clinical_validation", scores.get("disease_relevance", "")),
    }

def custom_target_from_request(req: TargetCatalogCustomRequest) -> Dict[str, Any]:
    pdb_id = str(req.pdb_id or "").strip().upper()
    return {
        "id": safe_name(req.target_id),
        "display_name": req.display_name,
        "short_name": req.display_name,
        "pathogen_scope": split_items(req.disease),
        "target_family": "custom_target",
        "mechanism": req.mechanism or "User supplied target; requires evidence review.",
        "recommendation": req.recommendation or "custom_review",
        "recommendation_reason": "User supplied custom target for project-level exploration.",
        "known_drugs": split_items(req.known_drugs),
        "representative_structures": ([{"pdb_id": pdb_id, "description": "User supplied structure", "ligand": req.reference_ligand, "use": "custom_review"}] if pdb_id else []),
        "binding_site": {
            "source": "user_supplied",
            "description": "User supplied binding-site placeholder; review before docking.",
            "reference_ligand": req.reference_ligand,
            "key_residues": [],
            "box_strategy": "Manual review required before real docking.",
        },
        "scores": {
            "disease_relevance": 0.5,
            "clinical_validation": 0.0,
            "structure_availability": 0.6 if pdb_id else 0.0,
            "pocket_confidence": 0.2,
            "known_ligand_data": 0.4 if req.reference_ligand else 0.0,
            "assay_feasibility": 0.2,
            "computational_feasibility": 0.4,
        },
        "source_urls": [],
    }

def cross_disease_catalog_targets() -> List[Dict[str, Any]]:
    return [
        {
            "id": "egfr_tyrosine_kinase",
            "display_name": "EGFR tyrosine kinase domain",
            "short_name": "EGFR TK",
            "pathogen_scope": ["oncology", "lung cancer", "solid tumor", "肿瘤", "肺癌"],
            "target_family": "human receptor tyrosine kinase",
            "mechanism": "EGFR kinase signaling drives subsets of non-small-cell lung cancer and other tumors.",
            "recommendation": "cross_disease_demo",
            "recommendation_reason": "Public structures, approved inhibitors, and clear ATP-site docking setup make EGFR a practical non-flu expansion target.",
            "known_drugs": ["gefitinib", "erlotinib", "afatinib", "osimertinib"],
            "representative_structures": [
                {"pdb_id": "1M17", "description": "EGFR kinase domain inhibitor complex", "ligand": "erlotinib-like inhibitor", "use": "oncology_demo_structure"},
                {"pdb_id": "4ZAU", "description": "EGFR T790M inhibitor complex family example", "ligand": "osimertinib analog", "use": "resistance_context"},
            ],
            "binding_site": {
                "source": "co_crystal",
                "description": "ATP-binding kinase hinge pocket with known small-molecule inhibitors.",
                "reference_ligand": "erlotinib",
                "key_residues": ["hinge region", "gatekeeper residue", "DFG region"],
                "box_strategy": "Use co-crystallized inhibitor centroid from a reviewed EGFR kinase complex.",
            },
            "scores": {"disease_relevance": 0.92, "clinical_validation": 0.96, "structure_availability": 0.90, "pocket_confidence": 0.88, "known_ligand_data": 0.92, "assay_feasibility": 0.78, "computational_feasibility": 0.84},
            "source_urls": ["https://www.rcsb.org/structure/1M17"],
        },
        {
            "id": "bace1_amyloid_beta_secretase",
            "display_name": "BACE1 beta-secretase",
            "short_name": "BACE1",
            "pathogen_scope": ["neurodegeneration", "alzheimer", "阿尔茨海默"],
            "target_family": "human aspartyl protease",
            "mechanism": "BACE1 cleaves APP and is a historically explored Alzheimer's disease drug target.",
            "recommendation": "cross_disease_demo",
            "recommendation_reason": "Rich public structural data and inhibitor controls are useful for method demonstration, while clinical translation should be framed cautiously.",
            "known_drugs": ["verubecestat", "lanabecestat", "atabecestat"],
            "representative_structures": [
                {"pdb_id": "2ZHV", "description": "BACE1 inhibitor complex family example", "ligand": "BACE1 inhibitor", "use": "neurodegeneration_demo_structure"}
            ],
            "binding_site": {
                "source": "co_crystal",
                "description": "Aspartyl protease active-site cleft with inhibitor co-crystal references.",
                "reference_ligand": "BACE1 inhibitor",
                "key_residues": ["catalytic Asp dyad", "flap region"],
                "box_strategy": "Use co-crystallized inhibitor centroid and review protonation state before docking.",
            },
            "scores": {"disease_relevance": 0.78, "clinical_validation": 0.55, "structure_availability": 0.92, "pocket_confidence": 0.88, "known_ligand_data": 0.86, "assay_feasibility": 0.72, "computational_feasibility": 0.82},
            "source_urls": ["https://www.rcsb.org/structure/2ZHV"],
        },
        {
            "id": "hiv1_protease",
            "display_name": "HIV-1 protease",
            "short_name": "HIV-1 PR",
            "pathogen_scope": ["hiv", "antiviral", "viral protease", "艾滋病"],
            "target_family": "viral aspartyl protease",
            "mechanism": "HIV-1 protease processes viral polyproteins and has many approved inhibitor controls.",
            "recommendation": "cross_disease_demo",
            "recommendation_reason": "Clinically validated target with abundant co-crystal structures and known inhibitors for benchmarking.",
            "known_drugs": ["saquinavir", "ritonavir", "lopinavir", "darunavir", "atazanavir"],
            "representative_structures": [
                {"pdb_id": "1HSG", "description": "HIV-1 protease inhibitor complex", "ligand": "indinavir-like inhibitor", "use": "antiviral_demo_structure"}
            ],
            "binding_site": {
                "source": "co_crystal",
                "description": "Dimeric aspartyl protease active site with flap dynamics; use controls and decoys carefully.",
                "reference_ligand": "protease inhibitor",
                "key_residues": ["Asp25", "Asp25'", "flap residues"],
                "box_strategy": "Use inhibitor centroid from a reviewed protease complex and preserve dimer context.",
            },
            "scores": {"disease_relevance": 0.95, "clinical_validation": 0.96, "structure_availability": 0.95, "pocket_confidence": 0.90, "known_ligand_data": 0.94, "assay_feasibility": 0.76, "computational_feasibility": 0.86},
            "source_urls": ["https://www.rcsb.org/structure/1HSG"],
        },
    ]

def expanded_target_catalog(query: str = "") -> Dict[str, Any]:
    catalog = load_target_catalog(None)
    base_targets = catalog.get("targets", []) if isinstance(catalog.get("targets"), list) else []
    merged = list(base_targets) + cross_disease_catalog_targets()
    expanded = dict(catalog)
    expanded["domain"] = "multi_disease_drug_discovery_targets"
    expanded["updated"] = "2026-05-25"
    expanded["targets"] = merged
    expanded["_expansion_note"] = "Includes influenza MVP targets plus oncology, neurodegeneration, and HIV demonstration targets."
    return expanded

def stage3_candidates_payload(name: str, project_dir: Path, req: Stage3CandidateRequest) -> Dict[str, Any]:
    ensure_project_dirs(project_dir)
    safe_round = max(1, int(req.round))
    source_json = None
    if str(req.source_mode or "").lower() == "text" and str(req.source_text or "").strip():
        source_path = project_dir / "stage3" / f"round_{safe_round}_source_text.txt"
        write_text(source_path, req.source_text)
        source_json = str(source_path)
    try:
        stage3_screen(
            argparse.Namespace(
                project=str(project_dir),
                round=safe_round,
                n=max(1, min(int(req.n), 200)),
                top=max(1, min(int(req.top), 50)),
                source_csv=req.source_csv if str(req.source_mode or "").lower() == "csv" else None,
                source_json=source_json,
                source_url=req.source_url or [],
                context_url=req.context_url or [],
                use_openai=bool(req.use_openai),
                openai_model=req.openai_model or "gpt-5.2",
                api_key=req.api_key,
                api_key_file=req.api_key_file,
                api_key_env=req.api_key_env or "OPENAI_API_KEY",
                prompt=req.prompt or "",
                timeout=max(1, min(int(req.timeout), 180)),
                max_url_bytes=max(1000, min(int(req.max_url_bytes), 2_000_000)),
                max_heavy_atoms=max(1, min(int(req.max_heavy_atoms), 200)),
                max_molecular_weight=max(100.0, min(float(req.max_molecular_weight), 2000.0)),
                max_lipinski_violations=max(0, min(int(req.max_lipinski_violations), 10)),
                allow_risk=bool(req.allow_risk),
                external_scores=req.external_scores,
                no_score=bool(req.no_score),
            )
        )
    except SystemExit as exc:
        http_from_cli_exit(exc)
    assets_path = project_dir / "stage3" / f"round_{safe_round}_stage3_assets.json"
    return {
        "stage": 3,
        "command": "stage3-candidates",
        "project": name,
        "round": safe_round,
        "assets": read_json_if_exists(assets_path),
        "api_key_security": {
            "api_key_input_received": bool(str(req.api_key or "").strip()),
            "api_key_file_requested": bool(str(req.api_key_file or "").strip()),
            "api_key_env": req.api_key_env or "OPENAI_API_KEY",
            "api_key_persisted": False,
            "policy": "OpenAI API keys are used in memory for the current request and are not written to project outputs.",
        },
        "raw_candidates": csv_to_dicts(project_dir / "stage3" / f"round_{safe_round}_raw_candidates.csv"),
        "filtered": csv_to_dicts(project_dir / "filtered" / f"round_{safe_round}_filtered.csv"),
        "candidates": csv_to_dicts(project_dir / "candidates" / f"round_{safe_round}_candidates.csv"),
        "ranked": csv_to_dicts(project_dir / "ranked" / f"round_{safe_round}_ranked.csv"),
        "files": {
            "assets": str(assets_path),
            "raw_candidates": str(project_dir / "stage3" / f"round_{safe_round}_raw_candidates.csv"),
            "filtered": str(project_dir / "filtered" / f"round_{safe_round}_filtered.csv"),
            "candidates": str(project_dir / "candidates" / f"round_{safe_round}_candidates.csv"),
            "ranked": str(project_dir / "ranked" / f"round_{safe_round}_ranked.csv"),
            "report": str(project_dir / "reports" / f"stage3_round_{safe_round}_report.md"),
        },
        "boundary": "Candidate generation/intake is computational-only and does not imply biological activity.",
    }

def stage3_status_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    assets_path = project_dir / "stage3" / f"round_{safe_round}_stage3_assets.json"
    report_path = project_dir / "reports" / f"stage3_round_{safe_round}_report.md"
    assets = read_json_if_exists(assets_path)
    openai = assets.get("openai", {}) if isinstance(assets.get("openai"), dict) else {}
    prompt_path = str(openai.get("prompt_file") or "")
    return {
        "stage": 3,
        "command": "stage3-status",
        "project": name,
        "round": safe_round,
        "has_assets": assets_path.exists(),
        "assets": assets,
        "api_key_security": {
            "api_key_input_received": False,
            "api_key_file_requested": False,
            "api_key_env": "OPENAI_API_KEY",
            "api_key_persisted": False,
            "policy": "OpenAI API keys are used in memory for the current request and are not written to project outputs.",
        },
        "raw_candidates": csv_to_dicts(project_dir / "stage3" / f"round_{safe_round}_raw_candidates.csv"),
        "filtered": csv_to_dicts(project_dir / "filtered" / f"round_{safe_round}_filtered.csv"),
        "candidates": csv_to_dicts(project_dir / "candidates" / f"round_{safe_round}_candidates.csv"),
        "ranked": csv_to_dicts(project_dir / "ranked" / f"round_{safe_round}_ranked.csv"),
        "feedback": read_json_if_exists(project_dir / "feedback" / f"round_{safe_round}_feedback.json"),
        "report": read_text_if_exists(report_path),
        "files": {
            "stage3_assets": str(assets_path),
            "assets": str(assets_path),
            "raw_candidates": str(project_dir / "stage3" / f"round_{safe_round}_raw_candidates.csv"),
            "filtered": str(project_dir / "filtered" / f"round_{safe_round}_filtered.csv"),
            "candidates": str(project_dir / "candidates" / f"round_{safe_round}_candidates.csv"),
            "ranked": str(project_dir / "ranked" / f"round_{safe_round}_ranked.csv"),
            "feedback": str(project_dir / "feedback" / f"round_{safe_round}_feedback.json"),
            "report": str(report_path),
            "prompt_file": prompt_path,
        },
        "boundary": "Candidate generation/intake is computational screening only and does not imply biological activity.",
    }

def ensure_stage3_assets_from_candidates(project_dir: Path, round_no: int, source_label: str = "stage8_demo_runner") -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    candidates = csv_to_dicts(project_dir / "candidates" / f"round_{safe_round}_candidates.csv")
    raw_path = project_dir / "stage3" / f"round_{safe_round}_raw_candidates.csv"
    filtered_path = project_dir / "filtered" / f"round_{safe_round}_filtered.csv"
    assets_path = project_dir / "stage3" / f"round_{safe_round}_stage3_assets.json"
    report_path = project_dir / "reports" / f"stage3_round_{safe_round}_report.md"
    if candidates and not raw_path.exists():
        write_csv(raw_path, candidates, ["round", "id", "smiles", "parent", "source", "rationale", "expected_interaction", "design_family", "risk_note"])
    if candidates and not filtered_path.exists():
        write_csv(filtered_path, candidates, ["round", "id", "smiles", "parent", "source", "rationale", "expected_interaction", "design_family", "risk_note"])
    assets = {
        "schema_version": "0.1",
        "generated_at": datetime.now().isoformat(),
        "round": safe_round,
        "sources": [source_label],
        "raw_candidates": len(candidates),
        "passed_filter": len(candidates),
        "failed_filter": 0,
        "openai": {"used": False, "model": "", "api_key_source": "not_used"},
        "files": {
            "raw_candidates": str(raw_path),
            "filtered_candidates": str(filtered_path),
            "scoring_candidates": str(project_dir / "candidates" / f"round_{safe_round}_candidates.csv"),
            "ranked": str(project_dir / "ranked" / f"round_{safe_round}_ranked.csv"),
        },
    }
    write_json(assets_path, assets)
    if not report_path.exists():
        write_text(
            report_path,
            "\n".join(
                [
                    f"# Stage 3 Candidate Intake · Round {safe_round}",
                    "",
                    f"- Source: `{source_label}`",
                    f"- Raw candidates: `{len(candidates)}`",
                    f"- Passed filter: `{len(candidates)}`",
                    "",
                    "## Boundary",
                    "",
                    "- This generated Stage 3 asset record documents the candidate intake files used by the product demo.",
                    "- It is computational-only and does not imply biological activity.",
                ]
            )
            + "\n",
        )
    return assets

def stage8_acceptance_check(check_id: str, label: str, passed: bool, evidence: str, next_action: str = "") -> Dict[str, Any]:
    return {
        "check_id": check_id,
        "label": label,
        "status": "passed" if passed else "failed",
        "evidence": evidence,
        "next_action": "" if passed else next_action,
    }

def render_stage8_acceptance_report(payload: Dict[str, Any]) -> str:
    lines = [
        f"# Stage 8 Acceptance Report · {payload.get('project', '')}",
        "",
        f"- Status: `{payload.get('status', '')}`",
        f"- Round: `{payload.get('round', '')}`",
        f"- Target: `{payload.get('target', {}).get('id', '')}`",
        f"- Generated: `{payload.get('generated_at', '')}`",
        "",
        "## Checks",
        "",
        "| check | status | evidence |",
        "|---|---|---|",
    ]
    for row in payload.get("checks", []):
        lines.append(f"| {row.get('label', '')} | {row.get('status', '')} | {row.get('evidence', '')} |")
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- Acceptance confirms that the computational demo workflow is reproducible from the web product surface.",
            "- It does not prove potency, efficacy, toxicity, dosing, safety, clinical benefit, or therapeutic usefulness.",
        ]
    )
    return "\n".join(lines) + "\n"

def create_stage8_acceptance_demo(req: Stage8AcceptanceDemoRequest) -> Dict[str, Any]:
    safe_name_value = safe_name(req.project_name or "flu_na_acceptance_demo")
    if not safe_name_value:
        raise HTTPException(400, "项目名不能为空")
    project_dir = PROJECTS_ROOT / safe_name_value
    if project_dir.exists() and not req.force:
        raise HTTPException(409, f"项目 {safe_name_value} 已存在，使用 force=true 覆盖")
    if project_dir.exists() and req.force:
        import shutil
        shutil.rmtree(project_dir)
    ensure_project_dirs(project_dir)
    cfg = default_config()
    target_cfg = cfg.get("target", {})
    if isinstance(target_cfg, dict):
        target_cfg["name"] = safe_name_value
    write_json(project_dir / "config.json", cfg)
    seed_rows = [{"id": sid, "smiles": smi, "note": note} for sid, smi, note in DEFAULT_SEEDS]
    write_csv(project_dir / "seeds" / "round_0_seeds.csv", seed_rows, ["id", "smiles", "note"])

    safe_round = max(1, int(req.round))
    runner_req = Stage8DemoRunnerRequest(
        round=safe_round,
        disease="甲流",
        target="influenza_a_h1n1_na",
        source_kind="url",
        source="https://www.rcsb.org/structure/3TI6",
        target_hint="neuraminidase oseltamivir pocket",
        candidates=max(1, min(int(req.candidates), 200)),
        top=max(1, min(int(req.top), 50)),
        decoys=max(0, min(int(req.decoys), 200)),
        run_docking=bool(req.run_docking),
        docking_backend="auto",
        title=f"{safe_name_value} Stage 8 Acceptance Demo",
    )
    runner = stage8_demo_runner(safe_name_value, project_dir, runner_req)
    ensure_stage3_assets_from_candidates(project_dir, safe_round, "stage8_acceptance_demo")
    closed_loop_report = render_stage8_closed_loop_report(
        safe_name_value,
        project_dir,
        safe_round,
        f"{safe_name_value} Stage 8 Acceptance Closed Loop Report",
    )
    preflight = stage8_preflight_payload(safe_name_value, project_dir, safe_round)
    command_center = stage8_payload(safe_name_value, project_dir, safe_round)
    target = command_center.get("target", {}) if isinstance(command_center.get("target"), dict) else {}
    downloads = command_center.get("download_links", {}) if isinstance(command_center.get("download_links"), dict) else {}
    candidate_rows = csv_to_dicts(project_dir / "candidates" / f"round_{safe_round}_candidates.csv")
    target_selection_rows = csv_to_dicts(project_dir / "targets" / "target_selection.csv")
    catalog = load_target_catalog(None)
    catalog_targets = sorted_targets(catalog, "甲流")
    target_ids = {str(row.get("id", "")) for row in catalog_targets}
    checks = [
        stage8_acceptance_check("target_catalog", "靶点目录", "influenza_a_h1n1_na" in target_ids, "influenza_a_h1n1_na present in catalog query 甲流", "Check influenza target catalog."),
        stage8_acceptance_check("target_selection", "靶点筛选", bool(target_selection_rows) and target_selection_rows[0].get("target_id") == "influenza_a_h1n1_na", str(project_dir / "targets" / "target_selection.csv"), "Run Stage 1 target selection."),
        stage8_acceptance_check("target_pack", "Target Pack", (project_dir / "target_pack.json").exists(), str(project_dir / "target_pack.json"), "Run Target Pack."),
        stage8_acceptance_check("candidate_intake", "候选入口", bool(candidate_rows) and (project_dir / "stage3" / f"round_{safe_round}_stage3_assets.json").exists(), f"{len(candidate_rows)} candidates", "Run unified candidate intake."),
        stage8_acceptance_check("stage4_assets", "Stage 4 真实库资产", (project_dir / "stage4" / f"round_{safe_round}_stage4_assets.json").exists(), str(project_dir / "stage4" / f"round_{safe_round}_stage4_assets.json"), "Run Stage 4."),
        stage8_acceptance_check("stage5_dashboard", "Stage 5 Dashboard", (project_dir / "stage5" / "dashboard_data.json").exists() and (project_dir / "stage5" / "index.html").exists(), str(project_dir / "stage5" / "index.html"), "Run Stage 5."),
        stage8_acceptance_check("stage6_validation", "Stage 6 验证运营", (project_dir / "stage6" / f"round_{safe_round}_validation_assets.json").exists(), str(project_dir / "stage6" / f"round_{safe_round}_validation_assets.json"), "Run Stage 6."),
        stage8_acceptance_check("stage7_delivery", "Stage 7 交付包", (project_dir / "stage7" / f"round_{safe_round}_delivery_manifest.json").exists(), str(project_dir / "stage7" / f"round_{safe_round}_delivery_manifest.json"), "Run Stage 7."),
        stage8_acceptance_check("stage8_preflight", "Stage 8 预检", preflight.get("status") == "ready", f"status={preflight.get('status')}", "Run Stage 8 repair."),
        stage8_acceptance_check("closed_loop_report", "闭环报告", bool(closed_loop_report.get("files", {}).get("report")) and Path(str(closed_loop_report.get("files", {}).get("report"))).exists(), str(closed_loop_report.get("files", {}).get("report", "")), "Generate Stage 8 report."),
        stage8_acceptance_check("download_materials", "下载材料", bool(downloads), f"{len(downloads)} download links", "Generate Stage 7/8 package."),
    ]
    status = "passed" if all(row["status"] == "passed" for row in checks) else "failed"
    payload = {
        "stage": 8,
        "command": "stage8-acceptance-demo",
        "project": safe_name_value,
        "round": safe_round,
        "status": status,
        "generated_at": datetime.now().isoformat(),
        "target": target,
        "checks": checks,
        "runner": runner,
        "preflight": preflight,
        "stage8": command_center,
        "files": {
            "acceptance_json": str(project_dir / "acceptance" / "stage8_acceptance_report.json"),
            "acceptance_md": str(project_dir / "acceptance" / "stage8_acceptance_report.md"),
            "closed_loop_report": str(closed_loop_report.get("files", {}).get("report", "")),
        },
        "boundary": [
            "Acceptance verifies product workflow reproducibility and artifact availability only.",
            "It does not prove potency, efficacy, toxicity, safety, dosing, or clinical usefulness.",
        ],
    }
    report_md = render_stage8_acceptance_report(payload)
    write_json(project_dir / "acceptance" / "stage8_acceptance_report.json", payload)
    write_text(project_dir / "acceptance" / "stage8_acceptance_report.md", report_md)
    payload["report_preview"] = report_md
    return payload

# ── 项目管理 ──────────────────────────────────────────────────────────────
@app.get("/api/projects")
def api_list_projects(): return list_projects()

@app.post("/api/projects")
def api_create_project(req: InitProjectRequest):
    if not req.name: raise HTTPException(400, "项目名不能为空")
    safe = safe_name(req.name)
    project_dir = PROJECTS_ROOT / safe
    if project_dir.exists() and not req.force:
        raise HTTPException(409, f"项目 {safe} 已存在，使用 force=true 覆盖")
    ensure_project_dirs(project_dir)
    cfg = default_config()
    target = cfg.get("target", {})
    if isinstance(target, dict): target["name"] = safe
    write_json(project_dir / "config.json", cfg)
    seed_rows = [{"id": sid, "smiles": smi, "note": note} for sid, smi, note in DEFAULT_SEEDS]
    write_csv(project_dir / "seeds" / "round_0_seeds.csv", seed_rows, ["id", "smiles", "note"])
    return {"name": safe, "message": "项目创建成功"}

@app.delete("/api/projects/{name}")
def api_delete_project(name: str):
    import shutil; shutil.rmtree(get_project_dir(name))
    return {"message": f"项目 {name} 已删除"}

@app.get("/api/projects/{name}")
def api_get_project(name: str):
    p = get_project_dir(name)
    cfg = {}; config_path = p / "config.json"
    if config_path.exists(): cfg = json.loads(config_path.read_text())
    return {"name": name, "config": cfg, "rounds": get_round_files(p),
            "has_brief": (p / "briefs" / "target_brief.json").exists()}

@app.get("/api/projects/{name}/artifact")
def api_project_artifact(name: str, path: str = Query(...)):
    p = get_project_dir(name)
    artifact = safe_project_artifact_path(p, path)
    return FileResponse(artifact)

# ── 配置 ──────────────────────────────────────────────────────────────────
@app.get("/api/projects/{name}/config")
def api_get_config(name: str):
    config_path = get_project_dir(name) / "config.json"
    if not config_path.exists(): raise HTTPException(404, "配置文件不存在")
    return json.loads(config_path.read_text())

@app.put("/api/projects/{name}/config")
def api_update_config(name: str, req: ConfigUpdateRequest):
    write_json(get_project_dir(name) / "config.json", req.config)
    return {"message": "配置已更新"}

# ── 数据查询 ──────────────────────────────────────────────────────────────
@app.get("/api/projects/{name}/candidates/{round_no}")
def api_get_candidates(name: str, round_no: int):
    return csv_to_dicts(get_project_dir(name) / "candidates" / f"round_{round_no}_candidates.csv")

@app.get("/api/projects/{name}/scores/{round_no}")
def api_get_scores(name: str, round_no: int):
    return csv_to_dicts(get_project_dir(name) / "scores" / f"round_{round_no}_scores.csv")

@app.get("/api/projects/{name}/ranked/{round_no}")
def api_get_ranked(name: str, round_no: int):
    return csv_to_dicts(get_project_dir(name) / "ranked" / f"round_{round_no}_ranked.csv")

@app.get("/api/projects/{name}/feedback/{round_no}")
def api_get_feedback(name: str, round_no: int):
    p = get_project_dir(name) / "feedback" / f"round_{round_no}_feedback.json"
    if not p.exists(): raise HTTPException(404, "反馈文件不存在")
    return json.loads(p.read_text())

@app.get("/api/projects/{name}/report/{round_no}")
def api_get_report(name: str, round_no: int):
    p = get_project_dir(name) / "reports" / f"round_{round_no}_summary.md"
    if not p.exists(): raise HTTPException(404, "报告不存在")
    return {"content": p.read_text(encoding="utf-8")}

# ── 靶点需求 Brief ───────────────────────────────────────────────────────
@app.post("/api/projects/{name}/brief")
def api_create_brief(name: str, req: BriefRequest):
    p = get_project_dir(name); ensure_project_dirs(p)
    class FakeArgs: pass
    args = FakeArgs()
    for field_name, field_value in req.model_dump().items(): setattr(args, field_name, field_value)
    brief = target_brief_from_args(args)
    prompt_text = render_generator_prompt(brief)
    brief_path = p / "briefs" / "target_brief.json"
    prompt_path = p / "prompts" / "generator_prompt.md"
    if (brief_path.exists() or prompt_path.exists()) and not req.force:
        raise HTTPException(409, "Target brief 已存在，勾选「覆盖已有」或使用 force")
    write_json(brief_path, brief)
    write_text(prompt_path, prompt_text)
    config_path = p / "config.json"
    cfg = json.loads(config_path.read_text()) if config_path.exists() else default_config()
    target = cfg.get("target", {})
    if not isinstance(target, dict): target = {}
    target.update({"name": req.target_name, "protein_pdb": req.protein_pdb,
        "reference_ligand_sdf": req.reference_ligand,
        "pocket": {"center": parse_triplet(req.center) or [0,0,0],
            "size": parse_triplet(req.size) or [20,20,20],
            "source": req.pocket_source, "description": req.pocket,
            "key_residues": split_items(req.key_residues)}})
    cfg["target"] = target
    cfg["target_brief_file"] = str(brief_path)
    cfg["generator_prompt_file"] = str(prompt_path)
    write_json(p / "config.json", cfg)
    return {"message": "靶点需求已创建", "brief": brief, "prompt": prompt_text}

@app.get("/api/projects/{name}/brief")
def api_get_brief(name: str):
    brief_path = get_project_dir(name) / "briefs" / "target_brief.json"
    if not brief_path.exists(): raise HTTPException(404, "Target brief 不存在，请先创建")
    return json.loads(brief_path.read_text())

@app.get("/api/projects/{name}/prompt")
def api_get_prompt(name: str):
    prompt_path = get_project_dir(name) / "prompts" / "generator_prompt.md"
    if not prompt_path.exists(): raise HTTPException(404, "Generator prompt 不存在，请先创建")
    return {"content": prompt_path.read_text(encoding="utf-8")}

@app.get("/api/prompt-examples")
def api_prompt_examples():
    return [{"title": t, "text": x} for t, x in PROMPT_EXAMPLES]

@app.get("/api/generator-adapters")
def api_generator_adapters():
    return generator_adapters_payload()

@app.get("/api/target-catalog")
def api_target_catalog(query: str = ""):
    try:
        catalog = expanded_target_catalog(query or "")
        targets = sorted_targets(catalog, query or "")
    except SystemExit as exc:
        http_from_cli_exit(exc)
    rows = []
    for target in targets:
        row = catalog_target_summary(target)
        row["score"] = target.get("_target_score", "")
        rows.append(row)
    return {
        "stage": 1,
        "command": "target-catalog",
        "query": query or "",
        "domain": catalog.get("domain", ""),
        "updated": catalog.get("updated", ""),
        "count": len(rows),
        "targets": rows,
        "boundary": "Target catalog entries are evidence-screening presets and must be reviewed before real experimental claims.",
    }

@app.post("/api/projects/{name}/target-intake")
def api_target_intake(name: str, req: TargetIntakeRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    catalog = load_target_catalog(None)
    target, matched_by = select_catalog_target(
        catalog,
        disease=req.disease or "influenza_a",
        target_id=req.target or "",
        target_hint=req.target_hint or "",
        source_text=req.source or "",
        pdb_id=extract_pdb_id_from_text(req.source or ""),
    )
    if not target:
        raise HTTPException(404, "未能从目标目录中解析出可用靶点")
    intake_result = stage1_persist_catalog_brief(
        p,
        target,
        req.disease or "influenza_a",
        req.source_kind or "text",
        req.source or "",
        req.target_hint or "",
        matched_by,
        force=True,
    )
    target_pack = stage1_target_pack_data(
        p,
        max(1, int(req.round)),
        target_id_override=str(target.get("id", "")),
        pdb_id_override=str(target.get("representative_structures", [{}])[0].get("pdb_id", "") if isinstance(target.get("representative_structures"), list) and target.get("representative_structures") else ""),
        reference_ligand_override=str(target.get("binding_site", {}).get("reference_ligand", "")) if isinstance(target.get("binding_site"), dict) else "",
    )
    return {
        "stage": 1,
        "command": "target-intake",
        "project": name,
        "source_kind": req.source_kind,
        "source": req.source,
        "matched_by": matched_by,
        "normalized_target": intake_result["normalized_target"],
        "brief": intake_result["brief"],
        "prompt": intake_result["prompt"],
        "target_pack": target_pack,
        "files": {
            "target_brief": str(intake_result["brief_path"]),
            "generator_prompt": str(intake_result["prompt_path"]),
            "config": str(intake_result["config_path"]),
            "target_pack_json": str(p / "target_pack.json"),
            "target_pack_report": str(p / "reports" / "target_pack_report.md"),
        },
    }

@app.get("/api/projects/{name}/target-pack")
def api_get_target_pack(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    path = p / "target_pack.json"
    return {
        "stage": 1,
        "project": name,
        "round": max(1, int(round_no)),
        "has_target_pack": path.exists(),
        "target_pack": read_json_if_exists(path),
        "files": {
            "target_pack_json": str(path),
            "target_pack_report": str(p / "reports" / "target_pack_report.md"),
        },
    }

@app.post("/api/projects/{name}/target-pack")
def api_build_target_pack(name: str, req: TargetPackRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    target_pack = stage1_target_pack_data(
        p,
        max(1, int(req.round)),
        target_id_override=req.target or "",
        pdb_id_override=req.pdb_id or "",
        reference_ligand_override=req.reference_ligand or "",
    )
    status = "ready" if target_pack.get("pocket", {}).get("center") else "ready_with_warnings"
    return {
        "stage": 1,
        "command": "target-pack",
        "project": name,
        "round": max(1, int(req.round)),
        "status": status,
        "target_pack": target_pack,
        "files": target_pack.get("files", {}),
    }

@app.post("/api/projects/{name}/target-catalog/custom")
def api_custom_target_catalog(name: str, req: TargetCatalogCustomRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    target = custom_target_from_request(req)
    if not target.get("id"):
        raise HTTPException(400, "target_id 不能为空")
    if not target.get("display_name"):
        raise HTTPException(400, "display_name 不能为空")

    catalog_path = p / "targets" / "custom_target_catalog.json"
    existing = read_json_if_exists(catalog_path)
    existing_targets = existing.get("targets", []) if isinstance(existing.get("targets"), list) else []
    merged_targets = [row for row in existing_targets if isinstance(row, dict) and row.get("id") != target["id"]]
    merged_targets.append(target)
    catalog = {
        "schema_version": existing.get("schema_version", "0.1"),
        "domain": existing.get("domain", "project_custom_targets"),
        "updated": datetime.now().date().isoformat(),
        "targets": merged_targets,
    }
    write_json(catalog_path, catalog)
    summary_path = p / "targets" / "custom_target_catalog.csv"
    write_csv(
        summary_path,
        [catalog_target_summary(row) for row in merged_targets],
        ["id", "display_name", "short_name", "recommendation", "target_family", "pathogen_scope", "known_drugs", "primary_pdb", "reference_ligand", "evidence_score"],
    )
    return {
        "stage": 1,
        "command": "custom-target-catalog",
        "project": name,
        "target": target,
        "count": len(merged_targets),
        "files": {
            "custom_catalog": str(catalog_path),
            "custom_catalog_summary": str(summary_path),
        },
        "boundary": "Custom targets are project-level hypotheses and require evidence review before docking or experimental claims.",
    }

@app.post("/api/projects/{name}/stage3/candidates")
def api_stage3_candidates(name: str, req: Stage3CandidateRequest):
    p = get_project_dir(name)
    return stage3_candidates_payload(name, p, req)

@app.get("/api/projects/{name}/stage3")
def api_stage3_status(name: str, round: int = 1):
    p = get_project_dir(name)
    return stage3_status_payload(name, p, round)

def catalog_target_map(catalog: Dict[str, object]) -> Dict[str, Dict[str, object]]:
    targets = catalog.get("targets", [])
    if not isinstance(targets, list):
        return {}
    return {
        str(target.get("id", "")): target
        for target in targets
        if isinstance(target, dict) and target.get("id")
    }

def catalog_target_for_pdb_id(catalog: Dict[str, object], pdb_id: str) -> Dict[str, object]:
    needle = str(pdb_id or "").strip().upper()
    if not needle:
        return {}
    targets = catalog.get("targets", [])
    if not isinstance(targets, list):
        return {}
    for target in targets:
        if not isinstance(target, dict):
            continue
        structures = target.get("representative_structures", [])
        if not isinstance(structures, list):
            continue
        for structure in structures:
            if not isinstance(structure, dict):
                continue
            if str(structure.get("pdb_id", "")).strip().upper() == needle:
                return target
    return {}

def extract_pdb_id_from_text(text: str) -> str:
    for match in re.findall(r"\b[0-9][A-Za-z0-9]{3}\b", str(text or "")):
        return match.upper()
    return ""

def select_catalog_target(
    catalog: Dict[str, object],
    disease: str = "",
    target_id: str = "",
    target_hint: str = "",
    source_text: str = "",
    pdb_id: str = "",
) -> tuple[Dict[str, object], str]:
    by_id = catalog_target_map(catalog)
    requested = str(target_id or "").strip()
    if requested and requested in by_id:
        return by_id[requested], "explicit_target"
    requested_pdb = str(pdb_id or "").strip().upper() or extract_pdb_id_from_text(source_text)
    if requested_pdb:
        pdb_target = catalog_target_for_pdb_id(catalog, requested_pdb)
        if pdb_target:
            return pdb_target, f"pdb_id:{requested_pdb}"
    hint = str(target_hint or "").strip()
    if hint:
        ranked = sorted_targets(catalog, hint)
        if ranked:
            return ranked[0], f"target_hint:{hint}"
    source_query = str(source_text or "").strip()
    if source_query:
        ranked = sorted_targets(catalog, source_query)
        if ranked:
            return ranked[0], "source_text"
    disease_query = str(disease or "").strip()
    if disease_query:
        ranked = sorted_targets(catalog, disease_query)
        if ranked:
            return ranked[0], f"disease:{disease_query}"
    targets = catalog.get("targets", [])
    if isinstance(targets, list):
        for target in targets:
            if isinstance(target, dict):
                return target, "catalog_first"
    return {}, "missing"

def stage1_normalized_target_summary(
    target: Dict[str, object],
    matched_by: str,
    source_kind: str,
    source: str,
    target_hint: str,
    disease: str,
) -> Dict[str, object]:
    structures = target.get("representative_structures", [])
    first_structure = structures[0] if isinstance(structures, list) and structures and isinstance(structures[0], dict) else {}
    site = target.get("binding_site", {}) if isinstance(target.get("binding_site"), dict) else {}
    return {
        "target_id": target.get("id", ""),
        "target_name": target.get("display_name", ""),
        "target_family": target.get("target_family", ""),
        "disease_context": disease,
        "source_kind": source_kind,
        "source": source,
        "target_hint": target_hint,
        "matched_by": matched_by,
        "pdb_id": first_structure.get("pdb_id", ""),
        "reference_ligand": site.get("reference_ligand", first_structure.get("ligand", "")),
        "positive_controls": target.get("known_drugs", []),
    }

def stage1_persist_catalog_brief(
    project_dir: Path,
    target: Dict[str, object],
    disease: str,
    source_kind: str,
    source: str,
    target_hint: str,
    matched_by: str,
    force: bool = True,
) -> Dict[str, Any]:
    brief_path = project_dir / "briefs" / "target_brief.json"
    prompt_path = project_dir / "prompts" / "generator_prompt.md"
    config_path = project_dir / "config.json"
    brief = target_to_brief(
        target,
        disease or "influenza_a",
        source or "",
        55,
        550.0,
    )
    brief["intake_source_kind"] = source_kind
    brief["intake_source"] = source
    brief["target_hint"] = target_hint
    brief["matched_by"] = matched_by
    if (brief_path.exists() or prompt_path.exists()) and not force:
        raise HTTPException(409, "Target brief or prompt already exists. Use force=true to overwrite.")
    write_json(brief_path, brief)
    prompt_text = render_generator_prompt(brief)
    write_text(prompt_path, prompt_text)
    config = read_json(config_path) if config_path.exists() else default_config()
    target_block = config.get("target", {})
    if not isinstance(target_block, dict):
        target_block = {}
    site = brief.get("binding_site", {})
    protein = brief.get("protein", {})
    if not isinstance(site, dict):
        site = {}
    if not isinstance(protein, dict):
        protein = {}
    target_block.update(
        {
            "name": brief.get("target_name", ""),
            "target_catalog_id": brief.get("target_catalog_id", ""),
            "protein_pdb": "",
            "reference_ligand_sdf": site.get("reference_ligand", ""),
            "pdb_id": protein.get("pdb_id", ""),
            "source_kind": source_kind,
            "source": source,
            "target_hint": target_hint,
            "matched_by": matched_by,
            "pocket": {
                "center": site.get("center", []),
                "size": site.get("size", []),
                "source": site.get("source", ""),
                "description": site.get("description", ""),
                "key_residues": site.get("key_residues", []),
                "box_strategy": site.get("box_strategy", ""),
            },
        }
    )
    config["target"] = target_block
    config["target_brief_file"] = str(brief_path)
    config["generator_prompt_file"] = str(prompt_path)
    write_json(config_path, config)
    normalized_target = stage1_normalized_target_summary(target, matched_by, source_kind, source, target_hint, disease)
    return {
        "brief_path": brief_path,
        "prompt_path": prompt_path,
        "config_path": config_path,
        "brief": brief,
        "prompt": prompt_text,
        "normalized_target": normalized_target,
        "config": config,
    }

def read_stage2_target_row(project_dir: Path, target_id: str = "") -> Dict[str, str]:
    rows = csv_to_dicts(project_dir / "evidence" / "stage2_target_sources.csv")
    if target_id:
        for row in rows:
            if str(row.get("target_id", "")) == str(target_id):
                return row
    return rows[0] if rows else {}

def stage1_target_pack_data(
    project_dir: Path,
    round_no: int,
    target_id_override: str = "",
    pdb_id_override: str = "",
    reference_ligand_override: str = "",
) -> Dict[str, Any]:
    brief_path = project_dir / "briefs" / "target_brief.json"
    config_path = project_dir / "config.json"
    target_pack_path = project_dir / "target_pack.json"
    report_path = project_dir / "reports" / "target_pack_report.md"
    brief = read_json_if_exists(brief_path)
    config = read_json_if_exists(config_path)
    target_block = config.get("target", {}) if isinstance(config.get("target"), dict) else {}
    catalog = load_target_catalog(None)
    target_id = (
        str(target_id_override or "").strip()
        or str(brief.get("target_catalog_id", "") or target_block.get("target_catalog_id", "") or "").strip()
    )
    target = {}
    if target_id:
        try:
            target = get_target_by_id(catalog, target_id)
        except SystemExit:
            target = {}
    if not target:
        target = select_catalog_target(
            catalog,
            str(brief.get("disease_context", "") or target_block.get("source_kind", "") or "influenza_a"),
            target_id=target_id,
            target_hint=str(target_block.get("target_hint", "") or brief.get("target_hint", "") or ""),
            source_text=str(brief.get("intake_source", "") or target_block.get("source", "") or ""),
            pdb_id=str(pdb_id_override or target_block.get("pdb_id", "") or ""),
        )[0]
    if not target:
        target = select_catalog_target(catalog, "influenza_a")[0]
    stage2_row = read_stage2_target_row(project_dir, target.get("id", ""))
    stage4_pack_path = project_dir / "stage4" / "pocket_pack.json"
    stage4_pack = read_json_if_exists(stage4_pack_path)
    site = brief.get("binding_site", {}) if isinstance(brief.get("binding_site"), dict) else {}
    target_site = target.get("binding_site", {}) if isinstance(target.get("binding_site"), dict) else {}
    positive_controls = split_items(stage2_row.get("positive_controls", ""))
    if not positive_controls:
        known_drugs = target.get("known_drugs", [])
        positive_controls = [str(item) for item in known_drugs] if isinstance(known_drugs, list) else split_items(str(known_drugs))
    reference_controls = split_items(stage2_row.get("reference_controls", ""))
    historical_controls = split_items(stage2_row.get("historical_controls", ""))
    pocket = stage4_pack.get("pocket", {}) if isinstance(stage4_pack.get("pocket"), dict) else {}
    if not pocket:
        pocket = {
            "status": "ready" if (site.get("center") or site.get("size") or target_block.get("pocket")) else "needs_review",
            "source": site.get("source", "") or target_site.get("source", "") or target_block.get("pocket", {}).get("source", ""),
            "center": site.get("center", []) or target_block.get("pocket", {}).get("center", []),
            "size": site.get("size", []) or target_block.get("pocket", {}).get("size", []),
            "description": site.get("description", "") or target_site.get("description", ""),
            "reference_ligand": reference_ligand_override or site.get("reference_ligand", "") or target_site.get("reference_ligand", ""),
            "key_residues": site.get("key_residues", []) or target_site.get("key_residues", []),
            "box_strategy": site.get("box_strategy", "") or target_site.get("box_strategy", ""),
        }
    evidence_score = stage2_row.get("evidence_score") or target.get("scores", {}).get("clinical_validation", "") or target.get("scores", {}).get("disease_relevance", "")
    confidence = evidence_score or target.get("scores", {}).get("pocket_confidence", "")
    target_pack = {
        "schema_version": "0.1",
        "stage": 1,
        "project": project_dir.name,
        "round": round_no,
        "generated_at": datetime.now().isoformat(),
        "target_id": target.get("id", ""),
        "target_name": target.get("display_name", ""),
        "disease_context": brief.get("disease_context", "") or target_block.get("disease_context", ""),
        "source_kind": brief.get("intake_source_kind", "") or target_block.get("source_kind", ""),
        "source": brief.get("intake_source", "") or target_block.get("source", ""),
        "target_hint": brief.get("target_hint", "") or target_block.get("target_hint", ""),
        "matched_by": brief.get("matched_by", "") or target_block.get("matched_by", ""),
        "pdb_id": str(pdb_id_override or brief.get("protein", {}).get("pdb_id", "") or target_block.get("pdb_id", "") or stage2_row.get("primary_pdb", "") or ""),
        "reference_ligand": reference_ligand_override or site.get("reference_ligand", "") or target_site.get("reference_ligand", "") or stage2_row.get("reference_ligand", ""),
        "evidence_score": evidence_score,
        "confidence": confidence,
        "controls": {
            "positive": positive_controls,
            "reference": reference_controls,
            "historical": historical_controls,
        },
        "pocket": pocket,
        "evidence_matrix": {
            "stage2_row": stage2_row,
            "stage2_assets": read_json_if_exists(project_dir / "evidence" / "stage2_closed_loop_assets.json"),
            "stage2_matrix": csv_to_dicts(project_dir / "evidence" / "stage2_target_sources.csv"),
        },
        "brief": brief,
        "boundary": [
            "Target pack combines target brief, evidence, and pocket metadata for computational screening only.",
            "It does not prove biological activity, therapeutic efficacy, safety, or clinical usefulness.",
            "Known drugs are controls and references, not new project hits.",
        ],
        "files": {
            "target_pack_json": str(target_pack_path),
            "target_pack_report": str(report_path),
            "target_brief": str(brief_path),
            "config": str(config_path),
        },
    }
    write_json(target_pack_path, target_pack)
    report_lines = [
        f"# Target Pack - {target.get('display_name', target_pack['target_id'])}",
        "",
        f"- Target ID: `{target_pack['target_id']}`",
        f"- Disease: `{target_pack['disease_context']}`",
        f"- PDB: `{target_pack['pdb_id']}`",
        f"- Reference ligand: `{target_pack['reference_ligand']}`",
        f"- Evidence score: `{target_pack['evidence_score']}`",
        f"- Confidence: `{target_pack['confidence']}`",
        f"- Pocket source: `{pocket.get('source', '')}`",
        f"- Pocket center: `{pocket.get('center', [])}`",
        f"- Pocket size: `{pocket.get('size', [])}`",
        "",
        "## Controls",
        "",
        f"- Positive: {', '.join(positive_controls) if positive_controls else '-'}",
        f"- Reference: {', '.join(reference_controls) if reference_controls else '-'}",
        f"- Historical: {', '.join(historical_controls) if historical_controls else '-'}",
        "",
        "## Boundary",
        "",
        "- Computational target pack only.",
        "- Not efficacy proof.",
    ]
    write_text(report_path, "\n".join(report_lines) + "\n")
    return target_pack

def stage4_build_pocket_pack(
    project_dir: Path,
    round_no: int,
    target_id: str = "",
    pdb_id: str = "",
    reference_ligand: str = "",
) -> Dict[str, Any]:
    config_path = project_dir / "config.json"
    config = read_json_if_exists(config_path)
    target_block = config.get("target", {}) if isinstance(config.get("target"), dict) else {}
    brief = read_json_if_exists(project_dir / "briefs" / "target_brief.json")
    brief_site = brief.get("binding_site", {}) if isinstance(brief.get("binding_site"), dict) else {}
    target_id = str(target_id or brief.get("target_catalog_id", "") or target_block.get("target_catalog_id", "") or "").strip()
    if not target_id:
        catalog = load_target_catalog(None)
        target, _matched = select_catalog_target(catalog, "influenza_a")
        target_id = str(target.get("id", ""))
    requested_pdb = str(pdb_id or brief.get("protein", {}).get("pdb_id", "") or target_block.get("pdb_id", "") or "").strip().upper()
    if not requested_pdb and target_id:
        catalog = load_target_catalog(None)
        target = {}
        try:
            target = get_target_by_id(catalog, target_id)
        except SystemExit:
            target = {}
        if isinstance(target, dict):
            structures = target.get("representative_structures", [])
            if isinstance(structures, list) and structures and isinstance(structures[0], dict):
                requested_pdb = str(structures[0].get("pdb_id", "")).upper()
    receptor_pdb = ""
    receptor_pdbqt = ""
    if requested_pdb:
        receptor_pdb, receptor_pdbqt = find_stage4_project_receptor(project_dir, requested_pdb)
    receptor_path = Path(receptor_pdb) if receptor_pdb else None
    pocket_data: Dict[str, Any] = {}
    pocket_source = ""
    if receptor_path and receptor_path.exists():
        extracted = extract_stage4_cocrystal_pocket(receptor_path, reference_ligand or str(brief_site.get("reference_ligand", "")))
        if extracted.get("status") == "extracted_from_receptor":
            pocket_source = "co_crystal"
            pocket_data = {
                "status": "ready",
                "source": pocket_source,
                "reference_ligand": reference_ligand or brief_site.get("reference_ligand", ""),
                "center": extracted.get("center", []),
                "size": extracted.get("size", []),
                "description": extracted.get("message", ""),
                "candidate_ligands": extracted.get("candidate_ligands", []),
                "detected_ligand": extracted.get("detected_ligand", {}),
                "box_padding_angstrom": extracted.get("box_padding_angstrom", ""),
                "box_min_size_angstrom": extracted.get("box_min_size_angstrom", ""),
            }
        else:
            pocket_data = {
                "status": "needs_review",
                "source": brief_site.get("source", "") or target_block.get("pocket", {}).get("source", "") or "curated",
                "reference_ligand": reference_ligand or brief_site.get("reference_ligand", ""),
                "center": brief_site.get("center", []) or target_block.get("pocket", {}).get("center", []),
                "size": brief_site.get("size", []) or target_block.get("pocket", {}).get("size", []),
                "description": extracted.get("message", "") or brief_site.get("description", ""),
                "candidate_ligands": extracted.get("candidate_ligands", []),
            }
    else:
        pocket_data = {
            "status": "needs_review" if (brief_site.get("center") or brief_site.get("size")) else "manual",
            "source": brief_site.get("source", "") or target_block.get("pocket", {}).get("source", "") or "manual",
            "reference_ligand": reference_ligand or brief_site.get("reference_ligand", ""),
            "center": brief_site.get("center", []) or target_block.get("pocket", {}).get("center", []),
            "size": brief_site.get("size", []) or target_block.get("pocket", {}).get("size", []),
            "description": brief_site.get("description", "") or target_block.get("pocket", {}).get("description", ""),
        }
    if not pocket_data.get("center"):
        pocket_data["center"] = []
    if not pocket_data.get("size"):
        pocket_data["size"] = []
    if pocket_data.get("center") and pocket_data.get("size") and pocket_data.get("status") in {"manual", "needs_review"}:
        pocket_data["status"] = "ready" if pocket_data.get("source") in {"co_crystal", "curated"} else "needs_review"
    pocket_pack_path = project_dir / "stage4" / "pocket_pack.json"
    pocket_pack = {
        "schema_version": "0.1",
        "stage": 4,
        "project": project_dir.name,
        "round": round_no,
        "generated_at": datetime.now().isoformat(),
        "target_id": target_id,
        "pdb_id": requested_pdb,
        "receptor": {
            "local_receptor_pdb": receptor_pdb,
            "local_receptor_pdbqt": receptor_pdbqt,
        },
        "pocket": pocket_data,
        "source": {
            "brief_file": str(project_dir / "briefs" / "target_brief.json"),
            "config_file": str(config_path),
        },
        "boundary": [
            "Pocket pack records computational docking coordinates only.",
            "It is not wet-lab validation or efficacy proof.",
        ],
        "files": {
            "pocket_pack": str(pocket_pack_path),
            "brief": str(project_dir / "briefs" / "target_brief.json"),
            "config": str(config_path),
        },
    }
    write_json(pocket_pack_path, pocket_pack)
    if requested_pdb or pocket_data.get("center") or pocket_data.get("size"):
        target_block = config.get("target", {}) if isinstance(config.get("target"), dict) else {}
        if not isinstance(target_block, dict):
            target_block = {}
        target_block["pocket"] = {
            "center": pocket_data.get("center", []),
            "size": pocket_data.get("size", []),
            "source": pocket_data.get("source", ""),
            "description": pocket_data.get("description", ""),
            "reference_ligand": pocket_data.get("reference_ligand", ""),
        }
        target_block["pdb_id"] = requested_pdb or target_block.get("pdb_id", "")
        if target_id:
            target_block["target_catalog_id"] = target_id
        config["target"] = target_block
        write_json(config_path, config)
    return pocket_pack

def stage8_export_bundle(project_dir: Path, round_no: int) -> Dict[str, Any]:
    target_pack = stage1_target_pack_data(project_dir, round_no)
    pocket_pack = stage4_build_pocket_pack(
        project_dir,
        round_no,
        target_id=str(target_pack.get("target_id", "") or ""),
        pdb_id=str(target_pack.get("pdb_id", "") or ""),
        reference_ligand=str(target_pack.get("reference_ligand", "") or ""),
    )
    stage7_path = project_dir / "stage7" / f"round_{round_no}_delivery_manifest.json"
    if not stage7_path.exists():
        stage7_package(argparse.Namespace(project=str(project_dir), round=round_no, title=f"{project_dir.name} Stage 7 Delivery Package"))
    stage5_path = project_dir / "stage5" / "dashboard_data.json"
    if not stage5_path.exists():
        stage5_dashboard(argparse.Namespace(project=str(project_dir), round=round_no, title=f"{project_dir.name} Stage 5 Dashboard"))
    stage6_path = project_dir / "stage6" / f"round_{round_no}_validation_assets.json"
    if not stage6_path.exists():
        stage6_validate(argparse.Namespace(project=str(project_dir), round=round_no, top=8))

    export_dir = project_dir / "exports"
    export_dir.mkdir(parents=True, exist_ok=True)
    zip_path = export_dir / f"round_{round_no}_export_bundle.zip"
    manifest_path = export_dir / f"round_{round_no}_export_manifest.json"

    include_map = {
        "target_pack_json": project_dir / "target_pack.json",
        "target_pack_report": project_dir / "reports" / "target_pack_report.md",
        "target_brief": project_dir / "briefs" / "target_brief.json",
        "generator_prompt": project_dir / "prompts" / "generator_prompt.md",
        "stage2_matrix": project_dir / "evidence" / "stage2_target_sources.csv",
        "stage2_assets": project_dir / "evidence" / "stage2_closed_loop_assets.json",
        "stage4_pocket_pack": project_dir / "stage4" / "pocket_pack.json",
        "stage4_receptor_package": project_dir / "stage4" / f"round_{round_no}_receptor_package.json",
        "stage4_assets": project_dir / "stage4" / f"round_{round_no}_stage4_assets.json",
        "stage4_report": project_dir / "reports" / f"stage4_round_{round_no}_report.md",
        "stage5_dashboard_data": project_dir / "stage5" / "dashboard_data.json",
        "stage5_dashboard_html": project_dir / "stage5" / "index.html",
        "stage6_assets": stage6_path,
        "stage6_quality_gates": project_dir / "stage6" / f"round_{round_no}_quality_gates.csv",
        "stage6_hit_triage": project_dir / "stage6" / f"round_{round_no}_hit_triage.csv",
        "stage6_assay_queue": project_dir / "stage6" / f"round_{round_no}_assay_queue.csv",
        "stage7_manifest": stage7_path,
        "stage7_summary": project_dir / "stage7" / f"round_{round_no}_executive_summary.md",
        "stage7_reproducibility": project_dir / "stage7" / f"round_{round_no}_reproducibility.md",
        "stage7_checklist": project_dir / "stage7" / f"round_{round_no}_investor_demo_checklist.csv",
        "stage7_spec": project_dir / "stage7" / "stage8_frontend_product_spec.md",
        "stage7_report": project_dir / "reports" / f"stage7_round_{round_no}_delivery_report.md",
        "stage8_command_center": None,
    }
    stage8_payload_data = stage8_payload(project_dir.name, project_dir, round_no)
    write_json(manifest_path, {
        "schema_version": "0.1",
        "stage": 8,
        "project": project_dir.name,
        "round": round_no,
        "generated_at": datetime.now().isoformat(),
        "included_files": [key for key, path in include_map.items() if path and Path(path).exists()],
        "boundary": stage8_payload_data.get("boundary", []),
        "package_status": "ready",
    })
    include_map["stage8_command_center"] = manifest_path

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for key, path in include_map.items():
            if not path:
                continue
            file_path = Path(path)
            if not file_path.exists() or not file_path.is_file():
                continue
            zf.write(file_path, arcname=file_path.relative_to(project_dir).as_posix())
        zf.writestr(
            "export_manifest.txt",
            "\n".join(
                [
                    f"Project: {project_dir.name}",
                    f"Round: {round_no}",
                    f"Target pack: {target_pack.get('target_id', '')}",
                    f"Pocket source: {pocket_pack.get('pocket', {}).get('source', '')}",
                    "Included files:",
                    *[f"- {key}" for key, path in include_map.items() if path and Path(path).exists()],
                ]
            )
            + "\n",
        )
        zf.writestr("stage8_command_center.json", json.dumps(stage8_payload_data, ensure_ascii=False, indent=2))

    included_files = [key for key, path in include_map.items() if path and Path(path).exists()]
    return {
        "stage": 8,
        "project": project_dir.name,
        "round": round_no,
        "package_status": "ready" if {"target_pack_json", "stage7_manifest"}.issubset(set(included_files)) else "ready_with_warnings",
        "included_files": included_files,
        "files": {
            "zip": str(zip_path),
            "manifest": str(manifest_path),
        },
        "target_pack": target_pack,
        "pocket_pack": pocket_pack,
        "stage8": stage8_payload_data,
    }

def product_check(check_id: str, label: str, status: str, evidence: str, remedy: str = "") -> Dict[str, Any]:
    return {
        "check_id": check_id,
        "label": label,
        "status": status,
        "evidence": evidence,
        "remedy": remedy,
    }

def product_section(label: str, checks: List[Dict[str, Any]]) -> Dict[str, Any]:
    if any(item.get("status") == "missing" for item in checks):
        status = "missing"
    elif any(item.get("status") == "warn" for item in checks):
        status = "warn"
    else:
        status = "ready"
    return {"label": label, "status": status, "checks": checks}

def product_system_health_payload(default_project: str = "") -> Dict[str, Any]:
    projects = list_projects()
    names = [str(p.get("name", "")) for p in projects]
    selected_project = default_project if default_project in names else (("flu_na_real_demo" if "flu_na_real_demo" in names else "") or (names[0] if names else ""))
    capabilities = stage4_capabilities()
    modules = capabilities.get("modules", {}) if isinstance(capabilities.get("modules"), dict) else {}
    executables = capabilities.get("executables", {}) if isinstance(capabilities.get("executables"), dict) else {}
    optional = capabilities.get("optional_libraries", {}) if isinstance(capabilities.get("optional_libraries"), dict) else {}
    real = capabilities.get("real_libraries", {}) if isinstance(capabilities.get("real_libraries"), dict) else {}
    docking_backend = capabilities.get("docking_backend", {}) if isinstance(capabilities.get("docking_backend"), dict) else {}

    def module_item(name: str, *legacy_groups: Dict[str, Any]) -> Dict[str, Any]:
        item = modules.get(name, {}) if isinstance(modules.get(name, {}), dict) else {}
        if item:
            return item
        for group in legacy_groups:
            legacy_item = group.get(name, {}) if isinstance(group.get(name, {}), dict) else {}
            if legacy_item:
                return legacy_item
        if name == "rdkit" and isinstance(capabilities.get("rdkit"), dict):
            return capabilities["rdkit"]
        return {}

    def exe_item(name: str) -> Dict[str, Any]:
        item = executables.get(name, {}) if isinstance(executables.get(name, {}), dict) else {}
        return item

    def module_ready(name: str, *legacy_groups: Dict[str, Any]) -> bool:
        return module_item(name, *legacy_groups).get("status") == "available"

    def exe_ready(name: str) -> bool:
        return exe_item(name).get("status") == "found"

    def module_evidence(name: str, *legacy_groups: Dict[str, Any]) -> str:
        item = module_item(name, *legacy_groups)
        return str(item.get("version") or item.get("origin") or item.get("status") or "not available")

    def exe_evidence(name: str) -> str:
        item = exe_item(name)
        return str(item.get("path") or item.get("source") or "not found")

    def exe_status(name: str) -> Dict[str, Any]:
        found = exe_ready(name)
        return product_check(
            name,
            name,
            "ready" if found else "warn",
            exe_evidence(name),
            f"Install or expose {name} on PATH, or continue with proxy/external-score mode." if not found else "",
        )

    def gnina_status() -> Dict[str, Any]:
        if exe_ready("gnina") or module_ready("gnina"):
            return product_check("gnina", "gnina", "ready", exe_evidence("gnina") if exe_ready("gnina") else module_evidence("gnina"), "")
        backends = docking_backend.get("available_backends", []) if isinstance(docking_backend.get("available_backends"), list) else []
        if "vina" in backends or exe_ready("vina") or module_ready("vina"):
            return product_check(
                "gnina",
                "gnina",
                "ready",
                "optional backend not installed; Vina backend is available",
                "Use Docker/Linux GNINA only when CNN rescoring is required.",
            )
        return product_check(
            "gnina",
            "gnina",
            "warn",
            "not found",
            "Install or expose gnina on PATH, or continue with proxy/external-score mode.",
        )

    def docking_backend_status() -> Dict[str, Any]:
        backends = docking_backend.get("available_backends", []) if isinstance(docking_backend.get("available_backends"), list) else []
        ready = docking_backend.get("status") == "available" or exe_ready("vina") or exe_ready("gnina") or module_ready("vina") or module_ready("gnina")
        evidence = ", ".join(str(item) for item in backends) if backends else str(docking_backend.get("message") or "no backend")
        if not evidence or evidence == "no backend":
            evidence = "vina/gnina not found"
        return product_check(
            "docking_backend",
            "Docking backend",
            "ready" if ready else "warn",
            evidence,
            "Install Vina or GNINA, or import an external true-score CSV." if not ready else "",
        )

    py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    runtime_checks = [
        product_check("python", "Python runtime", "ready", py_version, ""),
        product_check("projects_root", "Projects root", "ready" if PROJECTS_ROOT.exists() else "missing", str(PROJECTS_ROOT), "Create the webapp/projects directory."),
        product_check("project_count", "Registered projects", "ready" if projects else "warn", str(len(projects)), "Create or import a demo project such as flu_na_real_demo."),
    ]
    frontend_checks = [
        product_check("index_html", "Web entry", "ready" if (STATIC / "index.html").exists() else "missing", str(STATIC / "index.html"), "Restore webapp/static/index.html."),
        product_check("local_3dmol", "3Dmol.js local vendor", "ready" if (STATIC / "vendor" / "3Dmol-min.js").exists() else "missing", str(STATIC / "vendor" / "3Dmol-min.js"), "Restore the local 3Dmol vendor file so 3D structures work offline."),
    ]
    computational_checks = [
        product_check("rdkit", "RDKit", "ready" if module_ready("rdkit", real) else "warn", module_evidence("rdkit", real), "Install RDKit for real descriptors; proxy-only mode remains available." if not module_ready("rdkit", real) else ""),
        product_check("openbabel", "OpenBabel", "ready" if module_ready("openbabel", optional) or exe_ready("obabel") else "warn", module_evidence("openbabel", optional) if module_ready("openbabel", optional) else exe_evidence("obabel"), "Install OpenBabel for ligand/receptor preparation." if not (module_ready("openbabel", optional) or exe_ready("obabel")) else ""),
        product_check("meeko", "Meeko ligand prep", "ready" if module_ready("meeko", optional) or exe_ready("mk_prepare_ligand.py") else "warn", module_evidence("meeko", optional) if module_ready("meeko", optional) else exe_evidence("mk_prepare_ligand.py"), "Install Meeko for PDBQT ligand/receptor preparation." if not (module_ready("meeko", optional) or exe_ready("mk_prepare_ligand.py")) else ""),
        product_check("posebusters", "PoseBusters", "ready" if module_ready("posebusters", optional) or exe_ready("bust") else "warn", module_evidence("posebusters", optional) if module_ready("posebusters", optional) else exe_evidence("bust"), "Install PoseBusters for pose quality checks." if not (module_ready("posebusters", optional) or exe_ready("bust")) else ""),
        exe_status("vina"),
        docking_backend_status(),
        gnina_status(),
    ]
    delivery_checks = [
        product_check("default_project", "Default/demo project", "ready" if selected_project else "warn", selected_project or "none", "Create or restore flu_na_real_demo."),
        product_check("start_script", "One-command startup", "ready" if (ROOT / "start_web.sh").exists() else "missing", str(ROOT / "start_web.sh"), "Restore start_web.sh."),
        product_check("stage8_api", "Stage 8 API surface", "ready", "/api/projects/{name}/stage8/review-mode", ""),
    ]
    sections = {
        "runtime": product_section("Runtime and project storage", runtime_checks),
        "frontend": product_section("Frontend and visualization", frontend_checks),
        "computational_tools": product_section("Computational libraries and docking tools", computational_checks),
        "delivery": product_section("Demo and delivery readiness", delivery_checks),
    }
    statuses = [section["status"] for section in sections.values()]
    if any(status == "missing" for status in statuses[:2]):
        overall = "needs_setup"
    elif any(status in {"missing", "warn"} for status in statuses):
        overall = "ready_with_warnings"
    else:
        overall = "ready"
    return {
        "command": "system-health",
        "overall_status": overall,
        "generated_at": datetime.now().isoformat(),
        "service": "ai-molecule-screening-web",
        "url": "http://localhost:8765/",
        "projects_count": len(projects),
        "default_project": selected_project,
        "sections": sections,
        "recommended_commands": [
            "python3 webapp/server.py",
            "python3 health_check.py",
            "open http://localhost:8765/",
        ],
        "boundary": [
            "Computational screening environment readiness only.",
            "Tool availability does not validate potency, efficacy, toxicity, safety, dosing, or clinical usefulness.",
        ],
    }

def product_bootstrap_health_payload(project: str = "", round_no: int = 1) -> Dict[str, Any]:
    system_health = product_system_health_payload(default_project=project)
    selected = project or system_health.get("default_project", "")
    project_dir = PROJECTS_ROOT / selected if selected else PROJECTS_ROOT
    preflight = stage8_preflight_payload(selected, project_dir, round_no) if selected and project_dir.exists() else {}
    guide = stage8_demo_guide_payload(selected, project_dir, round_no) if selected and project_dir.exists() else {}
    stage8_ready = preflight.get("status") == "ready"
    startup_status = system_health.get("overall_status", "needs_setup")
    if startup_status == "ready" and not stage8_ready:
        startup_status = "ready_with_warnings"
    return {
        "stage": "product",
        "command": "product-bootstrap-health",
        "project": selected,
        "round": max(1, int(round_no)),
        "startup_status": startup_status,
        "system_health": system_health,
        "sections": {
            "startup": {
                "label": "Startup and service",
                "status": system_health.get("overall_status", "needs_setup"),
                "checks": [
                    product_check("start_script", "start_web.sh", "ready" if (ROOT / "start_web.sh").exists() else "missing", str(ROOT / "start_web.sh"), "Restore start_web.sh."),
                    product_check("server_py", "FastAPI server", "ready" if (ROOT / "webapp" / "server.py").exists() else "missing", str(ROOT / "webapp" / "server.py"), "Restore webapp/server.py."),
                    product_check("health_check", "health_check.py", "ready" if (ROOT / "health_check.py").exists() else "missing", str(ROOT / "health_check.py"), "Restore health_check.py."),
                ],
            },
            "demo": {
                "label": "Demo closed-loop readiness",
                "status": "ready" if stage8_ready else "warn",
                "checks": preflight.get("checks", [])[:12] if preflight else [],
                "next_primary_action": guide.get("next_primary_action", {}) if guide else {},
            },
            "delivery": {
                "label": "Delivery and handoff",
                "status": "ready" if selected else "warn",
                "checks": [
                    product_check("full_export_api", "Full export API", "ready", f"/api/projects/{selected}/stage8/full-export" if selected else "/api/projects/{name}/stage8/full-export", ""),
                    product_check("standard_delivery_script", "Standard delivery script", "ready" if (ROOT / "scripts" / "build_product_delivery.py").exists() else "warn", str(ROOT / "scripts" / "build_product_delivery.py"), "Run scripts/build_product_delivery.py after implementation."),
                ],
            },
        },
        "quick_start": [
            "./start_web.sh",
            "python3 health_check.py",
            "open http://localhost:8765/",
            f"POST /api/projects/{selected}/stage8/repair" if selected else "Create/select a project, then run Stage 8 repair.",
        ],
        "urls": {
            "app": "http://localhost:8765/",
            "health": "http://localhost:8765/api/health",
            "system_health": "http://localhost:8765/api/system/health",
        },
        "boundary": [
            "Computational startup and demo readiness only.",
            "No efficacy, potency, toxicity, safety, dosing, clinical benefit, or therapeutic claim is created.",
        ],
    }

def stage4_operator_guide_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    stage4 = stage4_payload(name, project_dir, round_no)
    capabilities = stage4_capabilities()
    plan = stage4.get("docking_plan", {}) if isinstance(stage4.get("docking_plan"), dict) else {}
    rec = stage4.get("receptor_package", {}) if isinstance(stage4.get("receptor_package"), dict) else {}
    results = stage4.get("docking_results", []) if isinstance(stage4.get("docking_results"), list) else []
    preflight = stage4.get("stage4_preflight", []) if isinstance(stage4.get("stage4_preflight"), list) else []
    executables = capabilities.get("executables", {}) if isinstance(capabilities.get("executables"), dict) else {}
    modules = capabilities.get("modules", {}) if isinstance(capabilities.get("modules"), dict) else {}

    def tool_status(tool: str, exe_key: str = "") -> Dict[str, str]:
        exe = executables.get(exe_key or tool, {}) if isinstance(executables.get(exe_key or tool, {}), dict) else {}
        mod = modules.get(tool, {}) if isinstance(modules.get(tool, {}), dict) else {}
        status = "available" if exe.get("status") == "found" or mod.get("status") == "available" else "missing"
        evidence = exe.get("path") or mod.get("version") or mod.get("status") or "not found"
        return {"tool": tool, "status": status, "evidence": str(evidence)}

    template_path = plan.get("expected_scores_csv") or str(project_dir / "stage4" / f"round_{round_no}_docking_scores_template.csv")
    workflow_steps = [
        {"step": "1", "label": "确认候选输入", "status": next((x.get("status") for x in preflight if x.get("step_id") == "candidate_input"), "missing"), "action": "检查 candidates/round_N_candidates.csv。"},
        {"step": "2", "label": "确认受体和口袋", "status": "ready" if rec.get("local_receptor_pdb") or rec.get("local_receptor_pdbqt") else "missing", "action": "应用靶点预设、填入 receptor 路径，或先生成 pocket pack。"},
        {"step": "3", "label": "生成 docking-ready 资产", "status": "ready" if stage4.get("has_assets") else "missing", "action": "运行 Stage 4 真实库校验。"},
        {"step": "4", "label": "接入真实分数", "status": "ready" if results else "pending", "action": f"Run Vina/GNINA or import external true-score CSV: {template_path}"},
        {"step": "5", "label": "Pose QC 与解释", "status": "ready" if any(str(r.get("pose_pass", "")).lower() in {"true", "pass", "passed", "yes", "1"} for r in results) else "pending", "action": "Run PoseBusters or import pose_pass/posebusters_report columns."},
        {"step": "6", "label": "同池校准", "status": "ready" if (project_dir / "stage4_5" / f"round_{round_no}_control_docking_scores.csv").exists() else "pending", "action": "Run Stage 4.5 with candidates, controls, and decoys."},
    ]
    return {
        "stage": 4,
        "command": "stage4-operator-guide",
        "project": name,
        "round": round_no,
        "docking_status": plan.get("status", "missing"),
        "receptor_status": rec.get("preparation_status", "missing"),
        "score_count": len(results),
        "score_csv_template": template_path,
        "workflow_steps": workflow_steps,
        "tool_readiness": [
            tool_status("rdkit"),
            tool_status("openbabel", "obabel"),
            tool_status("meeko", "mk_prepare_ligand.py"),
            tool_status("vina"),
            tool_status("gnina"),
            tool_status("posebusters", "bust"),
        ],
        "recommended_commands": (plan.get("commands", []) if isinstance(plan.get("commands"), list) else [])[:8],
        "boundary": [
            "Docking scores are computational ranking evidence only.",
            "External scores must be generated under reviewed receptor, ligand, pocket, and control settings.",
            "Pose QC and controls are required before using docking results for prioritization.",
        ],
    }

def stage8_action_plan_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    preflight = stage8_preflight_payload(name, project_dir, round_no)
    guide = stage8_demo_guide_payload(name, project_dir, round_no)
    command_center = stage8_payload(name, project_dir, round_no)
    missing_or_warn = [row for row in guide.get("steps", []) if row.get("status") != "ready"]
    primary = guide.get("next_primary_action", {}) or (missing_or_warn[0] if missing_or_warn else {})
    if primary:
        primary = dict(primary)
        primary["autofix"] = f"POST /api/projects/{name}/stage8/repair"
    actions = []
    if preflight.get("can_auto_repair"):
        actions.append({"label": "自动补齐缺失产物", "api": f"/api/projects/{name}/stage8/repair", "method": "POST"})
    actions.extend([
        {"label": "运行完整 Demo", "api": f"/api/projects/{name}/stage8/demo-runner", "method": "POST"},
        {"label": "生成闭环报告", "api": f"/api/projects/{name}/stage8/report", "method": "POST"},
        {"label": "导出证据包", "api": f"/api/projects/{name}/stage8/evidence-pack", "method": "POST"},
        {"label": "完整交付导出", "api": f"/api/projects/{name}/stage8/full-export", "method": "POST"},
    ])
    overall = "ready" if preflight.get("status") == "ready" else ("needs_repair" if preflight.get("can_auto_repair") else "ready_with_warnings")
    return {
        "stage": 8,
        "command": "stage8-action-plan",
        "project": name,
        "round": round_no,
        "overall_status": overall,
        "primary_action": primary,
        "missing_or_warn": missing_or_warn,
        "actions": actions,
        "preflight": preflight,
        "completion": guide.get("completion", {}),
        "stage8": command_center,
        "boundary": [
            "Action plan repairs computational workflow artifacts only.",
            "Generated outputs remain virtual-screening evidence and require experimental validation.",
        ],
    }

def generator_adapters_payload() -> Dict[str, Any]:
    def repo_status(rel: str) -> str:
        return "local_repo" if (ROOT / rel).exists() else "planned"

    adapters = [
        {"id": "proxy_smiles", "label": "Built-in proxy SMILES generator", "type": "local", "status": "available", "input": "target brief / prompt", "output": "candidate SMILES CSV", "integration": "ai_mol_loop.generate_candidates"},
        {"id": "openai_prompt", "label": "OpenAI prompt-based candidate intake", "type": "api", "status": "available", "input": "URL/text/prompt + optional API key", "output": "Stage 3 candidates", "integration": "POST /api/projects/{name}/stage3/candidates"},
        {"id": "reinvent4", "label": "REINVENT4", "type": "external_repo", "status": repo_status("ai_molecule_design_repos/REINVENT4"), "input": "scoring config / SMILES seeds", "output": "generated SMILES", "integration": "adapter planned through Stage 3 CSV import"},
        {"id": "drugex", "label": "DrugEx", "type": "external_repo", "status": repo_status("ai_molecule_design_repos/DrugEx"), "input": "training/inference config", "output": "generated SMILES", "integration": "adapter planned through Stage 3 CSV import"},
        {"id": "diffsbdd", "label": "DiffSBDD", "type": "external_repo", "status": "planned", "input": "protein pocket / receptor structure", "output": "structure-based ligand candidates", "integration": "adapter planned after receptor/pocket package stabilization"},
        {"id": "colabfold", "label": "ColabFold", "type": "external_tool", "status": "planned", "input": "protein sequence", "output": "predicted receptor structure", "integration": "adapter planned before Stage 4 receptor package"},
    ]
    return {
        "command": "generator-adapters",
        "adapters": adapters,
        "recommended_order": ["proxy_smiles", "openai_prompt", "reinvent4", "drugex", "diffsbdd", "colabfold"],
        "api_key_policy": "API keys are accepted per request when needed, used in memory, and not persisted to project files.",
        "boundary": [
            "Generator adapters produce candidate hypotheses only.",
            "Generated molecules require RDKit validation, controls, decoys, docking/QC, and experimental follow-up before any claim.",
        ],
    }

def product_deep_integration_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    paths = stage8_paths_for_round(project_dir, safe_round)
    catalog = expanded_target_catalog("")
    catalog_targets = catalog.get("targets", []) if isinstance(catalog.get("targets"), list) else []
    adapters = generator_adapters_payload()
    adapter_rows = adapters.get("adapters", []) if isinstance(adapters.get("adapters"), list) else []
    stage4_guide = stage4_operator_guide_payload(name, project_dir, safe_round)
    action_plan = stage8_action_plan_payload(name, project_dir, safe_round)
    preflight = action_plan.get("preflight", {}) if isinstance(action_plan.get("preflight"), dict) else {}
    science = scientific_readiness_payload(name, project_dir, safe_round)

    def area(
        area_id: str,
        label: str,
        status: str,
        evidence: str,
        route: str,
        api: str,
        actions: List[str],
        artifacts: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Any]:
        safe_status = status if status in {"ready", "warn", "missing"} else "warn"
        return {
            "area_id": area_id,
            "label": label,
            "status": safe_status,
            "evidence": evidence,
            "route": route,
            "api": api,
            "actions": actions,
            "artifacts": artifacts or {},
        }

    external_ready = [
        row for row in adapter_rows
        if row.get("id") in {"reinvent4", "drugex", "diffsbdd", "colabfold"}
        and row.get("status") in {"available", "local_repo"}
    ]
    target_pack_ready = paths["target_pack"].exists()
    stage2_ready = paths["stage2_matrix"].exists()
    stage4_ready = paths["stage4_assets"].exists()
    stage7_ready = paths["stage7_manifest"].exists()
    stage8_report_ready = paths["stage8_report"].exists()
    evidence_pack_ready = paths["evidence_pack_zip"].exists()
    delivery_script_ready = (ROOT / "scripts" / "build_product_delivery.py").exists()
    scientific_report = Path(str(science.get("files", {}).get("report", "")))
    science_level = str(science.get("readiness_level", "planning_only"))
    has_real_scores = int(stage4_guide.get("score_count", 0) or 0) > 0

    focus_areas = [
        area(
            "real_demo_hardening",
            "真实 Demo 加硬",
            "ready" if len(catalog_targets) >= 4 and target_pack_ready and stage2_ready else ("warn" if len(catalog_targets) >= 4 else "missing"),
            f"catalog_targets={len(catalog_targets)}; target_pack={target_pack_ready}; stage2={stage2_ready}",
            "brief",
            f"/api/target-catalog",
            [
                "把 EGFR、BACE1、HIV-1 protease 等非流感靶点做成可演示 Target Pack。",
                f"当前项目可先运行 POST /api/projects/{name}/target-pack 或 Stage 8 Demo Runner。",
            ],
            {"target_pack": str(paths["target_pack"]), "stage2_matrix": str(paths["stage2_matrix"])},
        ),
        area(
            "generator_adapters",
            "多生成器适配",
            "ready" if len(external_ready) >= 2 else "warn",
            f"available_or_local_external_adapters={len(external_ready)}; total_adapters={len(adapter_rows)}",
            "stage8",
            "/api/generator-adapters",
            [
                "优先把 REINVENT4 / DrugEx 输出作为 Stage 3 CSV 导入。",
                "DiffSBDD 和 ColabFold 保持结构包适配边界，先不声称模型训练或真实药效。",
            ],
            {"stage3_assets": str(paths["stage3_assets"]), "candidate_csv": str(paths["candidates"])},
        ),
        area(
            "stage4_scientific_trust",
            "Stage 4 科学可信度",
            "ready" if science_level in {"strong_computational", "moderate_computational"} and has_real_scores else ("warn" if stage4_ready else "missing"),
            f"stage4_assets={stage4_ready}; docking_status={stage4_guide.get('docking_status')}; scores={stage4_guide.get('score_count')}; readiness={science_level}",
            "stage4",
            f"/api/projects/{name}/stage4/operator-guide?round={safe_round}",
            [
                "补 receptor preparation、口袋坐标、真实 docking 分数、Pose QC、控药/decoy 同池校准。",
                f"查看科学就绪度 GET /api/projects/{name}/scientific-readiness?round={safe_round}。",
            ],
            {"stage4_assets": str(paths["stage4_assets"]), "scientific_report": str(scientific_report)},
        ),
        area(
            "stage8_command_center",
            "Stage 8 总控",
            "ready" if preflight.get("status") == "ready" else ("warn" if preflight.get("can_auto_repair") else "missing"),
            f"preflight={preflight.get('status', '-')}; missing={preflight.get('missing_count', '-')}; overall={action_plan.get('overall_status', '-')}",
            "stage8",
            f"/api/projects/{name}/stage8/action-plan?round={safe_round}",
            [
                "把一键预检、一键修复、一键 Demo、一键导出作为产品指挥台主路径。",
                f"缺失项可运行 POST /api/projects/{name}/stage8/repair。",
            ],
            {"stage8_report": str(paths["stage8_report"]), "evidence_pack_zip": str(paths["evidence_pack_zip"])},
        ),
        area(
            "delivery_reproducibility",
            "交付复现",
            "ready" if stage7_ready and evidence_pack_ready else ("warn" if delivery_script_ready else "missing"),
            f"stage7_manifest={stage7_ready}; evidence_pack={evidence_pack_ready}; standard_delivery_script={delivery_script_ready}",
            "stage7",
            f"/api/projects/{name}/stage8/full-export",
            [
                "使用 Stage 7 交付包和标准交付脚本沉淀可复现材料。",
                "正式交付前运行 scripts/build_product_delivery.py 生成标准 ZIP 和 manifest。",
            ],
            {"stage7_manifest": str(paths["stage7_manifest"]), "evidence_pack_zip": str(paths["evidence_pack_zip"])},
        ),
        area(
            "scientific_report_layer",
            "科学报告解释层",
            "ready" if stage8_report_ready and scientific_report.exists() else ("warn" if scientific_report.exists() else "missing"),
            f"stage8_report={stage8_report_ready}; scientific_report={scientific_report.exists()}; claim_level={science.get('claim_level', '-')}",
            "stage8",
            f"/api/projects/{name}/stage8/report",
            [
                "报告必须分清计算证据、公开文献证据、验证规划和不能声称的边界。",
                "面向答辩时优先展示评审模式、科学就绪度和边界声明。",
            ],
            {"stage8_report": str(paths["stage8_report"]), "scientific_report": str(scientific_report)},
        ),
    ]
    counts = {
        "ready": len([row for row in focus_areas if row["status"] == "ready"]),
        "warn": len([row for row in focus_areas if row["status"] == "warn"]),
        "missing": len([row for row in focus_areas if row["status"] == "missing"]),
        "total": len(focus_areas),
    }
    overall = "ready" if counts["missing"] == 0 and counts["warn"] == 0 else ("ready_with_warnings" if counts["missing"] == 0 else "needs_work")
    first_gap = next((row for row in focus_areas if row["status"] != "ready"), focus_areas[0])
    return {
        "stage": "product",
        "command": "product-deep-integration",
        "project": name,
        "round": safe_round,
        "generated_at": datetime.now().isoformat(),
        "overall_status": overall,
        "summary": {
            **counts,
            "completion_percent": round((counts["ready"] / counts["total"]) * 100, 1) if counts["total"] else 0,
            "primary_focus": first_gap["area_id"],
            "primary_action": first_gap["actions"][0] if first_gap.get("actions") else "",
        },
        "focus_areas": focus_areas,
        "canonical_urls": {
            "app": "http://localhost:8765/",
            "stage8_page": "http://localhost:8765/#stage8",
            "target_catalog": "http://localhost:8765/api/target-catalog",
            "generator_adapters": "http://localhost:8765/api/generator-adapters",
            "stage4_operator_guide": f"http://localhost:8765/api/projects/{name}/stage4/operator-guide?round={safe_round}",
            "stage8_action_plan": f"http://localhost:8765/api/projects/{name}/stage8/action-plan?round={safe_round}",
            "deep_integration": f"http://localhost:8765/api/projects/{name}/product/deep-integration?round={safe_round}",
        },
        "boundary": [
            "Deep integration tracks computational virtual-screening product readiness only.",
            "It does not prove biological activity, potency, efficacy, toxicity, safety, dosing, clinical benefit, or therapeutic usefulness.",
        ],
    }

def target_pack_validation_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    pack_path = project_dir / "target_pack.json"
    target_pack = read_json_if_exists(pack_path)
    if not target_pack:
        target_pack = stage1_target_pack_data(project_dir, round_no)
    pocket = target_pack.get("pocket", {}) if isinstance(target_pack.get("pocket"), dict) else {}
    controls = target_pack.get("controls", {}) if isinstance(target_pack.get("controls"), dict) else {}
    evidence = target_pack.get("evidence_matrix", {}) if isinstance(target_pack.get("evidence_matrix"), dict) else {}
    stage2_matrix = evidence.get("stage2_matrix", []) if isinstance(evidence.get("stage2_matrix"), list) else []
    positive_controls = controls.get("positive", []) if isinstance(controls.get("positive"), list) else []
    reference_controls = controls.get("reference", []) if isinstance(controls.get("reference"), list) else []

    def ready_if(value: Any, label: str, remedy: str) -> Dict[str, Any]:
        return product_check(label, label, "ready" if bool(value) else "missing", str(value or ""), remedy if not value else "")

    checks = {
        "target_id": ready_if(target_pack.get("target_id"), "target_id", "Select a target from the catalog or create a custom target."),
        "target_name": ready_if(target_pack.get("target_name"), "target_name", "Add a target display name."),
        "pdb_id": ready_if(target_pack.get("pdb_id"), "pdb_id", "Attach a PDB ID or protein structure source."),
        "reference_ligand": product_check(
            "reference_ligand",
            "reference_ligand",
            "ready" if target_pack.get("reference_ligand") else "warn",
            str(target_pack.get("reference_ligand") or ""),
            "Add a co-crystal/reference ligand when available.",
        ),
        "controls": product_check(
            "controls",
            "controls",
            "ready" if positive_controls or reference_controls else "missing",
            f"positive={len(positive_controls)}; reference={len(reference_controls)}",
            "Add known positive/reference controls for calibration.",
        ),
        "pocket_box": product_check(
            "pocket_box",
            "pocket_box",
            "ready" if pocket.get("center") and pocket.get("size") else "warn",
            f"center={pocket.get('center', [])}; size={pocket.get('size', [])}; source={pocket.get('source', '')}",
            "Extract a co-crystal ligand pocket or provide reviewed center/size coordinates.",
        ),
        "evidence_matrix": product_check(
            "evidence_matrix",
            "evidence_matrix",
            "ready" if stage2_matrix else "warn",
            f"rows={len(stage2_matrix)}",
            "Run Stage 2 evidence matrix for literature/PDB/control context.",
        ),
        "boundary": product_check("boundary", "boundary", "ready", "; ".join(str(x) for x in target_pack.get("boundary", [])), ""),
    }
    statuses = [item["status"] for item in checks.values()]
    status = "invalid" if any(checks[key]["status"] == "missing" for key in ["target_id", "pdb_id", "controls"]) else ("ready_with_warnings" if any(s != "ready" for s in statuses) else "ready")
    next_actions = [item["remedy"] for item in checks.values() if item.get("remedy")]
    exports = project_dir / "exports"
    exports.mkdir(parents=True, exist_ok=True)
    validation_path = exports / f"round_{round_no}_target_pack_validation.json"
    payload = {
        "stage": 1,
        "command": "target-pack-validate",
        "project": name,
        "round": round_no,
        "status": status,
        "checks": checks,
        "next_actions": next_actions or ["Target pack is ready for computational screening review."],
        "target_pack": target_pack,
        "files": {
            "target_pack_json": str(pack_path),
            "validation_json": str(validation_path),
        },
        "boundary": [
            "Target Pack validation checks data completeness only.",
            "It does not prove biological activity, efficacy, toxicity, safety, dosing, or clinical usefulness.",
        ],
    }
    write_json(validation_path, payload)
    return payload

def _float_or_none(value: Any) -> Optional[float]:
    try:
        if value in {"", None}:
            return None
        return float(value)
    except Exception:
        return None

def _best_score(rows: List[Dict[str, Any]], panel_types: set) -> Optional[float]:
    values = []
    for row in rows:
        if str(row.get("panel_type", "")).lower() in panel_types:
            score = _float_or_none(row.get("docking_score") or row.get("score") or row.get("affinity"))
            if score is not None:
                values.append(score)
    return min(values) if values else None

def scientific_readiness_payload(name: str, project_dir: Path, round_no: int) -> Dict[str, Any]:
    stage4 = stage4_payload(name, project_dir, round_no)
    stage45 = stage45_payload(name, project_dir, round_no)
    stage46 = stage46_payload(name, project_dir, round_no)
    scores = stage45.get("scores", []) if isinstance(stage45.get("scores"), list) else []
    validation = stage45.get("validation", {}) if isinstance(stage45.get("validation"), dict) else {}
    benchmark = stage46.get("benchmark", {}) if isinstance(stage46.get("benchmark"), dict) else {}
    bench_metrics = benchmark.get("metrics", {}) if isinstance(benchmark.get("metrics"), dict) else {}

    control_types = {"positive_control", "reference_control", "control", "known_control"}
    decoy_types = {"decoy"}
    candidate_types = {"candidate"}
    best_control = _best_score(scores, control_types)
    best_decoy = _best_score(scores, decoy_types)
    best_candidate = _best_score(scores, candidate_types)
    separation = (best_decoy - best_control) if best_control is not None and best_decoy is not None else None
    candidate_delta = (best_candidate - best_control) if best_candidate is not None and best_control is not None else None
    has_pose_pass = any(str(row.get("pose_pass", "")).lower() in {"true", "pass", "passed", "yes", "1"} for row in scores)
    auc = _float_or_none(bench_metrics.get("roc_auc"))
    has_stage2 = (project_dir / "evidence" / "stage2_target_sources.csv").exists()
    has_rdkit = bool(stage4.get("has_descriptors"))
    has_controls = best_control is not None
    has_decoys = best_decoy is not None
    has_scores = bool(scores)

    checks = {
        "target_evidence": product_check("target_evidence", "target_evidence", "ready" if has_stage2 else "warn", str(project_dir / "evidence" / "stage2_target_sources.csv"), "Run Stage 2 evidence matrix."),
        "rdkit_descriptors": product_check("rdkit_descriptors", "rdkit_descriptors", "ready" if has_rdkit else "warn", str(stage4.get("files", {}).get("real_descriptors", "")), "Run Stage 4 real-library assets."),
        "control_decoy_separation": product_check(
            "control_decoy_separation",
            "control_decoy_separation",
            "ready" if separation is not None and separation > 0 else ("warn" if has_controls and has_decoys else "missing"),
            f"best_control={best_control}; best_decoy={best_decoy}; separation={separation}",
            "Run Stage 4.5 with positive controls and decoys in the same docking settings.",
        ),
        "candidate_vs_control": product_check(
            "candidate_vs_control",
            "candidate_vs_control",
            "ready" if candidate_delta is not None and candidate_delta <= 0 else ("warn" if candidate_delta is not None else "missing"),
            f"best_candidate={best_candidate}; best_control={best_control}; delta={candidate_delta}",
            "Use candidate-vs-control delta as a prioritization signal, not an efficacy claim.",
        ),
        "pose_qc": product_check("pose_qc", "pose_qc", "ready" if has_pose_pass else ("warn" if has_scores else "missing"), f"pose_pass={has_pose_pass}", "Run PoseBusters or import pose QC results."),
        "benchmark_auc": product_check("benchmark_auc", "benchmark_auc", "ready" if auc is not None and auc >= 0.7 else ("warn" if auc is not None else "missing"), str(auc), "Run Stage 4.6 retrospective benchmark."),
    }
    if has_controls and has_decoys and has_scores and separation is not None and separation > 0 and auc is not None and auc >= 0.7:
        readiness_level = "strong_computational"
    elif has_controls and has_decoys and (has_rdkit or has_scores):
        readiness_level = "moderate_computational"
    else:
        readiness_level = "planning_only"

    exports = project_dir / "exports"
    exports.mkdir(parents=True, exist_ok=True)
    readiness_path = exports / f"round_{round_no}_scientific_readiness.json"
    report_path = exports / f"round_{round_no}_scientific_readiness.md"
    report_lines = [
        f"# Scientific Readiness · {name} round {round_no}",
        "",
        f"- Readiness level: `{readiness_level}`",
        f"- Claim level: `computational_screening_only`",
        f"- Best control score: `{best_control}`",
        f"- Best candidate score: `{best_candidate}`",
        f"- Best decoy score: `{best_decoy}`",
        f"- Control-decoy separation: `{separation}`",
        f"- ROC-AUC: `{auc}`",
        "",
        "## Boundary",
        "",
        "- Computational screening only.",
        "- No potency, efficacy, toxicity, safety, dosing, clinical benefit, or therapeutic claim.",
    ]
    payload = {
        "stage": "4-8",
        "command": "scientific-readiness",
        "project": name,
        "round": round_no,
        "generated_at": datetime.now().isoformat(),
        "readiness_level": readiness_level,
        "claim_level": "computational_screening_only",
        "metrics": {
            "stage45_scored": len(scores),
            "best_control_docking_score": best_control,
            "best_candidate_docking_score": best_candidate,
            "best_decoy_docking_score": best_decoy,
            "control_decoy_separation": separation,
            "candidate_vs_best_control_delta": candidate_delta,
            "roc_auc": auc,
            "pose_pass_count": validation.get("docking", {}).get("pose_pass_count", "") if isinstance(validation.get("docking"), dict) else "",
        },
        "controls": {"has_positive_control": has_controls, "best_score": best_control},
        "decoys": {"has_decoy": has_decoys, "best_score": best_decoy},
        "candidates": {"has_candidate_score": best_candidate is not None, "best_score": best_candidate},
        "checks": checks,
        "next_actions": [item["remedy"] for item in checks.values() if item.get("status") != "ready" and item.get("remedy")][:8],
        "files": {"readiness_json": str(readiness_path), "report": str(report_path)},
        "boundary": [
            "Computational readiness summarizes virtual-screening evidence only.",
            "It does not prove biological activity, potency, efficacy, toxicity, safety, dosing, or clinical usefulness.",
        ],
    }
    write_text(report_path, "\n".join(report_lines) + "\n")
    write_json(readiness_path, payload)
    return payload

def stage8_full_export_bundle(name: str, project_dir: Path, round_no: int, title: str = "") -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    title = title or f"{name} Full Product Export"
    ensure_project_dirs(project_dir)
    warnings = stage8_try_generate_pack_prereqs(name, project_dir, safe_round, title)
    try:
        closed_loop_report = render_stage8_closed_loop_report(name, project_dir, safe_round, title)
    except Exception as exc:
        closed_loop_report = {"files": {}, "content": ""}
        warnings.append(f"Stage 8 report not regenerated: {exc}")
    system_health = product_system_health_payload(default_project=name)
    target_validation = target_pack_validation_payload(name, project_dir, safe_round)
    scientific = scientific_readiness_payload(name, project_dir, safe_round)
    evidence_pack = stage8_evidence_pack(name, project_dir, safe_round, title)
    command_center = stage8_payload(name, project_dir, safe_round)

    export_dir = project_dir / "exports"
    export_dir.mkdir(parents=True, exist_ok=True)
    health_path = export_dir / f"round_{safe_round}_system_health.json"
    command_center_path = export_dir / f"round_{safe_round}_stage8_command_center.json"
    manifest_path = export_dir / f"round_{safe_round}_full_product_manifest.json"
    zip_path = export_dir / f"round_{safe_round}_full_product_export.zip"
    write_json(health_path, system_health)
    write_json(command_center_path, command_center)

    include_map = {
        "system_health": health_path,
        "target_pack_validation": Path(target_validation["files"]["validation_json"]),
        "scientific_readiness": Path(scientific["files"]["readiness_json"]),
        "scientific_readiness_report": Path(scientific["files"]["report"]),
        "stage8_command_center": command_center_path,
        "stage8_closed_loop_report": Path(str(closed_loop_report.get("files", {}).get("report", ""))),
        "evidence_pack_zip": Path(str(evidence_pack.get("files", {}).get("zip", ""))),
        "target_pack_json": project_dir / "target_pack.json",
        "target_brief": project_dir / "briefs" / "target_brief.json",
        "stage2_matrix": project_dir / "evidence" / "stage2_target_sources.csv",
        "raw_candidates": project_dir / "candidates" / f"round_{safe_round}_candidates.csv",
        "ranked_candidates": project_dir / "ranked" / f"round_{safe_round}_ranked.csv",
        "stage4_assets": project_dir / "stage4" / f"round_{safe_round}_stage4_assets.json",
        "stage45_scores": project_dir / "stage4_5" / f"round_{safe_round}_control_docking_scores.csv",
        "stage46_benchmark": project_dir / "stage4_6" / f"round_{safe_round}_retrospective_benchmark.json",
        "stage5_dashboard_data": project_dir / "stage5" / "dashboard_data.json",
        "stage6_assets": project_dir / "stage6" / f"round_{safe_round}_validation_assets.json",
        "stage7_manifest": project_dir / "stage7" / f"round_{safe_round}_delivery_manifest.json",
    }
    included_files = [key for key, path in include_map.items() if path and Path(path).exists() and Path(path).is_file()]
    manifest = {
        "schema_version": "0.2",
        "stage": 8,
        "command": "stage8-full-export",
        "project": name,
        "round": safe_round,
        "title": title,
        "generated_at": datetime.now().isoformat(),
        "package_status": "ready" if target_validation.get("status") == "ready" and scientific.get("readiness_level") != "planning_only" else "ready_with_warnings",
        "included_files": included_files,
        "warnings": warnings,
        "system_health_status": system_health.get("overall_status"),
        "target_pack_status": target_validation.get("status"),
        "scientific_readiness_level": scientific.get("readiness_level"),
        "boundary": [
            "Computational product export for review, reproduction, and handoff.",
            "No potency, efficacy, toxicity, safety, dosing, clinical benefit, or therapeutic claim is created.",
        ],
    }
    write_json(manifest_path, manifest)
    include_map["full_export_manifest"] = manifest_path
    included_files = [key for key, path in include_map.items() if path and Path(path).exists() and Path(path).is_file()]

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for key, path in include_map.items():
            if not path:
                continue
            file_path = Path(path)
            if not file_path.exists() or not file_path.is_file():
                continue
            try:
                arcname = file_path.resolve().relative_to(project_dir.resolve()).as_posix()
            except ValueError:
                arcname = f"external/{file_path.name}"
            zf.write(file_path, arcname=arcname)
        zf.writestr(
            "README_FULL_EXPORT.txt",
            "\n".join(
                [
                    f"Project: {name}",
                    f"Round: {safe_round}",
                    f"Title: {title}",
                    "This package is computational screening evidence only.",
                    "It does not prove potency, efficacy, toxicity, safety, dosing, or clinical usefulness.",
                    "",
                    "Included files:",
                    *[f"- {key}" for key in included_files],
                ]
            )
            + "\n",
        )

    download_links = {
        "full_export_zip": project_artifact_url(name, project_dir, zip_path),
        "full_export_manifest": project_artifact_url(name, project_dir, manifest_path),
    }
    closed_loop_report_path = str(closed_loop_report.get("files", {}).get("report", ""))
    if closed_loop_report_path and Path(closed_loop_report_path).exists():
        download_links["closed_loop_report"] = project_artifact_url(name, project_dir, Path(closed_loop_report_path))
    for key, url in stage8_download_links(name, project_dir, safe_round).items():
        download_links.setdefault(key, url)

    return {
        "stage": 8,
        "command": "stage8-full-export",
        "project": name,
        "round": safe_round,
        "package_status": manifest["package_status"],
        "included_files": included_files,
        "warnings": warnings,
        "files": {"zip": str(zip_path), "manifest": str(manifest_path)},
        "system_health": system_health,
        "target_pack_validation": target_validation,
        "scientific_readiness": scientific,
        "evidence_pack": evidence_pack,
        "download_links": download_links,
        "boundary": manifest["boundary"],
    }

def stage8_try_generate_pack_prereqs(name: str, project_dir: Path, round_no: int, title: str) -> List[str]:
    warnings: List[str] = []
    if not (project_dir / "target_pack.json").exists():
        try:
            stage1_target_pack_data(project_dir, round_no)
        except SystemExit as exc:
            warnings.append(f"Target Pack not regenerated: {exc}")
        except Exception as exc:
            warnings.append(f"Target Pack not regenerated: {exc}")
    if not (project_dir / "stage5" / "dashboard_data.json").exists():
        try:
            stage5_dashboard(argparse.Namespace(project=str(project_dir), round=round_no, title=title or f"{name} Stage 5 Dashboard"))
        except SystemExit as exc:
            warnings.append(f"Stage 5 dashboard not regenerated: {exc}")
        except Exception as exc:
            warnings.append(f"Stage 5 dashboard not regenerated: {exc}")
    if not (project_dir / "stage6" / f"round_{round_no}_validation_assets.json").exists():
        try:
            stage6_validate(argparse.Namespace(project=str(project_dir), round=round_no, top=8))
        except SystemExit as exc:
            warnings.append(f"Stage 6 validation not regenerated: {exc}")
        except Exception as exc:
            warnings.append(f"Stage 6 validation not regenerated: {exc}")
    if not (project_dir / "stage7" / f"round_{round_no}_delivery_manifest.json").exists():
        try:
            stage7_package(argparse.Namespace(project=str(project_dir), round=round_no, title=title or f"{name} Stage 7 Delivery Package"))
        except SystemExit as exc:
            warnings.append(f"Stage 7 delivery package not regenerated: {exc}")
        except Exception as exc:
            warnings.append(f"Stage 7 delivery package not regenerated: {exc}")
    return warnings

def stage8_evidence_pack(name: str, project_dir: Path, round_no: int, title: str = "") -> Dict[str, Any]:
    safe_round = max(1, int(round_no))
    ensure_project_dirs(project_dir)
    title = title or f"{name} Stage 8 Evidence Pack"
    warnings = stage8_try_generate_pack_prereqs(name, project_dir, safe_round, title)
    try:
        closed_loop_report = render_stage8_closed_loop_report(name, project_dir, safe_round, title)
    except Exception as exc:
        closed_loop_report = {"files": {}, "content": ""}
        warnings.append(f"Stage 8 closed-loop report not regenerated: {exc}")

    command_center = stage8_payload(name, project_dir, safe_round)
    preflight = stage8_preflight_payload(name, project_dir, safe_round)
    guide = stage8_demo_guide_payload(name, project_dir, safe_round)
    export_dir = project_dir / "exports"
    export_dir.mkdir(parents=True, exist_ok=True)
    zip_path = export_dir / f"round_{safe_round}_evidence_pack.zip"
    manifest_path = export_dir / f"round_{safe_round}_evidence_manifest.json"

    include_map = {
        "target_pack": project_dir / "target_pack.json",
        "target_brief": project_dir / "briefs" / "target_brief.json",
        "generator_prompt": project_dir / "prompts" / "generator_prompt.md",
        "stage2_matrix": project_dir / "evidence" / "stage2_target_sources.csv",
        "stage2_assets": project_dir / "evidence" / "stage2_closed_loop_assets.json",
        "raw_candidates": project_dir / "candidates" / f"round_{safe_round}_candidates.csv",
        "candidate_scores": project_dir / "scores" / f"round_{safe_round}_scores.csv",
        "ranked_candidates": project_dir / "ranked" / f"round_{safe_round}_ranked.csv",
        "feedback": project_dir / "feedback" / f"round_{safe_round}_feedback.json",
        "stage3_assets": project_dir / "stage3" / f"round_{safe_round}_stage3_assets.json",
        "stage3_report": project_dir / "reports" / f"stage3_round_{safe_round}_report.md",
        "stage4_assets": project_dir / "stage4" / f"round_{safe_round}_stage4_assets.json",
        "stage4_receptor_package": project_dir / "stage4" / f"round_{safe_round}_receptor_package.json",
        "stage4_docking_plan": project_dir / "stage4" / f"round_{safe_round}_docking_plan.json",
        "stage4_real_descriptors": project_dir / "stage4" / f"round_{safe_round}_real_descriptors.csv",
        "stage4_benchmark_panel": project_dir / "stage4" / f"round_{safe_round}_benchmark_panel.csv",
        "stage4_report": project_dir / "reports" / f"stage4_round_{safe_round}_report.md",
        "stage5_dashboard_data": project_dir / "stage5" / "dashboard_data.json",
        "stage5_dashboard_html": project_dir / "stage5" / "index.html",
        "stage6_assets": project_dir / "stage6" / f"round_{safe_round}_validation_assets.json",
        "stage6_quality_gates": project_dir / "stage6" / f"round_{safe_round}_quality_gates.csv",
        "stage6_hit_triage": project_dir / "stage6" / f"round_{safe_round}_hit_triage.csv",
        "stage6_assay_queue": project_dir / "stage6" / f"round_{safe_round}_assay_queue.csv",
        "stage6_risk_register": project_dir / "stage6" / f"round_{safe_round}_risk_register.csv",
        "stage6_report": project_dir / "reports" / f"stage6_round_{safe_round}_validation_report.md",
        "stage7_manifest": project_dir / "stage7" / f"round_{safe_round}_delivery_manifest.json",
        "stage7_summary": project_dir / "stage7" / f"round_{safe_round}_executive_summary.md",
        "stage7_reproducibility": project_dir / "stage7" / f"round_{safe_round}_reproducibility.md",
        "stage7_checklist": project_dir / "stage7" / f"round_{safe_round}_investor_demo_checklist.csv",
        "stage7_frontend_spec": project_dir / "stage7" / "stage8_frontend_product_spec.md",
        "stage7_report": project_dir / "reports" / f"stage7_round_{safe_round}_delivery_report.md",
        "acceptance_json": project_dir / "acceptance" / "stage8_acceptance_report.json",
        "acceptance_md": project_dir / "acceptance" / "stage8_acceptance_report.md",
        "closed_loop_report": Path(str(closed_loop_report.get("files", {}).get("report", ""))) if closed_loop_report.get("files") else project_dir / "reports" / f"stage8_round_{safe_round}_closed_loop_report.md",
        "evidence_manifest": manifest_path,
    }
    included_files = [key for key, path in include_map.items() if key != "evidence_manifest" and path and Path(path).exists()]
    missing_files = [key for key, path in include_map.items() if key != "evidence_manifest" and path and not Path(path).exists()]
    included_files.append("evidence_manifest")
    package_status = "ready" if {"ranked_candidates", "stage7_manifest", "closed_loop_report"}.issubset(set(included_files)) else "ready_with_warnings"
    manifest = {
        "schema_version": "0.2",
        "stage": 8,
        "command": "stage8-evidence-pack",
        "project": name,
        "round": safe_round,
        "title": title,
        "generated_at": datetime.now().isoformat(),
        "package_status": package_status,
        "included_files": included_files + ["stage8_command_center", "stage8_preflight", "stage8_demo_guide"],
        "missing_files": missing_files,
        "warnings": warnings,
        "preflight_status": preflight.get("status"),
        "boundary": [
            "Evidence pack exports computational screening artifacts for review, demo, and handoff.",
            "It does not prove biological activity, potency, efficacy, toxicity, safety, dosing, or clinical usefulness.",
            "Docking evidence is included only when generated by local tools or explicitly imported as external scores.",
        ],
    }
    write_json(manifest_path, manifest)

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for key, path in include_map.items():
            file_path = Path(path)
            if not file_path.exists() or not file_path.is_file():
                continue
            zf.write(file_path, arcname=file_path.relative_to(project_dir).as_posix())
        zf.writestr("stage8_command_center.json", json.dumps(command_center, ensure_ascii=False, indent=2))
        zf.writestr("stage8_preflight.json", json.dumps(preflight, ensure_ascii=False, indent=2))
        zf.writestr("stage8_demo_guide.json", json.dumps(guide, ensure_ascii=False, indent=2))
        zf.writestr(
            "README.md",
            "\n".join(
                [
                    f"# {title}",
                    "",
                    f"- Project: `{name}`",
                    f"- Round: `{safe_round}`",
                    f"- Package status: `{package_status}`",
                    "",
                    "This package contains computational screening workflow artifacts for review and product demonstration.",
                    "It does not establish efficacy, safety, toxicity, dosing, clinical benefit, or therapeutic usefulness.",
                ]
            )
            + "\n",
        )

    download_links = {
        "evidence_pack_zip": project_artifact_url(name, project_dir, zip_path),
        "evidence_manifest": project_artifact_url(name, project_dir, manifest_path),
    }
    if include_map["closed_loop_report"].exists():
        download_links["closed_loop_report"] = project_artifact_url(name, project_dir, include_map["closed_loop_report"])
    for key, url in stage8_download_links(name, project_dir, safe_round).items():
        download_links.setdefault(key, url)

    return {
        "stage": 8,
        "command": "stage8-evidence-pack",
        "project": name,
        "round": safe_round,
        "title": title,
        "generated_at": manifest["generated_at"],
        "package_status": package_status,
        "included_files": manifest["included_files"],
        "missing_files": missing_files,
        "warnings": warnings,
        "files": {
            "zip": str(zip_path),
            "manifest": str(manifest_path),
            "closed_loop_report": str(include_map["closed_loop_report"]) if include_map["closed_loop_report"].exists() else "",
        },
        "download_links": download_links,
        "preflight": preflight,
        "guide": guide,
        "stage8": command_center,
        "boundary": manifest["boundary"],
    }

# ── Stage 1: 靶点筛选与靶点 Brief ────────────────────────────────────────
@app.get("/api/projects/{name}/stage1")
def api_stage1_status(name: str):
    p = get_project_dir(name)
    target_selection_path = p / "targets" / "target_selection.csv"
    target_report_path = p / "reports" / "target_selection.md"
    brief_path = p / "briefs" / "target_brief.json"
    prompt_path = p / "prompts" / "generator_prompt.md"
    target_pack_path = p / "target_pack.json"
    target_pack_report_path = p / "reports" / "target_pack_report.md"
    return {
        "stage": 1,
        "project": name,
        "has_target_selection": target_selection_path.exists(),
        "has_target_report": target_report_path.exists(),
        "has_brief": brief_path.exists(),
        "has_prompt": prompt_path.exists(),
        "has_target_pack": target_pack_path.exists(),
        "target_selection": csv_to_dicts(target_selection_path),
        "target_report": read_text_if_exists(target_report_path),
        "brief": read_json(brief_path) if brief_path.exists() else {},
        "prompt": read_text_if_exists(prompt_path),
        "target_pack": read_json_if_exists(target_pack_path),
        "files": {
            "target_selection": str(target_selection_path),
            "target_report": str(target_report_path),
            "target_brief": str(brief_path),
            "generator_prompt": str(prompt_path),
            "target_pack_json": str(target_pack_path),
            "target_pack_report": str(target_pack_report_path),
        },
    }

@app.post("/api/projects/{name}/stage1/target-select")
def api_stage1_target_select(name: str, req: Stage1TargetSelectRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    try:
        target_select(
            argparse.Namespace(
                project=str(p),
                disease=req.disease or "influenza_a",
                target=req.target or None,
                top=max(1, min(int(req.top), 50)),
                catalog=req.catalog or None,
            )
        )
    except SystemExit as exc:
        http_from_cli_exit(exc)

    target_selection_path = p / "targets" / "target_selection.csv"
    target_report_path = p / "reports" / "target_selection.md"
    return {
        "stage": 1,
        "command": "target-select",
        "project": name,
        "disease": req.disease,
        "targets": csv_to_dicts(target_selection_path),
        "report": read_text_if_exists(target_report_path),
        "files": {
            "target_selection": str(target_selection_path),
            "target_report": str(target_report_path),
        },
    }

@app.post("/api/projects/{name}/stage1/brief-from-target")
def api_stage1_brief_from_target(name: str, req: Stage1BriefFromTargetRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    try:
        brief_from_target(
            argparse.Namespace(
                project=str(p),
                disease=req.disease or "influenza_a",
                target=req.target or None,
                catalog=req.catalog or None,
                free_text=req.free_text,
                max_heavy_atoms=req.max_heavy_atoms,
                max_molecular_weight=req.max_molecular_weight,
                force=req.force,
            )
        )
    except SystemExit as exc:
        http_from_cli_exit(exc)

    brief_path = p / "briefs" / "target_brief.json"
    prompt_path = p / "prompts" / "generator_prompt.md"
    config_path = p / "config.json"
    brief = read_json(brief_path)
    return {
        "stage": 1,
        "command": "brief-from-target",
        "project": name,
        "disease": req.disease,
        "selected_target": brief.get("target_catalog_id", ""),
        "brief": brief,
        "prompt": read_text_if_exists(prompt_path),
        "config": read_json(config_path) if config_path.exists() else {},
        "files": {
            "target_brief": str(brief_path),
            "generator_prompt": str(prompt_path),
            "config": str(config_path),
        },
    }

# ── Stage 2: 靶点证据矩阵与闭环资产 ───────────────────────────────────────
@app.get("/api/projects/{name}/stage2")
def api_stage2_status(name: str):
    p = get_project_dir(name)
    matrix_path = p / "evidence" / "stage2_target_sources.csv"
    assets_path = p / "evidence" / "stage2_closed_loop_assets.json"
    report_path = p / "reports" / "stage2_evidence_report.md"
    return {
        "stage": 2,
        "project": name,
        "has_matrix": matrix_path.exists(),
        "has_assets": assets_path.exists(),
        "has_report": report_path.exists(),
        "matrix": csv_to_dicts(matrix_path),
        "assets": read_json(assets_path) if assets_path.exists() else {},
        "report": read_text_if_exists(report_path),
        "files": {
            "matrix": str(matrix_path),
            "assets": str(assets_path),
            "report": str(report_path),
        },
    }

@app.post("/api/projects/{name}/stage2/evidence")
def api_stage2_evidence(name: str, req: Stage2EvidenceRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    try:
        evidence_stage2(
            argparse.Namespace(
                project=str(p),
                disease=req.disease or "influenza",
                target=req.target or None,
                top=max(1, min(int(req.top), 50)),
                catalog=req.catalog or None,
                sources=req.sources or None,
                evidence_dir=req.evidence_dir or None,
                refresh=bool(req.refresh),
                retmax=max(1, min(int(req.retmax), 100)),
                timeout=max(1, min(int(req.timeout), 120)),
                offline=bool(req.offline),
            )
        )
    except SystemExit as exc:
        http_from_cli_exit(exc)

    matrix_path = p / "evidence" / "stage2_target_sources.csv"
    assets_path = p / "evidence" / "stage2_closed_loop_assets.json"
    report_path = p / "reports" / "stage2_evidence_report.md"
    return {
        "stage": 2,
        "command": "evidence-stage2",
        "project": name,
        "disease": req.disease,
        "matrix": csv_to_dicts(matrix_path),
        "assets": read_json(assets_path) if assets_path.exists() else {},
        "report": read_text_if_exists(report_path),
        "files": {
            "matrix": str(matrix_path),
            "assets": str(assets_path),
            "report": str(report_path),
        },
    }

# ── Stage 4: 真实化学校验与对接准备 ───────────────────────────────────────
@app.get("/api/stage4/capabilities")
def api_stage4_capabilities():
    capabilities = stage4_capabilities()
    capabilities["stage"] = 4
    return capabilities

@app.get("/api/stage4/presets")
def api_stage4_presets():
    return {
        "stage": 4,
        "presets": stage4_target_presets(),
        "boundary": [
            "疾病靶点预设只提供靶点、结构和控药建议；真实 docking box 必须有可复核坐标来源。",
            "Vina 官方 1IEP 预设只用于验证本地工具链，不代表任何疾病项目的有效靶点。",
        ],
    }

@app.post("/api/stage4/smoke-test")
def api_stage4_smoke_test(req: Stage4SmokeTestRequest):
    capabilities = stage4_capabilities()
    executables = capabilities.get("executables", {}) if isinstance(capabilities.get("executables"), dict) else {}
    vina_status = executables.get("vina", {}) if isinstance(executables.get("vina", {}), dict) else {}
    vina_path = str(vina_status.get("path", "") or "")
    example_dir = DEFAULT_EVAL_TOOLS / "AutoDock-Vina" / "example" / "basic_docking" / "solution"
    receptor = example_dir / "1iep_receptor.pdbqt"
    ligand = example_dir / "1iep_ligand.pdbqt"
    box = {"center": [15.19, 53.903, 16.917], "size": [20.0, 20.0, 20.0]}
    base = {
        "stage": 4,
        "example": {
            "name": "AutoDock Vina official 1IEP example",
            "receptor": str(receptor),
            "ligand": str(ligand),
            "box": box,
            "boundary": "该烟测只验证本机 Vina 工具链，不代表疾病靶点或药效结论。",
        },
    }
    if not receptor.exists() or not ligand.exists():
        return {
            **base,
            "status": "missing_example_assets",
            "message": "AutoDock Vina 官方示例 receptor/ligand 文件不完整，无法烟测。",
        }
    if not vina_path:
        return {
            **base,
            "status": "not_available",
            "message": "未检测到 vina 可执行程序；可以先用预设生成 docking-ready 资产或导入外部真实分数 CSV。",
            "capabilities": capabilities,
        }
    command = [
        vina_path,
        "--receptor",
        str(receptor),
        "--ligand",
        str(ligand),
        "--center_x",
        str(box["center"][0]),
        "--center_y",
        str(box["center"][1]),
        "--center_z",
        str(box["center"][2]),
        "--size_x",
        str(box["size"][0]),
        "--size_y",
        str(box["size"][1]),
        "--size_z",
        str(box["size"][2]),
        "--score_only",
    ]
    try:
        proc = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=max(1, min(int(req.timeout), 120)),
        )
    except Exception as exc:
        return {
            **base,
            "status": "failed",
            "message": f"Vina 烟测启动失败: {exc}",
            "command": command,
        }
    output = (proc.stdout or "") + "\n" + (proc.stderr or "")
    match = re.search(r"affinity[:\s]+(-?\d+(?:\.\d+)?)", output, flags=re.IGNORECASE)
    if not match:
        match = re.search(r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+", output, flags=re.MULTILINE)
    return {
        **base,
        "status": "passed" if proc.returncode == 0 else "failed",
        "message": "Vina 官方 1IEP 示例烟测已执行。" if proc.returncode == 0 else "Vina 已找到，但官方示例烟测返回非零状态。",
        "command": command,
        "returncode": proc.returncode,
        "score": match.group(1) if match else "",
        "stdout_tail": (proc.stdout or "")[-1200:],
        "stderr_tail": (proc.stderr or "")[-1200:],
    }

@app.get("/api/projects/{name}/stage4")
def api_stage4_status(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage4_payload(name, p, safe_round)

@app.get("/api/projects/{name}/stage4/operator-guide")
def api_stage4_operator_guide(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage4_operator_guide_payload(name, p, safe_round)

@app.get("/api/projects/{name}/doctor")
def api_project_doctor(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage4_project_doctor(name, p, safe_round)

@app.post("/api/projects/{name}/stage4/repair")
def api_stage4_repair(name: str, req: Stage4RepairRequest):
    p = get_project_dir(name)
    safe_round = max(1, int(req.round))
    return stage4_repair_project(name, p, safe_round, persist=bool(req.persist))

@app.post("/api/projects/{name}/stage4/pocket-pack")
def api_stage4_pocket_pack(name: str, req: Stage4PocketPackRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    pocket_pack = stage4_build_pocket_pack(
        p,
        max(1, int(req.round)),
        target_id=req.target or "",
        pdb_id=req.pdb_id or "",
        reference_ligand=req.reference_ligand or "",
    )
    return {
        "stage": 4,
        "command": "stage4-pocket-pack",
        "project": name,
        "round": max(1, int(req.round)),
        "pocket": pocket_pack.get("pocket", {}),
        "receptor": pocket_pack.get("receptor", {}),
        "files": pocket_pack.get("files", {}),
    }

@app.get("/api/projects/{name}/stage45")
def api_stage45_status(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage45_payload(name, p, safe_round)

@app.get("/api/projects/{name}/stage46")
def api_stage46_status(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage46_payload(name, p, safe_round)

@app.post("/api/projects/{name}/stage4/real")
def api_stage4_real(name: str, req: Stage4RealRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    backend = req.docking_backend if req.docking_backend in {"auto", "vina", "gnina"} else "auto"
    config_path = p / "config.json"
    config = read_json(config_path) if config_path.exists() else default_config()
    target = config.get("target", {})
    if not isinstance(target, dict):
        target = {}
    center = parse_triplet(req.pocket_center) if req.pocket_center.strip() else []
    size = parse_triplet(req.pocket_size) if req.pocket_size.strip() else []
    if center or size or req.pocket_source:
        pocket = target.get("pocket", {})
        if not isinstance(pocket, dict):
            pocket = {}
        if center:
            pocket["center"] = center
        if size:
            pocket["size"] = size
        if req.pocket_source:
            pocket["source"] = req.pocket_source
        target["pocket"] = pocket
        config["target"] = target
        write_json(config_path, config)
    try:
        stage4_real(
            argparse.Namespace(
                project=str(p),
                round=safe_round,
                target=req.target or "",
                input_csv=req.input_csv or None,
                controls_csv=req.controls_csv or None,
                receptor_pdb=req.receptor_pdb or None,
                pdb_id=req.pdb_id or "",
                pocket_center=req.pocket_center,
                pocket_size=req.pocket_size,
                pocket_source=req.pocket_source,
                fetch_receptor=bool(req.fetch_receptor),
                top=max(1, min(int(req.top), 200)),
                decoys=max(0, min(int(req.decoys), 200)),
                max_conformers=max(0, min(int(req.max_conformers), 20)),
                seed=int(req.seed),
                no_sdf=bool(req.no_sdf),
                render_2d=bool(req.render_2d),
                docking_backend=backend,
                run_docking=bool(req.run_docking),
                docking_timeout=max(1, min(int(req.docking_timeout), 7200)),
                rescore=bool(req.rescore),
                rank_top=max(0, min(int(req.rank_top), 200)),
                feedback_top=max(0, min(int(req.feedback_top), 200)),
                external_scores=req.external_scores or None,
            )
        )
    except SystemExit as exc:
        http_from_cli_exit(exc)

    return stage4_payload(name, p, safe_round, command="stage4-real")

@app.post("/api/projects/{name}/stage45/validate")
def api_stage45_validate(name: str, req: Stage45ValidateRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    backend = req.docking_backend if req.docking_backend in {"auto", "vina", "gnina"} else "auto"
    try:
        stage45_validate_controls(
            argparse.Namespace(
                project=str(p),
                round=safe_round,
                target=req.target or "",
                top_candidates=max(1, min(int(req.top_candidates), 200)),
                decoys=max(0, min(int(req.decoys), 200)),
                controls_csv=req.controls_csv or None,
                docking_backend=backend,
                docking_timeout=max(1, min(int(req.docking_timeout), 7200)),
                seed=int(req.seed),
                no_docking=bool(req.no_docking),
            )
        )
    except SystemExit as exc:
        http_from_cli_exit(exc)

    return stage45_payload(name, p, safe_round, command="stage4-validate-controls")

@app.post("/api/projects/{name}/stage46/benchmark")
def api_stage46_benchmark(name: str, req: Stage46BenchmarkRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    try:
        stage46_retrospective_benchmark(
            argparse.Namespace(
                project=str(p),
                round=safe_round,
                positive_types=req.positive_types or "positive_control,reference_control,control",
                negative_types=req.negative_types or "decoy",
                top_k=req.top_k or "1,3,5,10",
            )
        )
    except SystemExit as exc:
        http_from_cli_exit(exc)

    return stage46_payload(name, p, safe_round, command="stage4-retrospective-benchmark")

@app.get("/api/projects/{name}/stage4/images/{round_no}/{filename}")
def api_stage4_image(name: str, round_no: int, filename: str):
    p = get_project_dir(name)
    safe_file = Path(filename).name
    if safe_file != filename or not safe_file.lower().endswith(".png"):
        raise HTTPException(400, "非法图片文件名")
    image_path = p / "stage4" / f"round_{max(1, int(round_no))}_2d" / safe_file
    if not image_path.exists():
        raise HTTPException(404, "Stage 4 图片不存在")
    return FileResponse(image_path, media_type="image/png")

# ── Stage 5: 产品化 Dashboard ───────────────────────────────────────────
@app.get("/api/projects/{name}/stage5")
def api_stage5_status(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    data_path = p / "stage5" / "dashboard_data.json"
    html_path = p / "stage5" / "index.html"
    report_path = p / "reports" / "stage5_dashboard_report.md"
    existing_dashboard = read_json_if_exists(data_path)
    title = str(existing_dashboard.get("title", "") or f"{name} Stage 5 Dashboard")
    dashboard = build_stage5_dashboard_data(p, safe_round, title)
    return {
        "stage": 5,
        "project": name,
        "round": safe_round,
        "has_dashboard": data_path.exists(),
        "has_html": html_path.exists(),
        "has_report": report_path.exists(),
        "dashboard": dashboard,
        "report": read_text_if_exists(report_path),
        "files": {
            "dashboard_data": str(data_path),
            "dashboard_html": str(html_path),
            "report": str(report_path),
        },
    }

@app.post("/api/projects/{name}/stage5/dashboard")
def api_stage5_dashboard(name: str, req: Stage5DashboardRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    title = req.title or f"{name} Stage 5 Dashboard"
    try:
        stage5_dashboard(argparse.Namespace(project=str(p), round=safe_round, title=title))
    except SystemExit as exc:
        http_from_cli_exit(exc)

    data_path = p / "stage5" / "dashboard_data.json"
    html_path = p / "stage5" / "index.html"
    report_path = p / "reports" / "stage5_dashboard_report.md"
    return {
        "stage": 5,
        "command": "stage5-dashboard",
        "project": name,
        "round": safe_round,
        "has_dashboard": data_path.exists(),
        "has_html": html_path.exists(),
        "has_report": report_path.exists(),
        "dashboard": read_json(data_path) if data_path.exists() else {},
        "report": read_text_if_exists(report_path),
        "files": {
            "dashboard_data": str(data_path),
            "dashboard_html": str(html_path),
            "report": str(report_path),
        },
    }

# ── Stage 6: 验证运营 ───────────────────────────────────────────────────
@app.get("/api/projects/{name}/stage6")
def api_stage6_status(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage6_payload(name, p, safe_round)

@app.post("/api/projects/{name}/stage6/validate")
def api_stage6_validate(name: str, req: Stage6ValidateRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    top = max(1, min(int(req.top), 50))
    try:
        stage6_validate(argparse.Namespace(project=str(p), round=safe_round, top=top))
    except SystemExit as exc:
        http_from_cli_exit(exc)

    return stage6_payload(name, p, safe_round, command="stage6-validate")

# ── Stage 7: 交付包 ─────────────────────────────────────────────────────
@app.get("/api/projects/{name}/stage7")
def api_stage7_status(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage7_payload(name, p, safe_round)

@app.post("/api/projects/{name}/stage7/package")
def api_stage7_package(name: str, req: Stage7PackageRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    title = req.title or f"{name} Stage 7 Delivery Package"
    try:
        stage7_package(argparse.Namespace(project=str(p), round=safe_round, title=title))
    except SystemExit as exc:
        http_from_cli_exit(exc)

    return stage7_payload(name, p, safe_round, command="stage7-package")

# ── Stage 8: 产品指挥台 ─────────────────────────────────────────────────
@app.get("/api/projects/{name}/demo-doctor")
def api_demo_doctor(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return demo_doctor_payload(name, p, safe_round)

@app.get("/api/projects/{name}/stage8/command-center")
def api_stage8_command_center(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage8_payload(name, p, safe_round)

@app.get("/api/projects/{name}/stage8/preflight")
def api_stage8_preflight(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage8_preflight_payload(name, p, safe_round)

@app.get("/api/projects/{name}/stage8/demo-guide")
def api_stage8_demo_guide(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage8_demo_guide_payload(name, p, safe_round)

@app.get("/api/projects/{name}/stage8/review-mode")
def api_stage8_review_mode(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage8_review_mode_payload(name, p, safe_round)

@app.get("/api/projects/{name}/stage8/action-plan")
def api_stage8_action_plan(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return stage8_action_plan_payload(name, p, safe_round)

@app.get("/api/projects/{name}/product/deep-integration")
def api_product_deep_integration(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    safe_round = max(1, int(round_no))
    return product_deep_integration_payload(name, p, safe_round)

@app.post("/api/stage8/acceptance-demo")
def api_stage8_acceptance_demo(req: Stage8AcceptanceDemoRequest):
    return create_stage8_acceptance_demo(req)

@app.post("/api/projects/{name}/stage8/demo-package")
def api_stage8_demo_package(name: str, req: Stage8DemoPackageRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    top = max(1, min(int(req.top), 50))
    title = req.title or f"{name} Stage 8 Demo Package"
    try:
        stage5_dashboard(argparse.Namespace(project=str(p), round=safe_round, title=title))
        stage6_validate(argparse.Namespace(project=str(p), round=safe_round, top=top))
        stage7_package(argparse.Namespace(project=str(p), round=safe_round, title=title))
    except SystemExit as exc:
        http_from_cli_exit(exc)

    payload = stage8_payload(name, p, safe_round)
    payload["command"] = "stage8-demo-package"
    payload["generated"] = {
        "stage5": (p / "stage5" / "index.html").exists(),
        "stage6": (p / "stage6" / f"round_{safe_round}_validation_assets.json").exists(),
        "stage7": (p / "stage7" / f"round_{safe_round}_delivery_manifest.json").exists(),
    }
    payload["download_links"] = stage8_download_links(name, p, safe_round)
    return payload

@app.post("/api/projects/{name}/stage8/demo-runner")
def api_stage8_demo_runner(name: str, req: Stage8DemoRunnerRequest):
    p = get_project_dir(name)
    return stage8_demo_runner(name, p, req)

@app.post("/api/projects/{name}/stage8/repair")
def api_stage8_repair(name: str, req: Stage8RepairRequest):
    p = get_project_dir(name)
    runner_req = Stage8DemoRunnerRequest(
        round=max(1, int(req.round)),
        disease=req.disease or "甲流",
        target=req.target or "influenza_a_h1n1_na",
        source_kind=req.source_kind or "url",
        source=req.source or "https://www.rcsb.org/structure/3TI6",
        target_hint=req.target_hint or "neuraminidase oseltamivir pocket",
        candidates=max(1, min(int(req.candidates), 200)),
        top=max(1, min(int(req.top), 50)),
        decoys=max(0, min(int(req.decoys), 200)),
        run_docking=bool(req.run_docking),
        docking_backend=req.docking_backend if req.docking_backend in {"auto", "vina", "gnina"} else "auto",
        docking_timeout=max(1, min(int(req.docking_timeout), 7200)),
        title=f"{name} Stage 8 Repair",
    )
    runner = stage8_demo_runner(name, p, runner_req)
    return {
        "stage": 8,
        "command": "stage8-repair",
        "project": name,
        "round": runner_req.round,
        "status": runner.get("status", "completed"),
        "runner": runner,
        "preflight": stage8_preflight_payload(name, p, runner_req.round),
        "boundary": [
            "Auto repair generates or refreshes computational workflow artifacts only.",
            "It does not prove potency, efficacy, toxicity, safety, dosing, or clinical benefit.",
        ],
    }

@app.post("/api/projects/{name}/stage8/report")
def api_stage8_report(name: str, req: Stage8ReportRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    return render_stage8_closed_loop_report(name, p, safe_round, req.title or "")

@app.post("/api/projects/{name}/stage8/evidence-pack")
def api_stage8_evidence_pack(name: str, req: Stage8EvidencePackRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    safe_round = max(1, int(req.round))
    return stage8_evidence_pack(name, p, safe_round, req.title or "")

@app.post("/api/projects/{name}/stage8/full-export")
def api_stage8_full_export(name: str, req: Stage8FullExportRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    return stage8_full_export_bundle(name, p, max(1, int(req.round)), req.title or "")

@app.post("/api/projects/{name}/export-package")
def api_export_package(name: str, req: ExportPackageRequest):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    return stage8_export_bundle(p, max(1, int(req.round)))

@app.get("/api/projects/{name}/target-pack/validate")
def api_target_pack_validate(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    return target_pack_validation_payload(name, p, max(1, int(round_no)))

@app.get("/api/projects/{name}/scientific-readiness")
def api_scientific_readiness(name: str, round_no: int = Query(1, alias="round")):
    p = get_project_dir(name)
    ensure_project_dirs(p)
    return scientific_readiness_payload(name, p, max(1, int(round_no)))

def product_job_snapshot(job_id: str) -> Dict[str, Any]:
    with job_lock:
        job = job_status.get(job_id)
        if not job:
            raise HTTPException(404, "Job not found")
        return dict(job)

def run_product_job(job_id: str, req: ProductJobRequest):
    with job_lock:
        if job_id in job_status:
            job_status[job_id]["status"] = "running"
            job_status[job_id]["started_at"] = datetime.now().isoformat()
            job_status[job_id].setdefault("events", []).append({"at": datetime.now().isoformat(), "message": "job started"})
    try:
        project_name = safe_name(req.project)
        if not project_name:
            raise HTTPException(400, "project is required")
        project_dir = get_project_dir(project_name)
        safe_round = max(1, int(req.round))
        payload = req.payload if isinstance(req.payload, dict) else {}
        if req.task == "stage8_full_export":
            result = stage8_full_export_bundle(project_name, project_dir, safe_round, str(payload.get("title") or ""))
        elif req.task == "scientific_readiness":
            result = scientific_readiness_payload(project_name, project_dir, safe_round)
        elif req.task == "target_pack_validate":
            result = target_pack_validation_payload(project_name, project_dir, safe_round)
        else:
            raise HTTPException(400, f"Unsupported job task: {req.task}")
        with job_lock:
            job_status[job_id]["status"] = "completed"
            job_status[job_id]["finished_at"] = datetime.now().isoformat()
            job_status[job_id]["result"] = result
            job_status[job_id].setdefault("events", []).append({"at": datetime.now().isoformat(), "message": "job completed"})
    except Exception as exc:
        with job_lock:
            job_status[job_id]["status"] = "failed"
            job_status[job_id]["finished_at"] = datetime.now().isoformat()
            job_status[job_id]["error"] = str(exc)
            job_status[job_id]["traceback"] = traceback.format_exc()
            job_status[job_id].setdefault("events", []).append({"at": datetime.now().isoformat(), "message": f"job failed: {exc}"})

@app.post("/api/jobs")
def api_create_product_job(req: ProductJobRequest):
    job_id = f"job_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    project_name = safe_name(req.project)
    with job_lock:
        job_status[job_id] = {
            "job_id": job_id,
            "task": req.task,
            "project": project_name,
            "round": max(1, int(req.round)),
            "status": "queued",
            "created_at": datetime.now().isoformat(),
            "events": [{"at": datetime.now().isoformat(), "message": "job queued"}],
        }
    threading.Thread(target=run_product_job, args=(job_id, req), daemon=True).start()
    return product_job_snapshot(job_id)

@app.get("/api/jobs")
def api_list_product_jobs():
    with job_lock:
        jobs = [dict(job) for job in job_status.values()]
    jobs.sort(key=lambda row: row.get("created_at", ""), reverse=True)
    return {"jobs": jobs[:100], "count": len(jobs)}

@app.get("/api/jobs/{job_id}")
def api_get_product_job(job_id: str):
    return product_job_snapshot(job_id)

# ── 工作流运行 ────────────────────────────────────────────────────────────
@app.post("/api/projects/{name}/run")
def api_run_workflow(name: str, req: RunRequest):
    p = get_project_dir(name)
    if not (p / "config.json").exists(): raise HTTPException(400, "请先初始化项目")
    run_id = f"{name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    start_run(run_id)

    def _run():
        try:
            cfg = json.loads((p / "config.json").read_text())
            # 清除旧轮次文件
            import glob
            for sd in ["candidates","scores","ranked","feedback","reports"]:
                d = p / sd
                if d.exists():
                    for f in d.glob("round_*"): f.unlink()
            sd2 = p / "seeds"
            if sd2.exists():
                for f in sd2.glob("round_*_seeds.csv"):
                    if f.stem != "round_0_seeds": f.unlink()
            for rn in range(1, req.rounds + 1):
                emit(run_id, "round_start", {"round": rn, "total": req.rounds})
                # generate
                emit(run_id, "step", {"round": rn, "step": "generate", "message": "生成候选分子..."})
                seeds = load_seed_smiles(p, rn)
                import random
                rng = random.Random(rn * 7919 + req.n)
                pool = []
                for parent in seeds:
                    for sub in SUBSTITUENTS:
                        pool.append((aromatic_variant(parent, sub), parent, "substitution"))
                        pool.append((chain_variant(parent, sub), parent, "extension"))
                for mol in TEMPLATE_MOLECULES: pool.append((mol, rng.choice(seeds), "template"))
                rng.shuffle(pool)
                seen = set(); candidates = []
                for smi, parent, source in pool:
                    if smi in seen: continue
                    seen.add(smi); idx = len(candidates) + 1
                    candidates.append({"round": rn, "id": f"r{rn:02d}_{idx:05d}", "smiles": smi, "parent": parent, "source": source})
                    if len(candidates) >= req.n: break
                write_csv(p / "candidates" / f"round_{rn}_candidates.csv", candidates, ["round","id","smiles","parent","source"])
                emit(run_id, "step", {"round": rn, "step": "generate", "message": f"生成 {len(candidates)} 个候选分子", "count": len(candidates)})
                # score
                emit(run_id, "step", {"round": rn, "step": "score", "message": "评分计算中..."})
                ext_by_id, ext_by_smiles = ({}, {}) if not req.external_scores else load_external_scores(req.external_scores)
                weights = cfg.get("weights", DEFAULT_WEIGHTS)
                if not isinstance(weights, dict): weights = dict(DEFAULT_WEIGHTS)
                scored = []
                for row in candidates:
                    smi = row["smiles"]; sc = estimate_scores(smi, seeds); note = "proxy"
                    external = ext_by_id.get(row["id"]) or ext_by_smiles.get(smi)
                    if external:
                        from ai_mol_loop import normalized_external_docking, parse_pose_value
                        dp = normalized_external_docking(external.get("docking_score") or external.get("affinity") or external.get("gnina_score") or external.get("vina_score"))
                        if dp: sc["docking"], sc["raw_docking_kcal_mol"] = dp; note = "proxy+external_docking"
                        pv = parse_pose_value(external.get("pose_pass") or external.get("pose_score") or external.get("posebusters_pass"))
                        if pv is not None: sc["pose"] = pv; note += "+external_pose"
                    total = weighted_total(sc, weights)
                    scored.append({"round": rn, "id": row["id"], "smiles": smi, "parent": row["parent"], "source": row["source"],
                        "validity_proxy": round(sc["validity"],4), "qed_proxy": round(sc["qed"],4),
                        "sa_proxy": round(sc["sa"],4), "lipinski_proxy": round(sc["lipinski"],4),
                        "docking_proxy": round(sc["docking"],4), "pose_proxy": round(sc["pose"],4),
                        "novelty_proxy": round(sc["novelty"],4),
                        "raw_docking_kcal_mol": round(sc.get("raw_docking_kcal_mol", 0), 3),
                        "heavy_atoms_proxy": int(sc["heavy_atoms_proxy"]),
                        "hetero_atoms_proxy": int(sc["hetero_atoms_proxy"]),
                        "total_proxy": round(total,4), "score_source": note})
                write_csv(p / "scores" / f"round_{rn}_scores.csv", scored,
                    ["round","id","smiles","parent","source","validity_proxy","qed_proxy","sa_proxy","lipinski_proxy",
                     "docking_proxy","pose_proxy","novelty_proxy","raw_docking_kcal_mol","heavy_atoms_proxy","hetero_atoms_proxy","total_proxy","score_source"])
                emit(run_id, "step", {"round": rn, "step": "score", "message": f"评分完成", "count": len(scored)})
                # rank
                emit(run_id, "step", {"round": rn, "step": "rank", "message": "排序筛选..."})
                thresholds = cfg.get("thresholds", {})
                advance_total = float(thresholds.get("advance_total", 0.64)) if isinstance(thresholds, dict) else 0.64
                pose_min = float(thresholds.get("pose_min", 0.50)) if isinstance(thresholds, dict) else 0.50
                rows = sorted(scored, key=lambda r2: (float(r2.get("total_proxy",0)), float(r2.get("docking_proxy",0)), float(r2.get("pose_proxy",0))), reverse=True)
                ranked = []; advanced_count = 0
                for idx, row in enumerate(rows, 1):
                    tv, pv = float(row.get("total_proxy",0)), float(row.get("pose_proxy",0))
                    decision = "advance" if idx <= req.top and tv >= advance_total and pv >= pose_min else "hold"
                    if decision == "advance": advanced_count += 1
                    row["rank"] = idx; row["decision"] = decision; ranked.append(row)
                fields = ["rank"] + list(scored[0].keys()) + ["decision"]
                write_csv(p / "ranked" / f"round_{rn}_ranked.csv", ranked, fields)
                summary = render_round_summary(rn, ranked, [r for r in ranked if r["decision"]=="advance"], advance_total, pose_min)
                write_text(p / "reports" / f"round_{rn}_summary.md", summary)
                emit(run_id, "step", {"round": rn, "step": "rank", "message": f"排序完成，{advanced_count} 个晋级", "advanced": advanced_count, "total": len(ranked)})
                # feedback
                emit(run_id, "step", {"round": rn, "step": "feedback", "message": "生成反馈种子..."})
                selected = [r for r in ranked if r["decision"]=="advance"][:req.top] or ranked[:req.top]
                seed_rows = [{"id": r["id"], "smiles": r["smiles"], "note": f"round_{rn}_rank_{r['rank']}_total_{r.get('total_proxy','')}"} for r in selected]
                write_csv(p / "seeds" / f"round_{rn}_seeds.csv", seed_rows, ["id","smiles","note"])
                fb = {"round": rn, "selected_count": len(seed_rows), "next_seed_file": str(p / "seeds" / f"round_{rn}_seeds.csv"), "selected_smiles": [r["smiles"] for r in seed_rows]}
                write_json(p / "feedback" / f"round_{rn}_feedback.json", fb)
                emit(run_id, "step", {"round": rn, "step": "feedback", "message": f"反馈完成，{len(seed_rows)} 个种子"})
                emit(run_id, "round_end", {"round": rn, "advanced": advanced_count, "total": len(ranked)})
            finish_run(run_id)
        except Exception as e:
            emit(run_id, "error", {"message": str(e), "traceback": traceback.format_exc()})
            finish_run(run_id, str(e))
    threading.Thread(target=_run, daemon=True).start()
    return {"run_id": run_id}

@app.get("/api/runs/{run_id}/stream")
async def api_run_stream(run_id: str):
    async def event_stream():
        with run_lock: q = run_queues.get(run_id)
        if not q:
            yield f"data: {json.dumps({'event': 'error', 'data': {'message': 'Run not found'}})}\n\n"; return
        while True:
            try:
                msg = q.get(timeout=30)
                yield f"data: {json.dumps(msg)}\n\n"
                if msg["event"] == "done": break
            except queue.Empty:
                yield f"data: {json.dumps({'event': 'heartbeat'})}\n\n"
    return StreamingResponse(event_stream(), media_type="text/event-stream")

@app.get("/api/runs/{run_id}/status")
def api_run_status(run_id: str):
    with run_lock: s = run_status.get(run_id)
    if not s: raise HTTPException(404, "Run not found")
    return s

# ── 诊断 ──────────────────────────────────────────────────────────────────
@app.get("/api/system/health")
def api_system_health(default_project: str = ""):
    return product_system_health_payload(default_project=default_project)

@app.get("/api/product/bootstrap-health")
def api_product_bootstrap_health(project: str = "", round_no: int = Query(1, alias="round")):
    safe_round = max(1, int(round_no))
    return product_bootstrap_health_payload(project=project, round_no=safe_round)

@app.get("/api/health")
def api_health():
    projects = list_projects()
    names = [str(p.get("name", "")) for p in projects]
    default_project = names[0] if names else ""
    default_dir = PROJECTS_ROOT / default_project if default_project else None
    return {
        "status": "ok",
        "service": "ai-molecule-screening-web",
        "url": "http://localhost:8765/",
        "default_project": default_project,
        "projects_count": len(projects),
        "checks": {
            "index_html": (STATIC / "index.html").exists(),
            "projects_root": PROJECTS_ROOT.exists(),
            "default_project": bool(default_project and default_dir and default_dir.exists()),
            "stage8_command_center": True,
        },
        "boundary": "Computational screening and validation planning only; no efficacy, potency, clinical, dosing, toxicity, or safety claim.",
    }

@app.get("/api/doctor")
def api_doctor():
    config = default_config()
    repos = {}
    for name, raw_path in config.get("tool_paths", {}).items():
        path = Path(str(raw_path)).expanduser()
        repos[name] = {"path": str(path), "present": path.exists()}
    capabilities = stage4_capabilities()
    stage4_executables = capabilities.get("executables", {}) if isinstance(capabilities.get("executables"), dict) else {}
    executables: Dict[str, Any] = {}
    for exe in ["reinvent", "drugex", "vina", "gnina", "bust", "openfe"]:
        status = stage4_executables.get(exe)
        if isinstance(status, dict):
            executables[exe] = status
        else:
            import shutil
            found = shutil.which(exe)
            executables[exe] = {"name": exe, "status": "found" if found else "not_found", "path": found or "", "source": "PATH" if found else ""}
    real_tool_ready = any(isinstance(v, dict) and v.get("status") == "found" for v in executables.values())
    return {
        "repos": repos,
        "executables": executables,
        "capabilities": capabilities,
        "proxy_only": not real_tool_ready,
        "boundary": "System doctor reports local computational tools only; availability is not validation of scientific results.",
    }

@app.get("/api/stats")
def api_stats():
    ps = list_projects()
    total_rounds = 0; total_candidates = 0
    for proj in ps:
        p = PROJECTS_ROOT / proj["name"]
        for f in (p / "candidates").glob("round_*.csv"):
            total_rounds += 1
            try: total_candidates += len(read_csv(f))
            except Exception: pass
    return {"projects": len(ps), "total_rounds": total_rounds, "total_candidates": total_candidates,
            "repos": sum(1 for v in default_config()["tool_paths"].values() if Path(str(v)).expanduser().exists())}

# ── 静态文件 ──────────────────────────────────────────────────────────────
@app.get("/")
def serve_index():
    return FileResponse(
        STATIC / "index.html",
        headers={
            "Cache-Control": "no-store, no-cache, must-revalidate, max-age=0",
            "Pragma": "no-cache",
        },
    )

@app.get("/favicon.ico")
def serve_favicon():
    return Response(status_code=204)

if __name__ == "__main__":
    print("  🧪 AI 分子设计闭环可视化平台 v1.2")
    print(f"  🌐 http://localhost:8765")
    uvicorn.run(app, host="0.0.0.0", port=8765, log_level="info")
