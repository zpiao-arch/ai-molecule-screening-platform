from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import sys
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from . import dataset
from .contracts import BenchmarkConfig
from . import layer_cache, metrics, provenance, verify_run


_L2_INSTANCES: dict[tuple[str, str], object] = {}


def _import_production(scoring_dir: Path):
    path = str(scoring_dir.resolve())
    if path not in sys.path:
        sys.path.insert(0, path)
    import scoring
    try:
        from scoring.l2_bindingdb import BindingDBFeature, Layer2BindingDB
    except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
        from l2_bindingdb import BindingDBFeature, Layer2BindingDB

    return scoring, BindingDBFeature, Layer2BindingDB


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_parquet(temporary, index=False, engine="pyarrow", compression="zstd")
    os.replace(temporary, path)


def _atomic_json(payload: dict, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    os.replace(temporary, path)


def _source_fingerprint(*paths: Path) -> dict[str, str]:
    return {str(path.resolve()): _sha256(path) for path in paths}


def _pool_manifest_fingerprint(target_manifest: pd.DataFrame) -> str:
    values = target_manifest.sort_values("target_index")["pool_sha256"].astype(str).tolist()
    return provenance.payload_sha256(values)


def _begin_stage(
    root: Path,
    stage: str,
    *,
    inputs: dict[str, object],
    resume: bool,
    existing_output: bool,
) -> None:
    checkpoint = provenance.validate_stage_checkpoint(root, stage, inputs=inputs)
    if not checkpoint["ok"] and checkpoint.get("reason") not in {"missing", "ok"}:
        raise RuntimeError(f"{stage} checkpoint {checkpoint['reason']}")
    if checkpoint.get("reason") == "missing" and resume and existing_output:
        raise RuntimeError(f"{stage} checkpoint missing; refusing unverified resume")
    provenance.write_stage_checkpoint(root, stage, inputs=inputs, status="in_progress")


def _complete_stage(root: Path, stage: str, *, inputs: dict[str, object], outputs: dict[str, object]) -> None:
    provenance.write_stage_checkpoint(
        root,
        stage,
        inputs=inputs,
        status="complete",
        outputs=outputs,
    )


def materialize_prepared_run(
    *,
    run_dir: str | Path,
    split: dataset.SplitResult,
    targets: list[dict[str, object]],
    library: pd.DataFrame,
    pool_size: int,
    seed: int,
) -> dict[str, object]:
    """Archive a deterministic target manifest and one independent pool per target."""
    root = Path(run_dir)
    root.mkdir(parents=True, exist_ok=True)
    pool_dir = root / "pool_manifest"
    pool_dir.mkdir(parents=True, exist_ok=True)
    if len({str(row["target_id"]) for row in targets}) != len(targets):
        raise ValueError("target IDs must be unique")

    library_smiles = set(library["canonical_smiles"].dropna().astype(str))
    manifest_rows = []
    total_pairs = 0
    total_fit_overlap = 0
    for index, target in enumerate(targets):
        target_id = str(target["target_id"])
        pool = dataset.build_target_pool(
            target_id=target_id,
            target_index=index,
            fit_rows=split.fit_rows,
            heldout_rows=split.eligible_heldout_rows,
            library=library,
            pool_size=pool_size,
            seed=seed,
        )
        fit_smiles = set(
            split.fit_rows.loc[split.fit_rows["target_id"] == target_id, "canonical_smiles"]
            .dropna()
            .astype(str)
        )
        fit_overlap = len(set(pool["canonical_smiles"]) & fit_smiles)
        if fit_overlap:
            raise RuntimeError(f"{target_id} pool overlaps {fit_overlap} fitted molecules")
        labeled = split.eligible_heldout_rows.loc[
            split.eligible_heldout_rows["target_id"] == target_id
        ]
        label_not_in_library = int((~labeled["canonical_smiles"].isin(library_smiles)).sum())

        pool_path = pool_dir / f"{target_id}.parquet"
        _atomic_parquet(pool, pool_path)
        counts = pool["label_role"].value_counts()
        manifest_row = dict(target)
        manifest_row.update(
            {
                "target_index": index,
                "pool_size": int(len(pool)),
                "pool_file": pool_path.relative_to(root).as_posix(),
                "pool_sha256": _sha256(pool_path),
                "n_pool_positive": int(counts.get("heldout_positive", 0)),
                "n_pool_negative": int(counts.get("heldout_negative", 0)),
                "n_pool_background": int(counts.get("unlabeled_background", 0)),
                "label_not_in_library": label_not_in_library,
                "fit_pair_overlap_count": fit_overlap,
            }
        )
        manifest_rows.append(manifest_row)
        total_pairs += len(pool)
        total_fit_overlap += fit_overlap

    target_manifest = pd.DataFrame(manifest_rows).sort_values("target_index")
    _atomic_parquet(target_manifest, root / "target_manifest.parquet")
    _atomic_parquet(split.fit_rows, root / "reconstructed_fit_rows.parquet")
    _atomic_parquet(split.heldout_rows, root / "reconstructed_heldout_rows.parquet")
    _atomic_parquet(split.eligible_heldout_rows, root / "eligible_heldout_rows.parquet")
    _atomic_parquet(
        library.drop_duplicates("canonical_smiles", keep="first").reset_index(drop=True),
        root / "library_manifest.parquet",
    )

    summary = {
        "n_targets": len(target_manifest),
        "pool_size": int(pool_size),
        "n_pairs": int(total_pairs),
        "fit_pair_overlap_count": int(total_fit_overlap),
        "exact_split_key_overlap_count": int(split.exact_key_overlap_count),
        "eligible_fit_pair_overlap_count": int(split.eligible_fit_pair_overlap_count),
        "seed": int(seed),
    }
    _atomic_json(summary, root / "prepare_summary.json")
    return summary


def _details(value) -> dict:
    if isinstance(value, dict):
        return value
    if isinstance(value, str) and value:
        parsed = json.loads(value)
        return parsed if isinstance(parsed, dict) else {}
    return {}


def _l2_instance(Layer2BindingDB, model_path: str | Path | None):
    key = (str(Path(model_path).resolve()) if model_path else "default", "mlp")
    if key not in _L2_INSTANCES:
        instance = Layer2BindingDB(
            model_path=str(model_path) if model_path else None,
            prefer="mlp",
        )
        instance._ensure_model()
        _L2_INSTANCES[key] = instance
    return _L2_INSTANCES[key]


def _factorized_l2(model, mol_features: np.ndarray, target_features: np.ndarray) -> np.ndarray:
    if not hasattr(model, "coefs_") or len(model.coefs_) != 3:
        combined = np.concatenate(
            [mol_features, np.repeat(target_features[None, :], len(mol_features), axis=0)],
            axis=1,
        )
        return np.asarray(model.predict_proba(combined)[:, 1], dtype=np.float64)

    weights = model.coefs_
    biases = model.intercepts_
    first_molecule = weights[0][: mol_features.shape[1], :]
    first_target = weights[0][mol_features.shape[1] :, :]
    hidden1 = mol_features @ first_molecule + target_features @ first_target + biases[0]
    hidden1 = np.maximum(hidden1, 0.0)
    hidden2 = np.maximum(hidden1 @ weights[1] + biases[1], 0.0)
    logits = (hidden2 @ weights[2] + biases[2]).reshape(-1)
    return 1.0 / (1.0 + np.exp(-logits))


def score_target_frame(
    cache: pd.DataFrame,
    *,
    target_id: str,
    target_text: str,
    scoring_dir: str | Path,
    l2_model_path: str | Path | None = None,
    strict: bool = True,
) -> pd.DataFrame:
    """Score one target pool using production-equivalent factorized L2 inference."""
    required = {"smiles", "mol_features", "l1", "l3", "l4"}
    missing = required - set(cache.columns)
    if missing:
        raise ValueError(f"cache missing required columns: {sorted(missing)}")
    if strict:
        for status_column in ("l1_status", "l3_status", "l4_status"):
            if status_column in cache and not cache[status_column].eq("ok").all():
                raise RuntimeError(f"strict scoring rejected non-ok {status_column}")

    scoring, BindingDBFeature, Layer2BindingDB = _import_production(Path(scoring_dir))
    l2 = _l2_instance(Layer2BindingDB, l2_model_path)
    mol_features = np.asarray(cache["mol_features"].tolist(), dtype=np.float32)
    if mol_features.ndim != 2 or mol_features.shape[1] != BindingDBFeature.FP_SIZE + BindingDBFeature.N_DESC:
        raise ValueError(f"unexpected molecular feature shape: {mol_features.shape}")
    target_features = BindingDBFeature().target_features(target_text)
    probabilities = _factorized_l2(l2._model, mol_features, target_features)
    normalized = np.asarray([l2.normalize(value) for value in probabilities], dtype=np.float64)

    result = cache.copy().reset_index(drop=True)
    result.insert(0, "target_text", target_text)
    result.insert(0, "target_id", target_id)
    result["l2_raw_probability"] = probabilities
    result["l2"] = normalized
    result["layer2_status"] = "ok"
    result["layer2_backend"] = f"BindingDB-L2-{l2.model_kind}-factorized"
    result["layer2_model_asset_id"] = l2.model_path.name
    result["final_score"] = scoring.Layer4Aggregator.combine(
        result["l1"].to_numpy(),
        result["l2"].to_numpy(),
        result["l3"].to_numpy(),
        result["l4"].to_numpy(),
    )

    gates = []
    reasons = []
    detail_values = result["l3_details"] if "l3_details" in result else [{}] * len(result)
    for l2_score, l3_detail in zip(result["l2"], detail_values):
        gate, reason = scoring.Layer4Aggregator.quality_gate(
            {"docking_normalized": l2_score, **_details(l3_detail)}
        )
        gates.append(gate)
        reasons.append(reason)
    result["gate_status"] = gates
    result["gate_reason"] = reasons
    return result


def _load_target_metadata(config: BenchmarkConfig) -> tuple[list[str], dict[str, dict[str, str]]]:
    metadata: dict[str, dict[str, str]] = {}
    with config.chembl_uniprot_mapping.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            uniprot, target_id = fields[:2]
            name = fields[2] if len(fields) > 2 else target_id
            target_type = fields[3] if len(fields) > 3 else ""
            metadata.setdefault(
                target_id,
                {"uniprot_id": uniprot, "target_name": name, "target_type": target_type},
            )
    return sorted(metadata), metadata


def _load_frozen_library(config: BenchmarkConfig) -> pd.DataFrame:
    source = config.frozen_library_smiles
    if not source.is_file():
        raise FileNotFoundError(
            f"frozen 100k library is missing: {source}; rebuild it with build_bigrun_library.py"
        )
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    rows = []
    seen = set()
    for index, raw in enumerate(source.read_text(encoding="utf-8").splitlines()):
        raw = raw.strip()
        if not raw:
            continue
        mol = Chem.MolFromSmiles(raw)
        if mol is None:
            continue
        canonical = Chem.MolToSmiles(mol, canonical=True)
        if canonical in seen:
            continue
        seen.add(canonical)
        rows.append({"library_index": index, "source_smiles": raw, "canonical_smiles": canonical})
    library = pd.DataFrame(rows)
    if len(library) < 10000:
        raise RuntimeError(f"frozen library has only {len(library)} unique valid molecules")
    return library


def prepare_run(
    config: BenchmarkConfig,
    run_dir: str | Path,
    *,
    n_targets: int,
    pool_size: int,
    seed: int,
) -> dict[str, object]:
    root = Path(run_dir)
    inputs = {
        "n_targets": int(n_targets),
        "pool_size": int(pool_size),
        "seed": int(seed),
        "bindingdb_sha256": _sha256(config.bindingdb_examples),
        "chembl_sha256": _sha256(config.chembl_examples),
        "frozen_library_sha256": _sha256(config.frozen_library_smiles),
        "source_code": _source_fingerprint(Path(__file__), Path(dataset.__file__), Path(config.scoring_dir / "pipeline_router.py")),
    }
    _begin_stage(root, "prepare", inputs=inputs, resume=False, existing_output=(root / "target_manifest.parquet").is_file())
    split = dataset.reconstruct_combined_split(config)
    fallback_ids, target_metadata = _load_target_metadata(config)
    targets = dataset.select_targets(
        split.eligible_heldout_rows,
        n_targets=n_targets,
        fallback_target_ids=fallback_ids,
    )

    scoring, _, _ = _import_production(config.scoring_dir)
    del scoring
    try:
        from scoring import pipeline_router
    except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
        import pipeline_router

    enriched = []
    for target in targets:
        row = dict(target)
        target_id = str(row["target_id"])
        metadata = target_metadata.get(target_id, {})
        if row["evaluation_tier"] == "unlabeled_fallback":
            row["target_text"] = " ".join(
                value for value in (target_id, metadata.get("target_name", ""), metadata.get("target_type", "")) if value
            )
        receptor = pipeline_router.lookup_receptor(str(row["target_text"]), chembl_id=target_id)
        receptor_ok, receptor_reason = pipeline_router.validate_receptor(receptor)
        row.update(
            {
                **metadata,
                "sequence_available": bool(metadata.get("uniprot_id")),
                "route_branch": "cascade" if receptor_ok else "library",
                "receptor_available": receptor_ok,
                "receptor": receptor.get("pdbqt") if receptor_ok else None,
                "box_center": receptor.get("box_center") if receptor_ok else None,
                "box_size": receptor.get("box_size") if receptor_ok else None,
                "route_reason": "valid_registered_receptor" if receptor_ok else receptor_reason,
            }
        )
        enriched.append(row)

    library = _load_frozen_library(config)
    summary = materialize_prepared_run(
        run_dir=root,
        split=split,
        targets=enriched,
        library=library,
        pool_size=pool_size,
        seed=seed,
    )
    run_config = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_targets": n_targets,
        "pool_size": pool_size,
        "seed": seed,
        "scoring_dir": str(config.scoring_dir),
        "bindingdb_sha256": _sha256(config.bindingdb_examples),
        "chembl_sha256": _sha256(config.chembl_examples),
        "frozen_library_sha256": _sha256(config.frozen_library_smiles),
        "design": "docs/2026-07-17-four-level-cli-1000x10000-design.md",
    }
    _atomic_json(run_config, root / "run_config.json")
    _complete_stage(
        root,
        "prepare",
        inputs=inputs,
        outputs={
            "target_manifest_sha256": _sha256(root / "target_manifest.parquet"),
            "pool_manifest_fingerprint": _pool_manifest_fingerprint(pd.read_parquet(root / "target_manifest.parquet")),
            "n_pairs": int(summary["n_pairs"]),
        },
    )
    return summary


def cache_run_layers(
    config: BenchmarkConfig,
    run_dir: str | Path,
    *,
    batch_size: int,
    shard_size: int,
    resume: bool,
    strict: bool,
    multi_process: bool = False,
    unimol_device: str = "mps",
) -> dict[str, object]:
    root = Path(run_dir)
    library = pd.read_parquet(root / "library_manifest.parquet")
    namespace, _ = layer_cache.cache_namespace(config.scoring_dir)
    inputs = {
        "library_manifest_sha256": _sha256(root / "library_manifest.parquet"),
        "n_molecules": int(len(library)),
        "batch_size": int(batch_size),
        "shard_size": int(shard_size),
        "multi_process": bool(multi_process),
        "unimol_device": str(unimol_device),
        "cache_namespace": namespace,
    }
    existing_output = bool(list((root / "layer_cache" / "layer_cache_shards").glob("part-*.parquet")))
    _begin_stage(root, "layer_cache", inputs=inputs, resume=resume, existing_output=existing_output)
    cache = layer_cache.build(
        library["canonical_smiles"].tolist(),
        root / "layer_cache",
        batch_size=batch_size,
        shard_size=shard_size,
        resume=resume,
        strict=strict,
        multi_process=multi_process,
        unimol_device=unimol_device,
        scoring_dir=config.scoring_dir,
    )
    summary = {
        "n_molecules": int(len(cache)),
        "l1_failures": int((cache["l1_status"] != "ok").sum()),
        "l3_failures": int((cache["l3_status"] != "ok").sum()),
        "l4_failures": int((cache["l4_status"] != "ok").sum()),
        "l4_std": float(pd.to_numeric(cache["l4"], errors="coerce").std()),
    }
    if strict and (summary["l1_failures"] or summary["l3_failures"] or summary["l4_failures"]):
        raise RuntimeError(f"strict layer cache contains failures: {summary}")
    _atomic_json(summary, root / "layer_cache_summary.json")
    _complete_stage(
        root,
        "layer_cache",
        inputs=inputs,
        outputs={
            "summary_sha256": _sha256(root / "layer_cache_summary.json"),
            "cache_metadata_sha256": _sha256(root / "layer_cache" / "cache_metadata.json"),
            "n_molecules": int(len(cache)),
        },
    )
    return summary


def _load_layer_cache(run_dir: Path) -> pd.DataFrame:
    paths = sorted((run_dir / "layer_cache" / "layer_cache_shards").glob("part-*.parquet"))
    if not paths:
        raise FileNotFoundError("no layer cache shards found")
    cache = pd.concat((pd.read_parquet(path) for path in paths), ignore_index=True)
    if not cache["smiles"].is_unique:
        raise RuntimeError("layer cache contains duplicate SMILES")
    return cache.set_index("smiles", drop=False)


def _candidate_evidence(scored: pd.DataFrame, pool: pd.DataFrame) -> pd.DataFrame:
    output = pd.DataFrame(
        {
            "target_id": scored["target_id"],
            "target_text": scored["target_text"],
            "candidate_index": np.arange(len(scored), dtype=np.int32),
            "canonical_smiles": scored["smiles"],
            "label_role": pool["label_role"].to_numpy(),
            "label": pool["label"].to_numpy(),
            "l1": scored["l1"],
            "l2": scored["l2"],
            "l2_raw_probability": scored["l2_raw_probability"],
            "l3": scored["l3"],
            "l4": scored["l4"],
            "final_score": scored["final_score"],
            "gate_status": scored["gate_status"],
            "gate_reason": scored["gate_reason"],
            "l1_status": scored["l1_status"],
            "layer2_status": scored["layer2_status"],
            "l3_status": scored["l3_status"],
            "l4_status": scored["l4_status"],
            "l1_backend": scored["l1_backend"],
            "layer2_backend": scored["layer2_backend"],
            "l3_backend": scored["l3_backend"],
            "l4_backend": scored["l4_backend"],
            "l1_model_asset_id": scored["l1_model_asset_id"],
            "layer2_model_asset_id": scored["layer2_model_asset_id"],
            "l3_model_asset_id": scored["l3_model_asset_id"],
            "l4_model_asset_id": scored["l4_model_asset_id"],
            "l4_pos_similarity": scored["l4_pos_similarity"],
        }
    )
    l1_details = [_details(value) for value in scored["l1_details"]]
    l3_details = [_details(value) for value in scored["l3_details"]]
    for column in ("mw", "logp", "tpsa", "hbd", "hba", "qed", "sa", "lipinski_violations"):
        output[column] = [detail.get(column) for detail in l1_details]
    for column in (
        "toxicity_count", "toxicity_severe", "toxicity_flags", "solubility_logS",
        "bbb_prob", "cyp3a4_risk", "cyp2d6_risk", "herg_risk",
        "hepatotoxicity_risk", "oral_bioavailability",
    ):
        output[column] = [detail.get(column) for detail in l3_details]
    return output


_PARITY_NUMERIC_FIELDS = ("l1", "l2", "l3", "l4", "final_score")
_PARITY_CATEGORICAL_FIELDS = (
    "gate_status",
    "l1_status", "layer2_status", "l3_status", "l4_status",
    "l1_model_asset_id", "layer2_model_asset_id", "l3_model_asset_id", "l4_model_asset_id",
)


def compare_parity_outputs(
    batch: pd.DataFrame,
    scalar: pd.DataFrame,
    *,
    tolerance: float = 1e-4,
) -> dict[str, object]:
    if len(batch) != len(scalar):
        return {
            "passed": False,
            "reason": "row_count_mismatch",
            "batch_rows": int(len(batch)),
            "scalar_rows": int(len(scalar)),
        }
    if batch["canonical_smiles"].astype(str).tolist() != scalar["canonical_smiles"].astype(str).tolist():
        return {"passed": False, "reason": "smiles_order_mismatch"}

    numeric = {}
    for field in _PARITY_NUMERIC_FIELDS:
        left = pd.to_numeric(batch[field], errors="coerce").to_numpy(dtype=float)
        right = pd.to_numeric(scalar[field], errors="coerce").to_numpy(dtype=float)
        finite_match = bool(np.array_equal(np.isfinite(left), np.isfinite(right)))
        difference = np.abs(left - right)
        max_abs_diff = float(np.nanmax(difference)) if difference.size else 0.0
        numeric[field] = {
            "max_abs_diff": max_abs_diff,
            "finite_mask_match": finite_match,
            "passed": finite_match and max_abs_diff <= tolerance,
        }

    categorical = {}
    for field in _PARITY_CATEGORICAL_FIELDS:
        left = batch[field].fillna("<null>").astype(str).tolist()
        right = scalar[field].fillna("<null>").astype(str).tolist()
        mismatches = [index for index, (a, b) in enumerate(zip(left, right)) if a != b]
        categorical[field] = {"mismatch_indices": mismatches, "passed": not mismatches}

    passed = all(item["passed"] for item in numeric.values()) and all(
        item["passed"] for item in categorical.values()
    )
    return {
        "passed": bool(passed),
        "n_molecules": int(len(batch)),
        "tolerance": float(tolerance),
        "numeric_fields": numeric,
        "categorical_fields": categorical,
    }


def _select_parity_fixture(
    cache: pd.DataFrame,
    run_dir: Path,
    target: pd.Series,
    *,
    n_molecules: int = 20,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Select a stable, varied valid fixture without using evaluation scores."""
    frame = cache.reset_index(drop=True)
    selected: list[int] = []
    categories: dict[str, str] = {}

    def add(index: int, category: str) -> None:
        if index not in selected:
            selected.append(index)
            categories[str(frame.loc[index, "smiles"])] = category

    # Include held-out positives for the one registered receptor when available.
    na_pool_path = run_dir / "pool_manifest" / "CHEMBL2051.parquet"
    if na_pool_path.is_file():
        na_pool = pd.read_parquet(na_pool_path)
        for smi in na_pool.loc[
            na_pool["label_role"].eq("heldout_positive"), "canonical_smiles"
        ].astype(str).head(2):
            matches = frame.index[frame["smiles"].astype(str).eq(smi)].tolist()
            if matches:
                add(matches[0], "CHEMBL2051_heldout_positive")

    details = [_details(value) for value in frame["l1_details"]]
    mw_values = [float(item.get("mw", 0.0) or 0.0) for item in details]
    if mw_values:
        add(int(np.argmax(mw_values)), "high_molecular_weight")
    l3_details = [_details(value) for value in frame["l3_details"]]
    toxic = [index for index, item in enumerate(l3_details) if int(item.get("toxicity_count", 0) or 0) > 0]
    if toxic:
        add(toxic[0], "toxicity_flagged")

    for index in range(len(frame)):
        if len(selected) >= n_molecules:
            break
        add(index, "coverage_normal")
    fixture = frame.iloc[selected[:n_molecules]].reset_index(drop=True)
    return fixture, categories


def _scalar_parity_frame(
    scoring,
    fixture: pd.DataFrame,
    *,
    target_text: str,
    l2_model_path: Path,
) -> pd.DataFrame:
    scorer = scoring.MoleculeScorer(
        use_unimol=True,
        l2_method="bindingdb",
        default_target_text=target_text,
        l2_model_path=str(l2_model_path),
        strict_backends=True,
        unimol_device="mps",
    )
    molecules = [(f"parity_{index}", smiles) for index, smiles in enumerate(fixture["smiles"].astype(str))]
    scored = scorer.score_batch(molecules, target_text=target_text)
    by_smiles = {str(row["smiles"]): row for row in scored}
    rows = []
    for smiles in fixture["smiles"].astype(str):
        result = by_smiles[smiles]
        rows.append(
            {
                "canonical_smiles": smiles,
                "l1": result["layer1_score"],
                "l2": result["docking_normalized"],
                "l3": result["admet_score"],
                "l4": result["unimol_score"],
                "final_score": result["final_score"],
                "gate_status": result["gate_status"],
                "l1_status": result["layer1_status"],
                "layer2_status": result["layer2_status"],
                "l3_status": result["layer3_status"],
                "l4_status": result["layer4_status"],
                "l1_model_asset_id": result["layer1_model_asset_id"],
                "layer2_model_asset_id": result["layer2_model_asset_id"],
                "l3_model_asset_id": result["layer3_model_asset_id"],
                "l4_model_asset_id": result["layer4_model_asset_id"],
            }
        )
    return pd.DataFrame(rows)


def _parity_gate(
    config: BenchmarkConfig,
    cache: pd.DataFrame,
    target: pd.Series,
    *,
    batch_model_path: str | Path | None = None,
    run_dir: str | Path | None = None,
) -> dict[str, object]:
    root = Path(run_dir) if run_dir is not None else Path.cwd()
    fixture, categories = _select_parity_fixture(cache, root, target)
    batch = score_target_frame(
        fixture,
        target_id=str(target["target_id"]),
        target_text=str(target["target_text"]),
        scoring_dir=config.scoring_dir,
        l2_model_path=batch_model_path,
    )
    scoring, _, _ = _import_production(config.scoring_dir)
    scalar = _scalar_parity_frame(
        scoring,
        fixture,
        target_text=str(target["target_text"]),
        l2_model_path=Path(batch_model_path or config.l2_model_path),
    )
    batch_canonical = batch.rename(columns={"smiles": "canonical_smiles"})
    comparison = compare_parity_outputs(batch_canonical, scalar, tolerance=1e-4)

    invalid_smiles = "not-a-valid-smiles"
    try:
        scoring.MoleculeScorer(
            use_unimol=True,
            l2_method="bindingdb",
            default_target_text=str(target["target_text"]),
            l2_model_path=str(config.l2_model_path),
            strict_backends=True,
        ).score_one(invalid_smiles, target_text=str(target["target_text"]))
        scalar_rejected = False
    except Exception:
        scalar_rejected = True
    try:
        layer_cache._canonicalize(invalid_smiles)
        batch_rejected = False
    except Exception:
        batch_rejected = True
    invalid_input = {
        "smiles": invalid_smiles,
        "scalar_rejected": scalar_rejected,
        "batch_rejected": batch_rejected,
        "passed": scalar_rejected and batch_rejected,
    }

    route_branch_match = True
    try:
        if str(target["target_id"]) == "CHEMBL2051":
            try:
                from scoring import pipeline_router
            except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
                import pipeline_router
            decision = pipeline_router.route("CHEMBL2051", chembl_id="CHEMBL2051")
            route_branch_match = decision.get("branch") == str(target["route_branch"])
    except Exception:
        route_branch_match = False

    evidence = {
        "n_molecules": int(len(fixture)),
        "tolerance": 1e-4,
        "numeric_fields": comparison.get("numeric_fields", {}),
        "categorical_fields": comparison.get("categorical_fields", {}),
        "passed": bool(comparison.get("passed") and route_branch_match and invalid_input["passed"]),
        "fixture_categories": categories,
        "route_branch_match": route_branch_match,
        "invalid_input": invalid_input,
    }
    if run_dir is not None:
        parity_rows = batch_canonical.copy()
        for column in scalar.columns:
            if column == "canonical_smiles":
                continue
            parity_rows[f"scalar_{column}"] = scalar[column].to_numpy()
        _atomic_parquet(parity_rows, Path(run_dir) / "parity_rows.parquet")
    return evidence


def score_run(
    config: BenchmarkConfig,
    run_dir: str | Path,
    *,
    resume: bool,
    limit_targets: int | None = None,
) -> dict[str, object]:
    root = Path(run_dir)
    targets = pd.read_parquet(root / "target_manifest.parquet").sort_values("target_index")
    if limit_targets is not None:
        targets = targets.head(limit_targets)
    cache = _load_layer_cache(root)
    target_manifest_path = root / "target_manifest.parquet"
    inputs = {
        "target_manifest_sha256": _sha256(target_manifest_path),
        "pool_manifest_fingerprint": _pool_manifest_fingerprint(targets),
        "cache_metadata_sha256": _sha256(root / "layer_cache" / "cache_metadata.json"),
        "l2_model_sha256": _sha256(config.l2_model_path),
        "source_code": _source_fingerprint(Path(__file__), Path(config.scoring_dir / "scoring.py"), Path(config.scoring_dir / "l2_bindingdb.py")),
    }
    existing_output = bool(list((root / "scores").glob("target_id=*/part-*.parquet")))
    _begin_stage(root, "score", inputs=inputs, resume=resume, existing_output=existing_output)
    parity = _parity_gate(
        config,
        cache.reset_index(drop=True),
        targets.iloc[0],
        batch_model_path=config.l2_model_path,
        run_dir=root,
    )
    _atomic_json(parity, root / "parity_evidence.json")
    if not parity["passed"]:
        raise RuntimeError(f"batch parity gate failed: {parity}")

    score_root = root / "scores"
    score_root.mkdir(exist_ok=True)
    per_target = []
    completed_rows = 0
    for completed, target in enumerate(targets.itertuples(index=False), 1):
        target_id = str(target.target_id)
        output_dir = score_root / f"target_id={target_id}"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "part-000000.parquet"
        if resume and output_path.is_file():
            evidence = pd.read_parquet(output_path)
            if len(evidence) != int(target.pool_size):
                raise RuntimeError(f"resume score row mismatch for {target_id}")
        else:
            pool = pd.read_parquet(root / str(target.pool_file)).reset_index(drop=True)
            try:
                cached = cache.loc[pool["canonical_smiles"].tolist()].reset_index(drop=True)
            except KeyError as exc:
                raise RuntimeError(f"pool molecule missing from layer cache for {target_id}") from exc
            scored = score_target_frame(
                cached,
                target_id=target_id,
                target_text=str(target.target_text),
                scoring_dir=config.scoring_dir,
                l2_model_path=config.l2_model_path,
            )
            evidence = _candidate_evidence(scored, pool)
            _atomic_parquet(evidence, output_path)
        target_metrics = metrics.compute_target_metrics(evidence)
        target_metrics.update(
            {
                "target_id": target_id,
                "evaluation_tier": str(target.evaluation_tier),
                "route_branch": str(target.route_branch),
            }
        )
        per_target.append(target_metrics)
        completed_rows += len(evidence)
        _atomic_json(
            {"completed_targets": completed, "completed_rows": completed_rows, "last_target_id": target_id},
            root / "score_progress.json",
        )

    metric_frame = pd.DataFrame(per_target)
    _atomic_parquet(metric_frame, root / "per_target_metrics.parquet")
    summary = {
        "n_targets": int(len(targets)),
        "n_pairs": int(completed_rows),
        "aggregates": metrics.aggregate_metrics(metric_frame, bootstrap_reps=2000, seed=42),
    }
    _atomic_json(summary, root / "score_summary.json")
    _complete_stage(
        root,
        "score",
        inputs=inputs,
        outputs={
            "summary_sha256": _sha256(root / "score_summary.json"),
            "n_targets": int(len(targets)),
            "n_pairs": int(completed_rows),
        },
    )
    return summary


def select_docking_rows(scores: pd.DataFrame, *, top_n: int) -> pd.DataFrame:
    if top_n <= 0:
        return scores.iloc[0:0].copy()
    result = scores.copy()
    if "candidate_index" not in result:
        result.insert(0, "candidate_index", np.arange(len(result), dtype=np.int32))
    result = result.sort_values(
        ["l2", "candidate_index"], ascending=[False, True], kind="mergesort"
    ).head(top_n).copy()
    result["selection_rank"] = np.arange(1, len(result) + 1, dtype=np.int32)
    return result.reset_index(drop=True)


def fuse_docking_results(
    scores: pd.DataFrame,
    docking: pd.DataFrame,
    *,
    weight: float = 0.30,
) -> tuple[pd.DataFrame, dict[str, dict[str, object]]]:
    """Apply LE correction only to successful selected rows and recompute metrics."""
    scoring_dir = Path(os.environ.get(
        "FOUR_LEVEL_SCORING_DIR",
        Path(__file__).resolve().parents[2] / "scoring",
    ))
    if str(scoring_dir) not in sys.path:
        sys.path.insert(0, str(scoring_dir))
    try:
        from scoring.dock_rerank import cascade_corrected_fusion
    except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
        from dock_rerank import cascade_corrected_fusion

    result = scores.copy().reset_index(drop=True)
    result["docking_status"] = "not_selected"
    result["docking_affinity"] = np.nan
    result["ligand_efficiency"] = np.nan
    result["docking_selection_rank"] = np.nan
    by_index = docking.set_index("candidate_index") if not docking.empty else pd.DataFrame()
    le = np.full(len(result), np.nan, dtype=float)
    for row_index, candidate_index in enumerate(result["candidate_index"]):
        if candidate_index not in by_index.index:
            continue
        record = by_index.loc[candidate_index]
        if isinstance(record, pd.DataFrame):
            record = record.iloc[0]
        status = str(record.get("status", "error:missing_status"))
        result.loc[row_index, "docking_status"] = status
        result.loc[row_index, "docking_affinity"] = record.get("affinity")
        result.loc[row_index, "ligand_efficiency"] = record.get("ligand_efficiency")
        result.loc[row_index, "docking_selection_rank"] = record.get("selection_rank")
        if status == "ok" and pd.notna(record.get("ligand_efficiency")):
            le[row_index] = float(record["ligand_efficiency"])
    result["fused_score"] = cascade_corrected_fusion(result["l2"].to_numpy(), le, w=weight)
    before = metrics.compute_target_metrics(result.drop(columns=["fused_score"]))
    after = metrics.compute_target_metrics(result.assign(final_score=result["fused_score"]))
    return result, {"before": before, "after": after, "weight": float(weight)}


def _dock_one_explicit(
    row: pd.Series,
    *,
    receptor: str,
    center: list[float],
    size: list[float],
    smina_bin: str,
    obabel_bin: str,
    workdir: Path,
    timeout: int,
    exhaustiveness: int,
    cpu: int,
    hac_max: int,
) -> dict[str, object]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    candidate_index = int(row["candidate_index"])
    smiles = str(row["canonical_smiles"])
    result: dict[str, object] = {
        "candidate_index": candidate_index,
        "canonical_smiles": smiles,
        "selection_rank": int(row["selection_rank"]),
        "status": "error:unknown",
        "affinity": None,
        "heavy_atoms": None,
        "ligand_efficiency": None,
        "obabel_command": None,
        "smina_command": None,
    }
    try:
        mol = Chem.MolFromSmiles(smiles)
        hac = Descriptors.HeavyAtomCount(mol) if mol is not None else None
        result["heavy_atoms"] = hac
        if hac is None:
            result["status"] = "invalid_smiles"
            return result
        if hac > hac_max:
            result["status"] = "skipped_hac"
            return result
        ligand_path = workdir / f"candidate_{candidate_index}.pdbqt"
        pose_path = workdir / f"candidate_{candidate_index}_pose.pdbqt"
        obabel_command = [
            obabel_bin, f"-:{smiles}", "-O", str(ligand_path), "--gen3d", "-p", "7.4"
        ]
        result["obabel_command"] = json.dumps(obabel_command, ensure_ascii=True)
        try:
            prep = subprocess.run(obabel_command, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            result["status"] = "prep_timeout"
            return result
        if prep.returncode != 0 or not ligand_path.is_file() or ligand_path.stat().st_size == 0:
            result["status"] = "prep_failed"
            return result

        smina_command = [
            smina_bin, "--receptor", receptor, "--ligand", str(ligand_path),
            "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
            "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
            "--exhaustiveness", str(exhaustiveness), "--cpu", str(cpu), "--num_modes", "3",
            "--seed", "42", "--out", str(pose_path),
        ]
        result["smina_command"] = json.dumps(smina_command, ensure_ascii=True)
        try:
            docked = subprocess.run(smina_command, capture_output=True, text=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            result["status"] = "dock_timeout"
            return result
        affinity = None
        for line in docked.stdout.splitlines():
            fields = line.split()
            if len(fields) >= 2 and fields[0] == "1":
                try:
                    affinity = float(fields[1])
                except ValueError:
                    pass
                break
        if affinity is None:
            result["status"] = "dock_no_score"
            return result
        result["affinity"] = round(affinity, 3)
        result["ligand_efficiency"] = round(affinity / hac, 4)
        result["status"] = "ok"
        return result
    except Exception as exc:
        result["status"] = f"error:{type(exc).__name__}"
        return result


def dock_run(
    config: BenchmarkConfig,
    run_dir: str | Path,
    *,
    target_id: str,
    top_n: int,
    resume: bool,
    workers: int,
    timeout: int = 90,
    exhaustiveness: int = 4,
    cpu: int = 2,
    hac_max: int = 50,
) -> dict[str, object]:
    root = Path(run_dir)
    target_manifest = pd.read_parquet(root / "target_manifest.parquet")
    target_rows = target_manifest.loc[target_manifest["target_id"] == target_id]
    if len(target_rows) != 1:
        raise ValueError(f"target {target_id} is not uniquely present in target manifest")
    target = target_rows.iloc[0]
    if not bool(target["receptor_available"]):
        raise RuntimeError(f"target {target_id} has no valid registered receptor")
    score_path = root / "scores" / f"target_id={target_id}" / "part-000000.parquet"
    inputs = {
        "target_id": target_id,
        "score_sha256": _sha256(score_path),
        "receptor_sha256": _sha256(Path(str(target["receptor"]))),
        "top_n": int(top_n),
        "workers": int(workers),
        "timeout": int(timeout),
        "exhaustiveness": int(exhaustiveness),
        "cpu": int(cpu),
        "hac_max": int(hac_max),
        "source_code": _source_fingerprint(Path(__file__), Path(config.scoring_dir / "dock_rerank.py")),
    }
    dock_dir_for_checkpoint = root / "docking" / target_id
    existing_output = (dock_dir_for_checkpoint / "selected_results.parquet").is_file()
    _begin_stage(root, "dock", inputs=inputs, resume=resume, existing_output=existing_output)
    scores = pd.read_parquet(score_path)
    selected = select_docking_rows(scores, top_n=top_n)
    docking_dir = root / "docking" / target_id
    docking_dir.mkdir(parents=True, exist_ok=True)
    result_path = docking_dir / "selected_results.parquet"
    existing = pd.read_parquet(result_path) if resume and result_path.is_file() else pd.DataFrame()
    done = set(existing["candidate_index"].astype(int)) if not existing.empty else set()
    pending = selected.loc[~selected["candidate_index"].isin(done)].copy()

    if pending.empty:
        docking_results = existing
    else:
        if str(os.environ.get("DOCK_BIN_DIR", "")):
            bin_dir = Path(os.environ["DOCK_BIN_DIR"])
            smina_bin, obabel_bin = str(bin_dir / "smina"), str(bin_dir / "obabel")
        else:
            scoring_dir = config.scoring_dir
            if str(scoring_dir) not in sys.path:
                sys.path.insert(0, str(scoring_dir))
            try:
                from scoring.dock_rerank import find_binary
            except ImportError:  # Direct script/PYTHONPATH=scoring compatibility.
                from dock_rerank import find_binary
            smina_bin, obabel_bin = find_binary("smina"), find_binary("obabel")
        if not smina_bin or not obabel_bin:
            raise RuntimeError("smina/obabel binaries are unavailable")
        if not Path(str(target["receptor"])).is_file():
            raise RuntimeError("registered receptor file is unavailable")
        workdir = docking_dir / "work"
        workdir.mkdir(parents=True, exist_ok=True)
        jobs = []
        for _, row in pending.iterrows():
            jobs.append(
                {
                    "row": row,
                    "receptor": str(target["receptor"]),
                    "center": list(target["box_center"]),
                    "size": list(target["box_size"]),
                    "smina_bin": smina_bin,
                    "obabel_bin": obabel_bin,
                    "workdir": workdir,
                    "timeout": timeout,
                    "exhaustiveness": exhaustiveness,
                    "cpu": cpu,
                    "hac_max": hac_max,
                }
            )
        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
            new_results = list(pool.map(lambda job: _dock_one_explicit(**job), jobs))
        docking_results = pd.concat([existing, pd.DataFrame(new_results)], ignore_index=True)
        docking_results = docking_results.drop_duplicates("candidate_index", keep="last").sort_values("candidate_index")
        _atomic_parquet(docking_results, result_path)

    fused, comparison = fuse_docking_results(scores, docking_results)
    fused_path = docking_dir / "fused_scores.parquet"
    _atomic_parquet(fused, fused_path)
    _atomic_json(comparison, docking_dir / "metrics_before_after.json")
    _atomic_json(
        {
            "target_id": target_id,
            "top_n_requested": top_n,
            "n_submitted": int(len(docking_results)),
            "status_counts": docking_results["status"].value_counts().to_dict(),
            "receptor": str(target["receptor"]),
            "receptor_sha256": _sha256(Path(str(target["receptor"]))),
            "smina_bin": smina_bin if pending.size else None,
            "obabel_bin": obabel_bin if pending.size else None,
            "timeout_seconds": timeout,
            "exhaustiveness": exhaustiveness,
            "cpu": cpu,
            "hac_max": hac_max,
        },
        docking_dir / "docking_summary.json",
    )
    _complete_stage(
        root,
        "dock",
        inputs=inputs,
        outputs={
            "summary_sha256": _sha256(docking_dir / "docking_summary.json"),
            "metrics_sha256": _sha256(docking_dir / "metrics_before_after.json"),
            "n_submitted": int(len(docking_results)),
        },
    )
    return {"target_id": target_id, "n_submitted": int(len(docking_results)), "comparison": comparison}


def _default_run_dir(config: BenchmarkConfig) -> Path:
    env = os.environ.get("FOUR_LEVEL_RUN_DIR")
    if env:
        return Path(env).resolve()
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return config.root / "scientific_validation" / "four_level_cli_1kx10k" / "runs" / run_id


def _parse_positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("value must be an integer") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return parsed


def _parse_nonnegative_int(value: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError("value must be an integer") from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return parsed


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Auditable four-level 1000x10000 benchmark CLI")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--run-dir", type=Path)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="reconstruct split and materialize target pools")
    prepare.add_argument("--targets", type=_parse_positive_int, default=1000)
    prepare.add_argument("--pool-size", type=_parse_positive_int, default=10000)
    prepare.add_argument("--seed", type=_parse_nonnegative_int, default=42)

    cache = subparsers.add_parser("cache-layers", help="build strict L1/L3/L4 molecular cache")
    cache.add_argument("--batch-size", type=_parse_positive_int, default=32)
    cache.add_argument("--shard-size", type=_parse_positive_int, default=512)
    cache.add_argument("--resume", action="store_true")
    cache.add_argument("--strict", action="store_true")
    cache.add_argument("--multi-process", action="store_true")
    cache.add_argument("--unimol-device", choices=["cpu", "mps"], default="mps")

    score = subparsers.add_parser("score", help="run 10M factorized L2 and four-level aggregation")
    score.add_argument("--resume", action="store_true")
    score.add_argument("--limit-targets", type=_parse_positive_int)

    dock = subparsers.add_parser("dock", help="dock registered-target L2 top-N candidates")
    dock.add_argument("--target", default="CHEMBL2051")
    dock.add_argument("--top-n", type=_parse_positive_int, default=300)
    dock.add_argument("--resume", action="store_true")
    dock.add_argument("--workers", type=_parse_positive_int, default=4)
    dock.add_argument("--timeout", type=_parse_positive_int, default=90)
    dock.add_argument("--exhaustiveness", type=_parse_positive_int, default=4)
    dock.add_argument("--cpu", type=_parse_positive_int, default=2)
    dock.add_argument("--hac-max", type=_parse_positive_int, default=50)
    subparsers.add_parser("report", help="recompute, render and seal the audit report")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    config = BenchmarkConfig(root=args.root)
    run_dir = (args.run_dir or _default_run_dir(config)).resolve()
    stage = {"prepare": "prepare", "cache-layers": "layer_cache", "score": "score", "dock": "dock"}.get(args.command)
    try:
        if args.command == "prepare":
            result = prepare_run(config, run_dir, n_targets=args.targets, pool_size=args.pool_size, seed=args.seed)
        elif args.command == "cache-layers":
            result = cache_run_layers(
                config, run_dir, batch_size=args.batch_size, shard_size=args.shard_size,
                resume=args.resume, strict=args.strict, multi_process=args.multi_process,
                unimol_device=args.unimol_device,
            )
        elif args.command == "score":
            result = score_run(config, run_dir, resume=args.resume, limit_targets=args.limit_targets)
        elif args.command == "report":
            from . import report

            summary_path, report_path = report.generate_report(run_dir)
            result = {
                "summary": str(summary_path),
                "report": str(report_path),
                "manifest": str(run_dir / verify_run.MANIFEST_NAME),
            }
        else:
            result = dock_run(
                config,
                run_dir,
                target_id=args.target,
                top_n=args.top_n,
                resume=args.resume,
                workers=args.workers,
                timeout=args.timeout,
                exhaustiveness=args.exhaustiveness,
                cpu=args.cpu,
                hac_max=args.hac_max,
            )
    except Exception as exc:
        if stage:
            provenance.fail_stage_checkpoint(run_dir, stage, exc)
        raise
    print(json.dumps({"run_dir": str(run_dir), **result}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
