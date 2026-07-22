from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from . import metrics


MANIFEST_NAME = "MANIFEST.sha256"


def render_canonical_json(payload: object) -> str:
    return json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True)


def verify_summary_bytes(run_dir: str | Path, recomputed: object) -> dict[str, object]:
    path = Path(run_dir) / "summary.json"
    if not path.is_file():
        return {"ok": False, "reason": "summary_missing"}
    expected = render_canonical_json(recomputed)
    actual = path.read_text(encoding="utf-8")
    return {
        "ok": actual == expected,
        "reason": "ok" if actual == expected else "summary_not_byte_identical",
    }


def verify_provenance(run_dir: str | Path) -> dict[str, object]:
    root = Path(run_dir)
    stages = ("prepare", "layer_cache", "score", "dock", "report")
    required = [
        "DESIGN.md",
        "run.log",
        "failure_manifest.json",
        *(f"checkpoints/{stage}.json" for stage in stages),
    ]
    missing = [relative for relative in required if not (root / relative).is_file()]
    invalid: list[str] = []
    checkpoints: dict[str, dict[str, object]] = {}

    design_path = root / "DESIGN.md"
    if design_path.is_file() and not design_path.read_text(encoding="utf-8").strip():
        invalid.append("DESIGN.md")

    failure_path = root / "failure_manifest.json"
    if failure_path.is_file():
        try:
            if not isinstance(json.loads(failure_path.read_text(encoding="utf-8")), dict):
                invalid.append("failure_manifest.json")
        except (OSError, UnicodeError, json.JSONDecodeError):
            invalid.append("failure_manifest.json")

    for stage in stages:
        relative = f"checkpoints/{stage}.json"
        path = root / relative
        if not path.is_file():
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError):
            invalid.append(relative)
            continue
        if (
            not isinstance(payload, dict)
            or payload.get("stage") != stage
            or payload.get("status") != "complete"
        ):
            invalid.append(relative)
            continue
        checkpoints[stage] = payload

    log_path = root / "run.log"
    if log_path.is_file():
        lines = [line for line in log_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        current_event_found = False
        for line in lines:
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            stage = event.get("stage") if isinstance(event, dict) else None
            checkpoint = checkpoints.get(str(stage))
            if not checkpoint:
                continue
            if (
                event.get("event") == "stage_checkpoint"
                and event.get("status") == checkpoint.get("status")
                and event.get("input_fingerprint") == checkpoint.get("input_fingerprint")
                and event.get("checkpoint") == f"checkpoints/{stage}.json"
            ):
                current_event_found = True
                break
        legacy_stage_logs = [
            path
            for path in root.glob("*.log")
            if path.name != "run.log" and path.is_file() and path.stat().st_size > 0
        ]
        if not current_event_found and not (lines and legacy_stage_logs):
            invalid.append("run.log")

    invalid = sorted(set(invalid))
    return {"ok": not missing and not invalid, "missing": missing, "invalid": invalid}


PARITY_NUMERIC_FIELDS = {"l1", "l2", "l3", "l4", "final_score"}
PARITY_CATEGORICAL_FIELDS = {
    "gate_status",
    "l1_status", "layer2_status", "l3_status", "l4_status",
    "l1_model_asset_id", "layer2_model_asset_id", "l3_model_asset_id", "l4_model_asset_id",
}


def verify_parity_evidence(run_dir: str | Path, *, strict: bool = True) -> dict[str, object]:
    path = Path(run_dir) / "parity_evidence.json"
    if not path.is_file():
        return {
            "ok": False,
            "reason": "parity_evidence_missing",
            "missing_numeric_fields": sorted(PARITY_NUMERIC_FIELDS),
            "missing_categorical_fields": sorted(PARITY_CATEGORICAL_FIELDS),
            "missing_checks": ["route_branch_match", "invalid_input"],
        }
    payload = json.loads(path.read_text(encoding="utf-8"))
    numeric = payload.get("numeric_fields", {})
    categorical = payload.get("categorical_fields", {})
    missing_numeric = sorted(PARITY_NUMERIC_FIELDS - set(numeric))
    missing_categorical = sorted(PARITY_CATEGORICAL_FIELDS - set(categorical))
    missing_checks = [name for name in ("route_branch_match", "invalid_input") if name not in payload]
    checks_ok = bool(payload.get("passed"))
    if strict:
        checks_ok = (
            checks_ok
            and int(payload.get("n_molecules", 0)) >= 20
            and not missing_numeric
            and not missing_categorical
            and not missing_checks
            and all(bool(item.get("passed")) for item in numeric.values())
            and all(bool(item.get("passed")) for item in categorical.values())
            and bool(payload.get("route_branch_match"))
            and bool(payload.get("invalid_input", {}).get("passed"))
        )
    return {
        "ok": checks_ok,
        "missing_numeric_fields": missing_numeric,
        "missing_categorical_fields": missing_categorical,
        "missing_checks": missing_checks,
        "evidence": payload,
    }


STRICT_SCORE_COLUMNS = {
    "target_id", "canonical_smiles", "label_role",
    "l1", "l2", "l3", "l4", "final_score", "gate_status",
    "l1_status", "layer2_status", "l3_status", "l4_status",
    "l1_backend", "layer2_backend", "l3_backend", "l4_backend",
    "l1_model_asset_id", "layer2_model_asset_id", "l3_model_asset_id", "l4_model_asset_id",
}


def verify_score_frame(
    frame: pd.DataFrame,
    *,
    target_id: str,
    pool: pd.DataFrame | None = None,
    strict: bool = True,
) -> dict[str, object]:
    failures = []
    required = {
        "target_id", "canonical_smiles", "label_role", "l1", "l2", "l3", "l4",
        "final_score", "l1_status", "layer2_status", "l3_status", "l4_status",
    }
    if strict:
        required |= STRICT_SCORE_COLUMNS
    missing = required - set(frame.columns)
    if missing:
        return {"ok": False, "failures": [f"missing_columns={sorted(missing)}"]}
    if not frame["target_id"].astype(str).eq(str(target_id)).all():
        failures.append("target_id_mismatch")
    if frame["canonical_smiles"].duplicated().any():
        failures.append("duplicate_smiles")
    if not frame["label_role"].isin(
        ["heldout_positive", "heldout_negative", "unlabeled_background"]
    ).all():
        failures.append("invalid_label_role")
    if pool is not None:
        if frame["canonical_smiles"].astype(str).tolist() != pool["canonical_smiles"].astype(str).tolist():
            failures.append("pool_smiles_mismatch")
        if frame["label_role"].astype(str).tolist() != pool["label_role"].astype(str).tolist():
            failures.append("pool_label_role_mismatch")
    if strict:
        for column in ("l1_status", "layer2_status", "l3_status", "l4_status"):
            if not frame[column].eq("ok").all():
                failures.append(f"non_ok_{column}")
        if not frame["gate_status"].isin(["PASS", "FAIL"]).all():
            failures.append("invalid_gate_status")
        for column in (
            "l1_backend", "layer2_backend", "l3_backend", "l4_backend",
            "l1_model_asset_id", "layer2_model_asset_id", "l3_model_asset_id", "l4_model_asset_id",
        ):
            if frame[column].isna().any() or frame[column].astype(str).str.len().eq(0).any():
                failures.append(f"missing_{column}")
    numeric = {}
    for column in ("l1", "l2", "l3", "l4", "final_score"):
        values = pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)
        numeric[column] = values
        if not np.isfinite(values).all():
            failures.append(f"nonfinite_{column}")
        if strict and np.isfinite(values).all() and ((values < 0.0) | (values > 1.0)).any():
            failures.append(f"out_of_range_{column}")
    expected = np.round(
        0.20 * numeric["l1"] + 0.50 * numeric["l2"]
        + 0.20 * numeric["l3"] + 0.10 * numeric["l4"],
        4,
    )
    if strict and not np.allclose(numeric["final_score"], expected, atol=1e-4, rtol=0, equal_nan=False):
        failures.append("final_formula_mismatch")
    return {"ok": not failures, "failures": sorted(set(failures))}


def _same_json_value(left: object, right: object) -> bool:
    return canonicalize_for_compare(left) == canonicalize_for_compare(right)


def canonicalize_for_compare(value: object) -> object:
    if isinstance(value, dict):
        return {key: canonicalize_for_compare(item) for key, item in sorted(value.items())}
    if isinstance(value, list):
        return [canonicalize_for_compare(item) for item in value]
    if isinstance(value, float) and pd.isna(value):
        return None
    return value


def verify_docking_target(
    run_dir: str | Path,
    *,
    target_id: str,
    expected_top_n: int = 300,
) -> dict[str, object]:
    root = Path(run_dir)
    dock_dir = root / "docking" / target_id
    failures = []
    required = {
        "scores": root / "scores" / f"target_id={target_id}" / "part-000000.parquet",
        "selected": dock_dir / "selected_results.parquet",
        "fused": dock_dir / "fused_scores.parquet",
        "summary": dock_dir / "docking_summary.json",
        "metrics": dock_dir / "metrics_before_after.json",
    }
    missing = [name for name, path in required.items() if not path.is_file()]
    if missing:
        return {"ok": False, "failures": [f"missing={missing}"]}
    scores = pd.read_parquet(required["scores"])
    selected = pd.read_parquet(required["selected"])
    fused = pd.read_parquet(required["fused"])
    summary = json.loads(required["summary"].read_text(encoding="utf-8"))
    stored_metrics = json.loads(required["metrics"].read_text(encoding="utf-8"))
    if len(selected) != expected_top_n:
        failures.append("docking_row_count_mismatch")
    if selected["candidate_index"].duplicated().any():
        failures.append("duplicate_docking_candidate")
    expected = scores.sort_values(
        ["l2", "candidate_index"], ascending=[False, True], kind="mergesort"
    ).head(expected_top_n)
    expected_indices = expected["candidate_index"].astype(int).tolist()
    actual_ranked = selected.sort_values("selection_rank")
    if actual_ranked["candidate_index"].astype(int).tolist() != expected_indices:
        failures.append("docking_selection_not_l2_topn")
    score_smiles = scores.set_index("candidate_index")["canonical_smiles"].astype(str).to_dict()
    if any(score_smiles.get(int(row.candidate_index)) != str(row.canonical_smiles) for row in selected.itertuples()):
        failures.append("docking_smiles_mismatch")
    if selected["status"].isna().any() or selected["status"].astype(str).str.len().eq(0).any():
        failures.append("missing_docking_status")
    status_counts = selected["status"].value_counts().to_dict()
    if int(summary.get("n_submitted", -1)) != expected_top_n:
        failures.append("docking_summary_count_mismatch")
    if summary.get("status_counts") != status_counts:
        failures.append("docking_status_counts_mismatch")
    if len(fused) != len(scores) or fused["candidate_index"].astype(int).tolist() != scores["candidate_index"].astype(int).tolist():
        failures.append("fused_score_alignment_mismatch")
    else:
        selected_status = selected.set_index("candidate_index")["status"].astype(str).to_dict()
        expected_status = [selected_status.get(int(index), "not_selected") for index in scores["candidate_index"]]
        if fused["docking_status"].astype(str).tolist() != expected_status:
            failures.append("fused_docking_status_mismatch")
        unchanged = np.asarray([status != "ok" for status in expected_status])
        if not np.allclose(
            pd.to_numeric(fused.loc[unchanged, "fused_score"], errors="coerce"),
            pd.to_numeric(scores.loc[unchanged, "l2"], errors="coerce"),
            atol=1e-12,
            rtol=0,
            equal_nan=False,
        ):
            failures.append("non_successful_docking_changed_baseline")
        try:
            weight = float(stored_metrics["weight"])
            base = pd.to_numeric(scores["l2"], errors="coerce").to_numpy(dtype=float)
            le_by_index = selected.set_index("candidate_index")["ligand_efficiency"].to_dict()
            le = np.asarray([
                float(le_by_index.get(int(index), np.nan))
                if selected_status.get(int(index)) == "ok" else np.nan
                for index in scores["candidate_index"]
            ])
            expected_fused = base.copy()
            mask = np.isfinite(le)
            if mask.any():
                standardized = -le[mask]
                standardized = (standardized - standardized.mean()) / (standardized.std() + 1e-9)
                expected_fused[mask] = base[mask] + weight * (base.std() + 1e-9) * standardized
            actual_fused = pd.to_numeric(fused["fused_score"], errors="coerce").to_numpy(dtype=float)
            if not np.allclose(actual_fused, expected_fused, atol=1e-12, rtol=0, equal_nan=False):
                failures.append("successful_docking_fusion_mismatch")
        except (KeyError, TypeError, ValueError):
            failures.append("docking_fusion_unverifiable")
    recomputed_before = metrics.compute_target_metrics(scores)
    recomputed_after = metrics.compute_target_metrics(fused.assign(final_score=fused["fused_score"]))
    if not _same_json_value(stored_metrics.get("before"), recomputed_before):
        failures.append("docking_before_metrics_mismatch")
    if not _same_json_value(stored_metrics.get("after"), recomputed_after):
        failures.append("docking_after_metrics_mismatch")
    return {
        "ok": not failures,
        "failures": sorted(set(failures)),
        "status_counts": status_counts,
        "n_submitted": int(len(selected)),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact_files(run_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in run_dir.rglob("*")
        if path.is_file() and path.name != MANIFEST_NAME and not path.name.endswith(".tmp")
    )


def write_manifest(run_dir: str | Path) -> Path:
    root = Path(run_dir).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"run directory not found: {root}")
    lines = [f"{_sha256(path)}  {path.relative_to(root).as_posix()}" for path in _artifact_files(root)]
    destination = root / MANIFEST_NAME
    temporary = root / f".{MANIFEST_NAME}.tmp"
    temporary.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    os.replace(temporary, destination)
    return destination


def verify_manifest(run_dir: str | Path) -> dict[str, object]:
    root = Path(run_dir).resolve()
    manifest = root / MANIFEST_NAME
    if not manifest.is_file():
        return {"ok": False, "mismatches": [MANIFEST_NAME], "missing": [MANIFEST_NAME], "unexpected": []}

    expected: dict[str, str] = {}
    malformed = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        parts = line.split("  ", 1)
        if len(parts) != 2 or len(parts[0]) != 64:
            malformed.append(line)
            continue
        expected[parts[1]] = parts[0]

    actual_paths = {path.relative_to(root).as_posix(): path for path in _artifact_files(root)}
    missing = sorted(set(expected) - set(actual_paths))
    unexpected = sorted(set(actual_paths) - set(expected))
    mismatches = sorted(
        relative
        for relative in set(expected) & set(actual_paths)
        if _sha256(actual_paths[relative]) != expected[relative]
    )
    mismatches.extend(missing)
    if malformed:
        mismatches.append(MANIFEST_NAME)
    mismatches = sorted(set(mismatches))
    return {
        "ok": not mismatches and not unexpected,
        "mismatches": mismatches,
        "missing": missing,
        "unexpected": unexpected,
        "malformed_lines": malformed,
        "n_files": len(actual_paths),
    }


def verify_artifacts(
    run_dir: str | Path,
    *,
    expected_targets: int | None = None,
    pool_size: int | None = None,
    strict: bool = True,
    require_docking: bool = False,
) -> dict[str, object]:
    """Recompute cardinality, status and contamination checks from candidate shards."""
    root = Path(run_dir).resolve()
    failures: list[str] = []
    checks: dict[str, object] = {}
    target_path = root / "target_manifest.parquet"
    if not target_path.is_file():
        return {"ok": False, "checks": {}, "failures": ["target_manifest_missing"]}
    targets = pd.read_parquet(target_path)
    target_ids = targets["target_id"].astype(str).tolist()
    checks["n_targets"] = len(target_ids)
    if expected_targets is not None and len(target_ids) != expected_targets:
        failures.append(f"n_targets={len(target_ids)} expected={expected_targets}")
    if strict and len(target_ids) != 1000:
        failures.append(f"strict_n_targets={len(target_ids)} expected=1000")
    if len(set(target_ids)) != len(target_ids):
        failures.append("duplicate_target_ids")

    required_columns = {
        "target_id", "canonical_smiles", "label_role", "l1", "l2", "l3", "l4",
        "final_score", "l1_status", "layer2_status", "l3_status", "l4_status",
    }
    score_paths = sorted((root / "scores").glob("target_id=*/part-*.parquet"))
    total_rows = 0
    per_target_rows: dict[str, int] = {}
    labeled_pairs: set[tuple[str, str]] = set()
    all_pairs: set[tuple[str, str]] = set()
    l4_values: list[float] = []
    target_lookup = targets.set_index("target_id")
    for path in score_paths:
        frame = pd.read_parquet(path)
        target_id = path.parent.name.split("=", 1)[-1]
        per_target_rows[target_id] = per_target_rows.get(target_id, 0) + len(frame)
        total_rows += len(frame)
        pool = None
        if target_id in target_lookup.index and "pool_file" in target_lookup.columns:
            pool_file = target_lookup.loc[target_id, "pool_file"]
            if pd.notna(pool_file) and str(pool_file):
                pool_path = root / str(pool_file)
                if pool_path.is_file():
                    pool = pd.read_parquet(pool_path).reset_index(drop=True)
                    expected_hash = str(target_lookup.loc[target_id].get("pool_sha256", ""))
                    actual_hash = _sha256(pool_path)
                    if expected_hash and expected_hash != actual_hash:
                        failures.append(f"{target_id}:pool_sha256_mismatch")
        elif strict:
            failures.append(f"{target_id}:pool_file_missing")
        frame_check = verify_score_frame(frame, target_id=target_id, pool=pool, strict=strict)
        failures.extend(f"{target_id}:{failure}" for failure in frame_check.get("failures", []))
        if strict and "l4" in frame:
            values = pd.to_numeric(frame["l4"], errors="coerce").to_numpy(dtype=float)
            l4_values.extend(values[np.isfinite(values)].tolist())
        labeled = frame[frame["label_role"].isin(["heldout_positive", "heldout_negative"])]
        labeled_pairs.update(zip(labeled["target_id"].astype(str), labeled["canonical_smiles"].astype(str)))
        all_pairs.update(zip(frame["target_id"].astype(str), frame["canonical_smiles"].astype(str)))

    checks["n_score_partitions"] = len(score_paths)
    checks["n_score_rows"] = total_rows
    if set(per_target_rows) != set(target_ids):
        failures.append("score_target_set_mismatch")
    expected_pool = pool_size
    if expected_pool is None and "pool_size" in targets:
        values = targets["pool_size"].dropna().astype(int).unique().tolist()
        expected_pool = values[0] if len(values) == 1 else None
    if expected_pool is None:
        failures.append("pool_size_unknown")
    else:
        for target_id in target_ids:
            if per_target_rows.get(target_id) != expected_pool:
                failures.append(
                    f"{target_id}:rows={per_target_rows.get(target_id, 0)} expected={expected_pool}"
                )
    checks["pool_size"] = expected_pool
    checks["expected_total_rows"] = len(target_ids) * expected_pool if expected_pool else None
    if checks["expected_total_rows"] != total_rows:
        failures.append(f"n_score_rows={total_rows} expected={checks['expected_total_rows']}")

    fit_path = root / "reconstructed_fit_rows.parquet"
    contamination = 0
    if fit_path.is_file():
        fit = pd.read_parquet(fit_path).dropna(subset=["target_id", "canonical_smiles"])
        fit_pairs = set(zip(fit["target_id"].astype(str), fit["canonical_smiles"].astype(str)))
        contamination = len(fit_pairs & labeled_pairs)
        all_contamination = len(fit_pairs & all_pairs)
        if contamination:
            failures.append(f"fit_pair_contamination={contamination}")
        if all_contamination:
            failures.append(f"fit_pair_contamination_all={all_contamination}")
    else:
        failures.append("reconstructed_fit_rows_missing")
        all_contamination = None
    checks["fit_pair_contamination"] = contamination
    checks["fit_pair_contamination_all"] = all_contamination
    if strict and (len(l4_values) < 2 or float(np.std(l4_values)) <= 0.0):
        failures.append("l4_constant_or_missing")

    parity_path = root / "parity_evidence.json"
    if parity_path.is_file():
        parity_result = verify_parity_evidence(root, strict=strict)
        checks["parity_passed"] = bool(parity_result["ok"])
        if not parity_result["ok"]:
            failures.append("parity_gate_failed")
    else:
        checks["parity_passed"] = False
        failures.append("parity_evidence_missing")

    if require_docking:
        for target in targets.itertuples(index=False):
            if bool(getattr(target, "receptor_available", False)):
                docking_result = verify_docking_target(root, target_id=str(target.target_id), expected_top_n=300)
                failures.extend(f"{target.target_id}:{failure}" for failure in docking_result.get("failures", []))
    checks["docking_required"] = require_docking
    return {"ok": not failures, "checks": checks, "failures": failures}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify a frozen four-level CLI run")
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)
    artifact_result = verify_artifacts(args.run_dir, strict=args.strict, require_docking=args.strict)
    manifest_result = verify_manifest(args.run_dir)
    result = {
        "ok": bool(artifact_result["ok"] and manifest_result["ok"]),
        "artifacts": artifact_result,
        "manifest": manifest_result,
    }
    provenance_result = verify_provenance(args.run_dir)
    result["provenance"] = provenance_result
    if not provenance_result["ok"]:
        result["ok"] = False
    if args.strict and artifact_result["ok"]:
        try:
            from . import report
            recomputed = report.recompute_summary(args.run_dir)
            summary_result = verify_summary_bytes(args.run_dir, recomputed)
        except Exception as exc:
            summary_result = {"ok": False, "reason": f"recompute_error:{type(exc).__name__}:{exc}"}
        result["summary"] = summary_result
        if not summary_result["ok"]:
            result["ok"] = False
    print(result)
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
