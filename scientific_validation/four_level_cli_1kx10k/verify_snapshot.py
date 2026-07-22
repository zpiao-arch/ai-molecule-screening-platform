"""Verify the compact evidence snapshot shipped in the source archive."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd

from . import verify_run


REQUIRED_CHECKPOINTS = ("prepare", "layer_cache", "score", "dock", "report")


def _read_json(root: Path, relative: str, failures: list[str]) -> dict[str, Any] | None:
    path = root / relative
    if not path.is_file():
        failures.append(f"missing:{relative}")
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        failures.append(f"invalid_json:{relative}")
        return None
    if not isinstance(payload, dict):
        failures.append(f"invalid_object:{relative}")
        return None
    return payload


def verify_snapshot(snapshot_dir: str | Path) -> dict[str, object]:
    root = Path(snapshot_dir).expanduser().resolve()
    failures: list[str] = []
    snapshot = _read_json(root, "snapshot.json", failures)
    if snapshot is None:
        return {
            "ok": False,
            "schema": None,
            "claims_row_level_replay": False,
            "failures": failures,
        }

    if snapshot.get("schema") != "four-level-cli.compact-snapshot.v1":
        failures.append("snapshot_schema_mismatch")
    full_manifest = snapshot.get("complete_run_manifest_sha256")
    if not isinstance(full_manifest, str) or len(full_manifest) < 8:
        failures.append("complete_run_manifest_sha256_missing")
    expected = snapshot.get("expected")
    if not isinstance(expected, dict):
        failures.append("snapshot_expected_missing")
        expected = {}

    manifest = verify_run.verify_manifest(root)
    if not manifest["ok"]:
        failures.append("compact_manifest_invalid")

    summary = _read_json(root, "summary.json", failures)
    score_summary = _read_json(root, "score_summary.json", failures)
    prepare = _read_json(root, "prepare_summary.json", failures)
    layer_cache = _read_json(root, "layer_cache_summary.json", failures)
    failure_manifest = _read_json(root, "failure_manifest.json", failures)
    parity = _read_json(root, "parity_evidence.json", failures)

    n_targets = int(expected.get("n_targets", 0) or 0)
    n_pairs = int(expected.get("n_pairs", 0) or 0)
    partitions = int(expected.get("score_partitions", 0) or 0)
    if n_targets <= 0 or n_pairs <= 0 or partitions <= 0:
        failures.append("snapshot_expected_counts_invalid")

    if summary is not None:
        if summary.get("n_targets") != n_targets:
            failures.append("summary_n_targets_mismatch")
        if summary.get("n_pairs") != n_pairs:
            failures.append("summary_n_pairs_mismatch")
        if summary.get("fit_pair_contamination") != 0:
            failures.append("summary_fit_pair_contamination_nonzero")
        layer_failures = summary.get("layer_failures") or {}
        if any(layer_failures.get(layer, 1) != 0 for layer in ("l1", "l2", "l3", "l4")):
            failures.append("summary_layer_failures_nonzero")
    if score_summary is not None:
        if score_summary.get("n_targets") != n_targets:
            failures.append("score_summary_n_targets_mismatch")
        if score_summary.get("n_pairs") != n_pairs:
            failures.append("score_summary_n_pairs_mismatch")
    if prepare is not None:
        if prepare.get("n_targets") != n_targets:
            failures.append("prepare_summary_n_targets_mismatch")
        if prepare.get("n_pairs") != n_pairs:
            failures.append("prepare_summary_n_pairs_mismatch")
        if int(prepare.get("pool_size", 0) or 0) <= 0:
            failures.append("prepare_summary_pool_size_invalid")
    if layer_cache is not None and any(
        int(layer_cache.get(key, 1) or 0) != 0 for key in ("l1_failures", "l3_failures", "l4_failures")
    ):
        failures.append("layer_cache_failures_nonzero")
    if failure_manifest is not None:
        layer_failures = failure_manifest.get("layer_failures") or {}
        if any(int(layer_failures.get(layer, 1) or 0) != 0 for layer in ("l1", "l3", "l4")):
            failures.append("failure_manifest_layer_failures_nonzero")
        interpretation = failure_manifest.get("interpretation") or {}
        for key in ("docking_failures_are_explicit", "pair_heldout_is_not_cold_start", "unlabeled_background_is_not_confirmed_negative"):
            if interpretation.get(key) is not True:
                failures.append(f"failure_manifest_interpretation_missing:{key}")
    if parity is not None:
        if parity.get("passed") is not True:
            failures.append("parity_failed")
        if parity.get("route_branch_match") is not True:
            failures.append("parity_route_branch_mismatch")
        if (parity.get("invalid_input") or {}).get("passed") is not True:
            failures.append("parity_invalid_input_missing")

    metrics_path = root / "per_target_metrics.parquet"
    metrics = None
    if not metrics_path.is_file():
        failures.append("missing:per_target_metrics.parquet")
    else:
        try:
            metrics = pd.read_parquet(metrics_path)
        except Exception:
            failures.append("invalid:per_target_metrics.parquet")
    if metrics is not None:
        if len(metrics) != partitions:
            failures.append("metrics_partition_count_mismatch")
        if "target_id" not in metrics or metrics["target_id"].astype(str).nunique() != len(metrics):
            failures.append("metrics_target_ids_invalid")
        if "n_pool" not in metrics or metrics["n_pool"].astype(int).nunique() != 1:
            failures.append("metrics_pool_sizes_invalid")

    for stage in REQUIRED_CHECKPOINTS:
        checkpoint = _read_json(root, f"checkpoints/{stage}.json", failures)
        if checkpoint is not None and (
            checkpoint.get("stage") != stage or checkpoint.get("status") != "complete"
        ):
            failures.append(f"checkpoint_not_complete:{stage}")

    docking = _read_json(root, "docking/CHEMBL2051/docking_summary.json", failures)
    if docking is not None:
        if int(docking.get("n_submitted", 0) or 0) <= 0:
            failures.append("docking_summary_empty")
        if not isinstance(docking.get("status_counts"), dict):
            failures.append("docking_summary_status_counts_missing")

    return {
        "ok": not failures,
        "schema": snapshot.get("schema"),
        "claims_row_level_replay": False,
        "expected": expected,
        "manifest": manifest,
        "failures": sorted(set(failures)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Verify compact four-level CLI evidence")
    parser.add_argument("--snapshot-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    result = verify_snapshot(args.snapshot_dir)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
