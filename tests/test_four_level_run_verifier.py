from scientific_validation.four_level_cli_1kx10k import verify_run
from scientific_validation.four_level_cli_1kx10k import metrics
import json
import numpy as np
import pandas as pd
import pytest


def test_manifest_detects_score_shard_tampering(tmp_path):
    score_dir = tmp_path / "scores" / "target_id=CHEMBL1"
    score_dir.mkdir(parents=True)
    shard = score_dir / "part-000000.parquet"
    shard.write_bytes(b"original-score-shard")
    (tmp_path / "summary.json").write_text("{}", encoding="utf-8")

    verify_run.write_manifest(tmp_path)
    assert verify_run.verify_manifest(tmp_path)["ok"] is True

    shard.write_bytes(b"tampered-score-shard")
    result = verify_run.verify_manifest(tmp_path)

    assert result["ok"] is False
    assert "scores/target_id=CHEMBL1/part-000000.parquet" in result["mismatches"]


def test_manifest_ignores_stale_atomic_temporary_files(tmp_path):
    (tmp_path / "summary.json").write_text("{}", encoding="utf-8")
    (tmp_path / ".MANIFEST.sha256.tmp").write_text("stale", encoding="utf-8")

    verify_run.write_manifest(tmp_path)

    result = verify_run.verify_manifest(tmp_path)
    assert result["ok"] is True
    assert ".MANIFEST.sha256.tmp" not in (tmp_path / "MANIFEST.sha256").read_text(encoding="utf-8")


def test_run_verifier_checks_cardinality_status_and_fit_contamination(tmp_path):
    target_manifest = pd.DataFrame(
        {
            "target_id": ["CHEMBL1", "CHEMBL2"],
            "pool_size": [3, 3],
            "route_branch": ["library", "library"],
            "receptor_available": [False, False],
        }
    )
    target_manifest.to_parquet(tmp_path / "target_manifest.parquet", index=False)
    pd.DataFrame(
        {"target_id": ["CHEMBL1"], "canonical_smiles": ["FIT"]}
    ).to_parquet(tmp_path / "reconstructed_fit_rows.parquet", index=False)
    (tmp_path / "parity_evidence.json").write_text('{"passed": true}', encoding="utf-8")
    required = {
        "target_id": ["CHEMBL1"] * 3,
        "canonical_smiles": ["P", "N", "U"],
        "label_role": ["heldout_positive", "heldout_negative", "unlabeled_background"],
        "l1": [0.5, 0.4, 0.3], "l2": [0.6, 0.2, 0.1], "l3": [0.8, 0.8, 0.8],
        "l4": [0.7, 0.6, 0.5], "final_score": [0.6, 0.4, 0.3],
        "l1_status": ["ok"] * 3, "layer2_status": ["ok"] * 3,
        "l3_status": ["ok"] * 3, "l4_status": ["ok"] * 3,
    }
    score_dir = tmp_path / "scores" / "target_id=CHEMBL1"
    score_dir.mkdir(parents=True)
    pd.DataFrame(required).to_parquet(score_dir / "part-000000.parquet", index=False)
    score_dir = tmp_path / "scores" / "target_id=CHEMBL2"
    score_dir.mkdir(parents=True)
    required["target_id"] = ["CHEMBL2"] * 3
    pd.DataFrame(required).to_parquet(score_dir / "part-000000.parquet", index=False)

    result = verify_run.verify_artifacts(tmp_path, expected_targets=2, pool_size=3, strict=False)

    assert result["ok"] is True
    assert result["checks"]["n_score_rows"] == 6


def test_stage_checkpoint_rejects_changed_inputs(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import provenance

    provenance.write_stage_checkpoint(
        tmp_path,
        "score",
        inputs={"target_manifest_sha256": "aaa", "model_sha256": "bbb"},
        status="in_progress",
    )
    assert provenance.validate_stage_checkpoint(
        tmp_path,
        "score",
        inputs={"target_manifest_sha256": "aaa", "model_sha256": "bbb"},
    )["ok"] is True

    with pytest.raises(RuntimeError, match="input fingerprint mismatch"):
        provenance.require_stage_checkpoint(
            tmp_path,
            "score",
            inputs={"target_manifest_sha256": "changed", "model_sha256": "bbb"},
        )


def test_summary_verifier_requires_byte_identical_recomputation(tmp_path):
    payload = {"n_targets": 2, "aggregates": {"auc": {"median": 0.75}}}
    (tmp_path / "summary.json").write_text(
        verify_run.render_canonical_json(payload), encoding="utf-8"
    )

    assert verify_run.verify_summary_bytes(tmp_path, payload)["ok"] is True
    changed = {"n_targets": 2, "aggregates": {"auc": {"median": 0.76}}}
    assert verify_run.verify_summary_bytes(tmp_path, changed)["ok"] is False


def test_strict_provenance_verifier_requires_archived_design_logs_and_checkpoints(tmp_path):
    result = verify_run.verify_provenance(tmp_path)

    assert result["ok"] is False
    assert set(result["missing"]) >= {
        "DESIGN.md", "run.log", "failure_manifest.json",
        "checkpoints/prepare.json", "checkpoints/layer_cache.json",
        "checkpoints/score.json", "checkpoints/dock.json", "checkpoints/report.json",
    }


def test_report_prepare_does_not_fabricate_run_log(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import report

    report.prepare_delivery_artifacts(tmp_path)

    assert not (tmp_path / "run.log").exists()


def test_stage_checkpoint_appends_structured_run_event(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import provenance

    checkpoint = provenance.write_stage_checkpoint(
        tmp_path,
        "score",
        inputs={"target_manifest_sha256": "abc"},
        status="complete",
        outputs={"n_rows": 10},
    )

    events = [json.loads(line) for line in (tmp_path / "run.log").read_text(encoding="utf-8").splitlines()]
    assert events == [{
        "checkpoint": "checkpoints/score.json",
        "event": "stage_checkpoint",
        "input_fingerprint": json.loads(checkpoint.read_text(encoding="utf-8"))["input_fingerprint"],
        "stage": "score",
        "status": "complete",
        "timestamp_utc": json.loads(checkpoint.read_text(encoding="utf-8"))["updated_at_utc"],
    }]


def test_failed_stage_checkpoint_records_error_and_event(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import provenance

    provenance.write_stage_checkpoint(
        tmp_path,
        "prepare",
        inputs={"targets": 0},
        status="in_progress",
    )
    checkpoint = provenance.fail_stage_checkpoint(
        tmp_path,
        "prepare",
        ValueError("targets must be greater than zero"),
    )

    payload = json.loads(checkpoint.read_text(encoding="utf-8"))
    assert payload["status"] == "failed"
    assert payload["metadata"]["error_type"] == "ValueError"
    assert "greater than zero" in payload["metadata"]["error_message"]
    assert provenance.validate_stage_checkpoint(
        tmp_path,
        "prepare",
        inputs={"targets": 0},
        require_complete=True,
    )["reason"] == "checkpoint is not complete"
    events = [json.loads(line) for line in (tmp_path / "run.log").read_text(encoding="utf-8").splitlines()]
    assert events[-1]["status"] == "failed"
    assert events[-1]["stage"] == "prepare"


def test_batch_main_marks_stage_failed_when_stage_raises(tmp_path, monkeypatch):
    from scientific_validation.four_level_cli_1kx10k import batch_cli, provenance

    def fail_prepare(config, run_dir, *, n_targets, pool_size, seed):
        del config, n_targets, pool_size, seed
        provenance.write_stage_checkpoint(
            run_dir,
            "prepare",
            inputs={"test": True},
            status="in_progress",
        )
        raise RuntimeError("synthetic prepare failure")

    monkeypatch.setattr(batch_cli, "prepare_run", fail_prepare)

    with pytest.raises(RuntimeError, match="synthetic prepare failure"):
        batch_cli.main(["--run-dir", str(tmp_path), "prepare"])

    payload = json.loads((tmp_path / "checkpoints" / "prepare.json").read_text(encoding="utf-8"))
    assert payload["status"] == "failed"
    assert payload["metadata"]["error_type"] == "RuntimeError"


def _write_compact_snapshot_fixture(root, *, tamper_summary=False):
    from scientific_validation.four_level_cli_1kx10k import verify_run

    root.mkdir(parents=True, exist_ok=True)
    expected_targets = 2
    expected_pairs = 6
    (root / "snapshot.json").write_text(
        json.dumps(
            {
                "schema": "four-level-cli.compact-snapshot.v1",
                "complete_run_manifest_sha256": "full-run-manifest",
                "expected": {
                    "n_targets": expected_targets,
                    "n_pairs": expected_pairs,
                    "score_partitions": expected_targets,
                },
                "omitted_artifact_classes": ["scores", "pool_manifest", "layer_cache_shards"],
            }
        ),
        encoding="utf-8",
    )
    summary_targets = expected_targets + 1 if tamper_summary else expected_targets
    (root / "summary.json").write_text(
        json.dumps(
            {
                "n_targets": summary_targets,
                "n_pairs": expected_pairs,
                "fit_pair_contamination": 0,
                "layer_failures": {"l1": 0, "l2": 0, "l3": 0, "l4": 0},
            }
        ),
        encoding="utf-8",
    )
    (root / "score_summary.json").write_text(
        json.dumps({"n_targets": expected_targets, "n_pairs": expected_pairs}),
        encoding="utf-8",
    )
    (root / "prepare_summary.json").write_text(
        json.dumps({"n_targets": expected_targets, "n_pairs": expected_pairs, "pool_size": 3}),
        encoding="utf-8",
    )
    (root / "layer_cache_summary.json").write_text(
        json.dumps({"n_molecules": 3, "l1_failures": 0, "l3_failures": 0, "l4_failures": 0}),
        encoding="utf-8",
    )
    (root / "failure_manifest.json").write_text(
        json.dumps(
            {
                "layer_failures": {"l1": 0, "l3": 0, "l4": 0},
                "interpretation": {
                    "docking_failures_are_explicit": True,
                    "pair_heldout_is_not_cold_start": True,
                    "unlabeled_background_is_not_confirmed_negative": True,
                },
            }
        ),
        encoding="utf-8",
    )
    (root / "parity_evidence.json").write_text(
        json.dumps({"passed": True, "n_molecules": 20, "route_branch_match": True, "invalid_input": {"passed": True}}),
        encoding="utf-8",
    )
    metrics_frame = pd.DataFrame(
        {
            "target_id": ["T1", "T2"],
            "n_pool": [3, 3],
            "evaluation_tier": ["strict_labeled", "strict_labeled"],
            "route_branch": ["library", "library"],
        }
    )
    metrics_frame.to_parquet(root / "per_target_metrics.parquet", index=False)
    docking = root / "docking" / "CHEMBL2051"
    docking.mkdir(parents=True)
    (docking / "docking_summary.json").write_text(
        json.dumps({"target_id": "CHEMBL2051", "n_submitted": 300, "status_counts": {"ok": 1}}),
        encoding="utf-8",
    )
    checkpoints = root / "checkpoints"
    checkpoints.mkdir()
    for stage in ("prepare", "layer_cache", "score", "dock", "report"):
        (checkpoints / f"{stage}.json").write_text(
            json.dumps({"stage": stage, "status": "complete"}),
            encoding="utf-8",
        )
    verify_run.write_manifest(root)


def test_compact_snapshot_verifier_accepts_consistent_evidence(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import verify_snapshot

    _write_compact_snapshot_fixture(tmp_path)

    result = verify_snapshot.verify_snapshot(tmp_path)

    assert result["ok"] is True
    assert result["claims_row_level_replay"] is False


def test_compact_snapshot_verifier_rejects_cross_document_count_drift(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import verify_snapshot

    _write_compact_snapshot_fixture(tmp_path, tamper_summary=True)

    result = verify_snapshot.verify_snapshot(tmp_path)

    assert result["ok"] is False
    assert "summary_n_targets_mismatch" in result["failures"]


def test_strict_provenance_rejects_placeholder_files(tmp_path):
    required = [
        "DESIGN.md", "run.log", "failure_manifest.json",
        "checkpoints/prepare.json", "checkpoints/layer_cache.json",
        "checkpoints/score.json", "checkpoints/dock.json", "checkpoints/report.json",
    ]
    for relative in required:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("placeholder\n", encoding="utf-8")

    result = verify_run.verify_provenance(tmp_path)

    assert result["ok"] is False
    assert "run.log" in result["invalid"]


def test_strict_provenance_accepts_legacy_run_index_with_real_stage_log(tmp_path):
    (tmp_path / "DESIGN.md").write_text("design\n", encoding="utf-8")
    (tmp_path / "failure_manifest.json").write_text("{}\n", encoding="utf-8")
    (tmp_path / "run.log").write_text(
        "Four-level CLI frozen run log\nDetailed stdout/stderr: pipeline.log\n",
        encoding="utf-8",
    )
    (tmp_path / "pipeline.log").write_text("stage output\n", encoding="utf-8")
    for stage in ("prepare", "layer_cache", "score", "dock", "report"):
        path = tmp_path / "checkpoints" / f"{stage}.json"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps({"stage": stage, "status": "complete"}),
            encoding="utf-8",
        )

    result = verify_run.verify_provenance(tmp_path)

    assert result == {"ok": True, "missing": [], "invalid": []}


def test_report_prepare_resolves_design_relative_to_packaged_repo(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import report

    repo = tmp_path / "repo"
    (repo / "docs").mkdir(parents=True)
    (repo / "docs" / "design.md").write_text("design", encoding="utf-8")
    run_dir = repo / "runs" / "run-1"
    run_dir.mkdir(parents=True)
    (run_dir / "run_config.json").write_text(
        json.dumps({"design": "docs/design.md", "scoring_dir": str(repo / "scoring")}),
        encoding="utf-8",
    )

    report.prepare_delivery_artifacts(run_dir)

    assert (run_dir / "DESIGN.md").read_text(encoding="utf-8") == "design"


def test_report_generation_seals_manifest(tmp_path, monkeypatch):
    from scientific_validation.four_level_cli_1kx10k import report

    for relative in (
        "target_manifest.parquet",
        "score_summary.json",
        "docking/CHEMBL2051/metrics_before_after.json",
    ):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("placeholder", encoding="utf-8")
    monkeypatch.setattr(report, "prepare_delivery_artifacts", lambda root: None)
    monkeypatch.setattr(report, "recompute_summary", lambda root: {"n_targets": 1})

    report.generate_report(tmp_path)

    assert (tmp_path / "MANIFEST.sha256").is_file()


def test_strict_parity_evidence_requires_full_fields_route_and_invalid_input(tmp_path):
    incomplete = {
        "n_molecules": 20,
        "l2_max_abs_diff": 0.0,
        "tolerance": 1e-4,
        "passed": True,
    }
    (tmp_path / "parity_evidence.json").write_text(json.dumps(incomplete), encoding="utf-8")

    result = verify_run.verify_parity_evidence(tmp_path, strict=True)

    assert result["ok"] is False
    assert set(result["missing_numeric_fields"]) == {"l1", "l2", "l3", "l4", "final_score"}
    assert "route_branch_match" in result["missing_checks"]
    assert "invalid_input" in result["missing_checks"]


def test_score_frame_verifier_detects_formula_drift_nonfinite_values_and_pool_mismatch():
    pool = pd.DataFrame(
        {
            "canonical_smiles": ["CCO", "CCN"],
            "label_role": ["heldout_positive", "unlabeled_background"],
        }
    )
    frame = pd.DataFrame(
        {
            "target_id": ["CHEMBL1", "CHEMBL1"],
            "canonical_smiles": ["CCO", "CCN"],
            "label_role": ["heldout_positive", "unlabeled_background"],
            "l1": [0.8, 0.7], "l2": [0.6, 0.5], "l3": [0.9, 0.8], "l4": [0.4, 0.3],
            "final_score": [0.68, 0.58],
            "gate_status": ["PASS", "PASS"],
            "l1_status": ["ok", "ok"], "layer2_status": ["ok", "ok"],
            "l3_status": ["ok", "ok"], "l4_status": ["ok", "ok"],
            "l1_backend": ["RDKit"] * 2, "layer2_backend": ["BindingDB"] * 2,
            "l3_backend": ["ADMET"] * 2, "l4_backend": ["UniMol"] * 2,
            "l1_model_asset_id": ["rdkit"] * 2, "layer2_model_asset_id": ["l2"] * 2,
            "l3_model_asset_id": ["admet"] * 2, "l4_model_asset_id": ["unimol"] * 2,
        }
    )
    assert verify_run.verify_score_frame(frame, target_id="CHEMBL1", pool=pool, strict=True)["ok"] is True

    frame.loc[0, "final_score"] = 0.0
    frame.loc[1, "l4"] = np.nan
    frame.loc[1, "label_role"] = "heldout_negative"
    result = verify_run.verify_score_frame(frame, target_id="CHEMBL1", pool=pool, strict=True)

    assert result["ok"] is False
    assert set(result["failures"]) >= {"nonfinite_l4", "final_formula_mismatch", "pool_label_role_mismatch"}


def test_score_frame_verifier_rejects_out_of_range_layer_values():
    frame = pd.DataFrame(
        {
            "target_id": ["CHEMBL1"],
            "canonical_smiles": ["CCO"],
            "label_role": ["heldout_positive"],
            "l1": [1.2], "l2": [0.6], "l3": [0.8], "l4": [0.4],
            "final_score": [0.74],
            "gate_status": ["PASS"],
            "l1_status": ["ok"], "layer2_status": ["ok"],
            "l3_status": ["ok"], "l4_status": ["ok"],
            "l1_backend": ["RDKit"], "layer2_backend": ["BindingDB"],
            "l3_backend": ["ADMET"], "l4_backend": ["UniMol"],
            "l1_model_asset_id": ["rdkit"], "layer2_model_asset_id": ["l2"],
            "l3_model_asset_id": ["admet"], "l4_model_asset_id": ["unimol"],
        }
    )

    result = verify_run.verify_score_frame(frame, target_id="CHEMBL1", strict=True)

    assert result["ok"] is False
    assert "out_of_range_l1" in result["failures"]


def test_docking_verifier_recomputes_topn_status_counts_and_fusion(tmp_path):
    target_id = "CHEMBL2051"
    score_dir = tmp_path / "scores" / f"target_id={target_id}"
    dock_dir = tmp_path / "docking" / target_id
    score_dir.mkdir(parents=True)
    dock_dir.mkdir(parents=True)
    scores = pd.DataFrame(
        {
            "candidate_index": [0, 1, 2, 3],
            "canonical_smiles": ["C", "CC", "CCC", "CCCC"],
            "label_role": ["heldout_positive", "heldout_negative", "unlabeled_background", "unlabeled_background"],
            "l2": [0.9, 0.8, 0.7, 0.6],
            "final_score": [0.85, 0.75, 0.65, 0.55],
        }
    )
    scores.to_parquet(score_dir / "part-000000.parquet", index=False)
    selected = pd.DataFrame(
        {
            "candidate_index": [0, 1],
            "canonical_smiles": ["C", "CC"],
            "selection_rank": [1, 2],
            "status": ["ok", "dock_timeout"],
            "affinity": [-7.0, np.nan],
            "heavy_atoms": [1, 2],
            "ligand_efficiency": [-7.0, np.nan],
        }
    )
    selected.to_parquet(dock_dir / "selected_results.parquet", index=False)
    fused = scores.copy()
    fused["docking_status"] = ["ok", "dock_timeout", "not_selected", "not_selected"]
    fused["fused_score"] = [0.9, 0.8, 0.7, 0.6]
    fused.to_parquet(dock_dir / "fused_scores.parquet", index=False)
    before = metrics.compute_target_metrics(scores)
    after = metrics.compute_target_metrics(fused.assign(final_score=fused["fused_score"]))
    (dock_dir / "metrics_before_after.json").write_text(
        json.dumps({"before": before, "after": after, "weight": 0.3}), encoding="utf-8"
    )
    (dock_dir / "docking_summary.json").write_text(
        json.dumps({"target_id": target_id, "n_submitted": 2, "top_n_requested": 2, "status_counts": {"ok": 1, "dock_timeout": 1}}),
        encoding="utf-8",
    )

    result = verify_run.verify_docking_target(tmp_path, target_id=target_id, expected_top_n=2)
    assert result["ok"] is True

    fused.loc[0, "fused_score"] = 0.95
    fused.to_parquet(dock_dir / "fused_scores.parquet", index=False)
    tampered_after = metrics.compute_target_metrics(fused.assign(final_score=fused["fused_score"]))
    (dock_dir / "metrics_before_after.json").write_text(
        json.dumps({"before": before, "after": tampered_after, "weight": 0.3}), encoding="utf-8"
    )
    tampered = verify_run.verify_docking_target(tmp_path, target_id=target_id, expected_top_n=2)
    assert tampered["ok"] is False
    assert "successful_docking_fusion_mismatch" in tampered["failures"]

    fused.loc[0, "fused_score"] = 0.9
    fused.to_parquet(dock_dir / "fused_scores.parquet", index=False)
    (dock_dir / "metrics_before_after.json").write_text(
        json.dumps({"before": before, "after": after, "weight": 0.3}), encoding="utf-8"
    )

    selected.loc[1, "candidate_index"] = 3
    selected.to_parquet(dock_dir / "selected_results.parquet", index=False)
    assert verify_run.verify_docking_target(tmp_path, target_id=target_id, expected_top_n=2)["ok"] is False


def test_audit_report_states_evaluation_boundaries():
    from scientific_validation.four_level_cli_1kx10k import report

    markdown = report.render_markdown(
        {
            "n_targets": 1000,
            "n_pairs": 10_000_000,
            "tier_counts": {"strict_labeled": 167, "limited_positive": 658, "unlabeled_fallback": 175},
            "fit_pair_contamination": 0,
            "layer_failures": {"l1": 0, "l2": 0, "l3": 0, "l4": 0},
            "aggregates": {"recall_at_1pct_l2": {"n": 825, "mean": 0.9, "median": 0.8, "q1": 0.7, "q3": 1.0, "median_bootstrap_95ci": [0.7, 0.9]}},
            "strict_aggregates": {"recall_at_1pct_l2": {"n": 167, "mean": 0.1, "median": 0.0, "q1": 0.0, "q3": 0.1, "median_bootstrap_95ci": [0.0, 0.1]}},
            "docking": {
                "target_id": "CHEMBL2051",
                "n_submitted": 300,
                "status_counts": {"ok": 157, "skipped_hac": 82, "dock_timeout": 39, "dock_no_score": 21, "prep_timeout": 1},
                "before": {"auc_l2_labeled": 0.844444, "auc_final4_labeled": 0.8},
                "after": {"auc_fused_labeled": 0.844444},
            },
        }
    )

    assert "pair-heldout" in markdown
    assert "PU" in markdown
    assert "10,000,000" in markdown
    assert "CHEMBL2051" in markdown
    assert "0.935" in markdown
    assert "不是本次通用能力结论" in markdown
    assert "成功 157/300" in markdown
    assert "0.844444" in markdown
    assert "均值 0.1" in markdown
    assert "均值 0.9" not in markdown
