import numpy as np
import pandas as pd

from scientific_validation.four_level_cli_1kx10k import metrics


def _frame_with_background():
    rows = []
    for i in range(5):
        rows.append({
            "canonical_smiles": f"P{i}",
            "label_role": "heldout_positive",
            "l2": 0.9 - i * 0.01,
            "final_score": 0.85 - i * 0.01,
            "fused_score": 0.2 + i * 0.01,
        })
    for i in range(5):
        rows.append({
            "canonical_smiles": f"N{i}",
            "label_role": "heldout_negative",
            "l2": 0.1 + i * 0.01,
            "final_score": 0.15 + i * 0.01,
            "fused_score": 0.8 - i * 0.01,
        })
    for i in range(90):
        rows.append({
            "canonical_smiles": f"U{i}",
            "label_role": "unlabeled_background",
            "l2": 0.99 if i < 2 else 0.3,
            "final_score": 0.98 if i < 2 else 0.35,
            "fused_score": 0.4,
        })
    return pd.DataFrame(rows)


def test_labeled_auc_excludes_unlabeled_background():
    result = metrics.compute_target_metrics(_frame_with_background(), top_fraction=0.01)

    assert result["auc_l2_labeled"] == 1.0
    assert result["auc_final4_labeled"] == 1.0
    assert result["n_labeled_positive"] == 5
    assert result["n_labeled_negative"] == 5


def test_top_one_percent_recall_and_pu_enrichment_use_full_pool():
    frame = _frame_with_background()
    frame.loc[frame["canonical_smiles"] == "P0", "l2"] = 1.0

    result = metrics.compute_target_metrics(frame, top_fraction=0.01)

    assert result["top_k_1pct"] == 1
    assert result["recall_at_1pct_l2"] == 0.2
    assert result["observed_enrichment_at_1pct_l2"] == 20.0


def test_fused_metrics_are_recomputed_from_fused_scores():
    result = metrics.compute_target_metrics(_frame_with_background(), top_fraction=0.01)

    assert result["auc_fused_labeled"] == 0.0
    assert result["auc_fused_labeled"] != result["auc_l2_labeled"]


def test_no_labeled_negative_has_no_strict_auc():
    frame = _frame_with_background()
    frame = frame.loc[frame["label_role"] != "heldout_negative"]

    result = metrics.compute_target_metrics(frame)

    assert result["auc_l2_labeled"] is None
    assert result["auc_final4_labeled"] is None


def test_bootstrap_summary_is_deterministic():
    per_target = pd.DataFrame(
        {"auc_l2_labeled": np.linspace(0.5, 0.9, 20), "auc_final4_labeled": np.linspace(0.55, 0.95, 20)}
    )

    first = metrics.aggregate_metrics(per_target, bootstrap_reps=200, seed=42)
    second = metrics.aggregate_metrics(per_target, bootstrap_reps=200, seed=42)

    assert first == second
    assert first["auc_l2_labeled"]["n"] == 20
    assert len(first["auc_l2_labeled"]["median_bootstrap_95ci"]) == 2
