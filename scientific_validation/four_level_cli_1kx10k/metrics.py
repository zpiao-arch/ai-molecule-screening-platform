from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score


POSITIVE_ROLE = "heldout_positive"
NEGATIVE_ROLE = "heldout_negative"


def _strict_auc(frame: pd.DataFrame, score_column: str) -> float | None:
    labeled = frame.loc[frame["label_role"].isin([POSITIVE_ROLE, NEGATIVE_ROLE])]
    if labeled.empty:
        return None
    labels = labeled["label_role"].eq(POSITIVE_ROLE).astype(int).to_numpy()
    if len(np.unique(labels)) != 2:
        return None
    scores = pd.to_numeric(labeled[score_column], errors="coerce")
    valid = scores.notna().to_numpy()
    if len(np.unique(labels[valid])) != 2:
        return None
    return round(float(roc_auc_score(labels[valid], scores.to_numpy()[valid])), 6)


def _early_retrieval(
    frame: pd.DataFrame,
    score_column: str,
    *,
    top_k: int,
    n_positive: int,
) -> tuple[float | None, float | None]:
    if n_positive == 0:
        return None, None
    scores = pd.to_numeric(frame[score_column], errors="coerce").to_numpy()
    sortable = np.where(np.isfinite(scores), scores, -np.inf)
    order = np.argsort(-sortable, kind="mergesort")[:top_k]
    roles = frame["label_role"].to_numpy()
    hits = int(np.sum(roles[order] == POSITIVE_ROLE))
    recall = hits / n_positive
    observed_rate = hits / top_k
    baseline = n_positive / len(frame)
    enrichment = observed_rate / baseline if baseline else None
    return round(float(recall), 6), round(float(enrichment), 6) if enrichment is not None else None


def compute_target_metrics(frame: pd.DataFrame, top_fraction: float = 0.01) -> dict[str, object]:
    if frame.empty:
        raise ValueError("cannot compute metrics for an empty target pool")
    if not 0 < top_fraction <= 1:
        raise ValueError("top_fraction must be in (0, 1]")
    required = {"label_role", "l2", "final_score"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"score frame missing columns: {sorted(missing)}")

    n_positive = int(frame["label_role"].eq(POSITIVE_ROLE).sum())
    n_negative = int(frame["label_role"].eq(NEGATIVE_ROLE).sum())
    top_k = max(1, int(round(len(frame) * top_fraction)))
    result: dict[str, object] = {
        "n_pool": int(len(frame)),
        "n_labeled_positive": n_positive,
        "n_labeled_negative": n_negative,
        "n_unlabeled_background": int(frame["label_role"].eq("unlabeled_background").sum()),
        "top_k_1pct": top_k,
        "auc_l2_labeled": _strict_auc(frame, "l2"),
        "auc_final4_labeled": _strict_auc(frame, "final_score"),
    }
    for column, suffix in (("l2", "l2"), ("final_score", "final4")):
        recall, enrichment = _early_retrieval(
            frame, column, top_k=top_k, n_positive=n_positive
        )
        result[f"recall_at_1pct_{suffix}"] = recall
        result[f"observed_enrichment_at_1pct_{suffix}"] = enrichment

    if "fused_score" in frame:
        result["auc_fused_labeled"] = _strict_auc(frame, "fused_score")
        recall, enrichment = _early_retrieval(
            frame, "fused_score", top_k=top_k, n_positive=n_positive
        )
        result["recall_at_1pct_fused"] = recall
        result["observed_enrichment_at_1pct_fused"] = enrichment

    for status_column in ("l1_status", "layer2_status", "l3_status", "l4_status"):
        if status_column in frame:
            result[f"{status_column}_failure_rate"] = round(
                float(1.0 - frame[status_column].eq("ok").mean()), 8
            )
    if "gate_status" in frame:
        result["gate_pass_rate"] = round(float(frame["gate_status"].eq("PASS").mean()), 8)
    return result


def _summary(values: Iterable[float], bootstrap_reps: int, rng: np.random.RandomState) -> dict[str, object]:
    array = np.asarray(list(values), dtype=float)
    array = array[np.isfinite(array)]
    if not len(array):
        return {"n": 0, "mean": None, "median": None, "q1": None, "q3": None,
                "median_bootstrap_95ci": [None, None]}
    bootstrap = np.empty(bootstrap_reps, dtype=float)
    for index in range(bootstrap_reps):
        bootstrap[index] = np.median(rng.choice(array, size=len(array), replace=True))
    return {
        "n": int(len(array)),
        "mean": round(float(np.mean(array)), 6),
        "median": round(float(np.median(array)), 6),
        "q1": round(float(np.quantile(array, 0.25)), 6),
        "q3": round(float(np.quantile(array, 0.75)), 6),
        "median_bootstrap_95ci": [
            round(float(np.quantile(bootstrap, 0.025)), 6),
            round(float(np.quantile(bootstrap, 0.975)), 6),
        ],
    }


def aggregate_metrics(
    per_target: pd.DataFrame,
    *,
    bootstrap_reps: int = 2000,
    seed: int = 42,
) -> dict[str, dict[str, object]]:
    if bootstrap_reps < 1:
        raise ValueError("bootstrap_reps must be positive")
    rng = np.random.RandomState(seed)
    metric_columns = [
        column
        for column in per_target.columns
        if column.startswith(("auc_", "recall_at_", "observed_enrichment_at_", "gate_pass_rate"))
    ]
    return {
        column: _summary(pd.to_numeric(per_target[column], errors="coerce").dropna(), bootstrap_reps, rng)
        for column in metric_columns
    }
