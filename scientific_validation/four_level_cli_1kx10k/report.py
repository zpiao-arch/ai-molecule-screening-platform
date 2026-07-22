from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd

from . import metrics, verify_run
from . import provenance


def prepare_delivery_artifacts(run_dir: str | Path) -> None:
    root = Path(run_dir).resolve()
    run_config_path = root / "run_config.json"
    if run_config_path.is_file():
        run_config = json.loads(run_config_path.read_text(encoding="utf-8"))
        design_value = str(run_config.get("design", ""))
        design = Path(design_value)
        if not design.is_absolute():
            scoring_dir = Path(str(run_config.get("scoring_dir", "")))
            design = (scoring_dir.parent / design).resolve() if scoring_dir else (Path.cwd() / design).resolve()
        if design.is_file() and not (root / "DESIGN.md").is_file():
            shutil.copyfile(design, root / "DESIGN.md")
    docking_summary = {}
    docking_path = root / "docking" / "CHEMBL2051" / "docking_summary.json"
    if docking_path.is_file():
        docking_summary = json.loads(docking_path.read_text(encoding="utf-8"))
    layer_summary = {}
    layer_path = root / "layer_cache_summary.json"
    if layer_path.is_file():
        layer_summary = json.loads(layer_path.read_text(encoding="utf-8"))
    failure_manifest = {
        "layer_failures": {
            "l1": layer_summary.get("l1_failures"),
            "l3": layer_summary.get("l3_failures"),
            "l4": layer_summary.get("l4_failures"),
        },
        "docking_status_counts": docking_summary.get("status_counts", {}),
        "interpretation": {
            "unlabeled_background_is_not_confirmed_negative": True,
            "pair_heldout_is_not_cold_start": True,
            "docking_failures_are_explicit": True,
        },
    }
    (root / "failure_manifest.json").write_text(
        json.dumps(failure_manifest, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )


def _metric_line(values: dict, name: str) -> str:
    value = values.get(name, {}) if isinstance(values, dict) else {}
    if not value or value.get("n", 0) == 0:
        return "无可报告值"
    ci = value.get("median_bootstrap_95ci", [None, None])
    return f"中位数 {value.get('median')}，均值 {value.get('mean')}，IQR [{value.get('q1')}, {value.get('q3')}]，bootstrap 95% CI [{ci[0]}, {ci[1]}]（n={value.get('n')}）"


def render_markdown(summary: dict[str, object]) -> str:
    aggregates = summary.get("aggregates", {})
    strict = summary.get("strict_aggregates", aggregates)
    docking = summary.get("docking") or {}
    tier_counts = summary.get("tier_counts", {})
    failures = summary.get("layer_failures", {})
    n_pairs = int(summary.get("n_pairs", 0))
    docking_statuses = docking.get("status_counts", {}) if isinstance(docking, dict) else {}
    docking_ok = int(docking_statuses.get("ok", 0))
    docking_submitted = int(docking.get("n_submitted", 0) or 0) if isinstance(docking, dict) else 0
    docking_before = docking.get("before", {}) if isinstance(docking, dict) else {}
    docking_after = docking.get("after", {}) if isinstance(docking, dict) else {}
    return f"""# 四级 CLI 1000x10000 闭环审计报告

## 结论范围

本运行完成工程规模 **{int(summary.get('n_targets', 0)):,} 靶点 × 10,000 候选 = {n_pairs:,} 对**。工程吞吐与科研效果分开报告。效果标签采用 **pair-heldout**：评测对不出现在模型拟合对中，但不等同于冷靶点、冷分子、骨架或时间外推。

随机背景是 **PU（positive-unlabeled）未标注背景**，不是确认阴性；严格 AUC 只使用留出确认正例与留出确认阴性。靶点层级：{json.dumps(tier_counts, ensure_ascii=False, sort_keys=True)}。

## 严格留出效果

- L2 AUC：{_metric_line(strict, 'auc_l2_labeled')}
- 四级综合 AUC：{_metric_line(strict, 'auc_final4_labeled')}
- L2 top-1% recall：{_metric_line(strict, 'recall_at_1pct_l2')}
- 四级 top-1% recall：{_metric_line(strict, 'recall_at_1pct_final4')}
- PU 观察富集（L2 top-1%）：{_metric_line(strict, 'observed_enrichment_at_1pct_l2')}

上述 AUC 的拟合对直接污染计数为 **{summary.get('fit_pair_contamination', '未记录')}**。层失败计数：{json.dumps(failures, ensure_ascii=False, sort_keys=True)}。

## 受体级联

本地注册表当前只有 **{docking.get('target_id', '无')}** 的有效受体；本次只对其 L2 top-{docking.get('n_submitted', '未运行')} 运行真实 obabel+smina，并按逐分子状态重新计算融合指标。其余靶点没有 3D 对接结论。

- 对接成功 {docking_ok}/{docking_submitted}；逐状态计数：{json.dumps(docking_statuses, ensure_ascii=False, sort_keys=True)}。
- 对接前 L2 AUC：{docking_before.get('auc_l2_labeled')}；对接前四级 AUC：{docking_before.get('auc_final4_labeled')}；融合后 AUC：{docking_after.get('auc_fused_labeled')}。

CHEMBL2051 只有 {docking_after.get('n_labeled_positive', docking_before.get('n_labeled_positive', '未知'))} 个确认正例和 {docking_after.get('n_labeled_negative', docking_before.get('n_labeled_negative', '未知'))} 个确认阴性，属于 limited-positive 层级，不进入 167 靶点严格聚合。

历史材料中的 **0.935** 是人工平衡集/特定级联展示值，**不是本次通用能力结论**；本报告不把它替换为 1000 靶点泛化性能。

## 复算与限制

报告由候选级 `scores/` 分区和运行配置重新计算；`MANIFEST.sha256`、环境/资产哈希、分片 checkpoint 与对接命令参数用于复核。该 pair-heldout 结果仍可能受同靶点/同分子跨划分影响，不能替代冷启动、时间外或湿实验验证。
"""


def recompute_summary(run_dir: str | Path) -> dict[str, object]:
    root = Path(run_dir).resolve()
    verification = verify_run.verify_artifacts(
        root,
        strict=True,
        require_docking=True,
    )
    if not verification["ok"]:
        raise RuntimeError(f"run verification failed: {verification['failures']}")
    targets = pd.read_parquet(root / "target_manifest.parquet")
    metric_rows = []
    layer_failures = {"l1": 0, "l2": 0, "l3": 0, "l4": 0}
    n_pairs = 0
    for target in targets.itertuples(index=False):
        path = root / "scores" / f"target_id={target.target_id}" / "part-000000.parquet"
        frame = pd.read_parquet(path)
        row = metrics.compute_target_metrics(frame)
        row.update({"target_id": str(target.target_id), "evaluation_tier": str(target.evaluation_tier)})
        metric_rows.append(row)
        n_pairs += len(frame)
        for key, column in (("l1", "l1_status"), ("l2", "layer2_status"), ("l3", "l3_status"), ("l4", "l4_status")):
            layer_failures[key] += int((frame[column] != "ok").sum())

    metric_frame = pd.DataFrame(metric_rows)
    strict_frame = metric_frame[metric_frame["evaluation_tier"] == "strict_labeled"]
    docking = {}
    docking_dir = root / "docking" / "CHEMBL2051"
    if (docking_dir / "docking_summary.json").is_file():
        docking.update(json.loads((docking_dir / "docking_summary.json").read_text(encoding="utf-8")))
        if (docking_dir / "metrics_before_after.json").is_file():
            docking.update(json.loads((docking_dir / "metrics_before_after.json").read_text(encoding="utf-8")))

    return {
        "n_targets": int(len(targets)),
        "n_pairs": int(n_pairs),
        "tier_counts": targets["evaluation_tier"].value_counts().to_dict(),
        "fit_pair_contamination": verification["checks"]["fit_pair_contamination"],
        "layer_failures": layer_failures,
        "aggregates": metrics.aggregate_metrics(metric_frame, bootstrap_reps=2000, seed=42),
        "strict_aggregates": metrics.aggregate_metrics(strict_frame, bootstrap_reps=2000, seed=42),
        "docking": docking,
        "verification": verification,
    }


def generate_report(run_dir: str | Path) -> tuple[Path, Path]:
    root = Path(run_dir).resolve()
    prepare_delivery_artifacts(root)
    summary = recompute_summary(root)
    summary_path = root / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    report_path = root / "四级CLI_1000x10000_审计报告.md"
    report_path.write_text(render_markdown(summary), encoding="utf-8")
    report_inputs = {
        "target_manifest_sha256": verify_run._sha256(root / "target_manifest.parquet"),
        "score_summary_sha256": verify_run._sha256(root / "score_summary.json"),
        "docking_metrics_sha256": verify_run._sha256(root / "docking" / "CHEMBL2051" / "metrics_before_after.json"),
        "report_code_sha256": verify_run._sha256(Path(__file__)),
    }
    provenance.write_stage_checkpoint(
        root,
        "report",
        inputs=report_inputs,
        status="complete",
        outputs={
            "summary_sha256": verify_run._sha256(summary_path),
            "report_sha256": verify_run._sha256(report_path),
        },
    )
    verify_run.write_manifest(root)
    return summary_path, report_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate and seal the four-level audit report")
    parser.add_argument("--run-dir", required=True, type=Path)
    args = parser.parse_args(argv)
    summary_path, report_path = generate_report(args.run_dir)
    print(json.dumps({
        "run_dir": str(Path(args.run_dir).resolve()),
        "summary": str(summary_path),
        "report": str(report_path),
        "manifest": str(Path(args.run_dir).resolve() / verify_run.MANIFEST_NAME),
    }, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
