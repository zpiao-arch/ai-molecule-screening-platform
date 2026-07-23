#!/usr/bin/env python3
"""Structured Python boundary for Open Molecule Lab CSV operations."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


SOURCE_ROOT = Path(__file__).resolve().parents[3]
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))


def _error_field(path: Path, message: str) -> str:
    if "表头" in message:
        return "columns"
    if "重复 id" in message:
        return "id"
    if "不能为空" in message:
        try:
            with path.open("r", encoding="utf-8-sig", newline="") as handle:
                for row in csv.DictReader(handle):
                    if not str(row.get("id") or "").strip():
                        return "id"
                    if not str(row.get("smiles") or "").strip():
                        return "smiles"
        except OSError:
            pass
    return "csvText"


def validate_molecule_set(path: Path, max_rows: int) -> dict[str, object]:
    from scoring.scoring import read_molecules_csv

    try:
        molecules = read_molecules_csv(path)
    except (OSError, UnicodeError, ValueError) as exc:
        return {
            "ok": False,
            "error": str(exc),
            "field": _error_field(path, str(exc)),
        }
    if not molecules:
        return {
            "ok": False,
            "error": "输入 CSV 不包含分子行",
            "field": "rows",
        }
    if len(molecules) > max_rows:
        return {
            "ok": False,
            "error": f"分子行数超过上限 {max_rows}",
            "field": "rows",
        }
    return {
        "ok": True,
        "nRows": len(molecules),
        "columns": ["id", "smiles"],
        "sample": [{"id": mol_id, "smiles": smiles} for mol_id, smiles in molecules[:3]],
    }


def _json_safe(value):
    import math

    if value is None:
        return None
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if hasattr(value, "item"):
        return _json_safe(value.item())
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    return value


def summarize_results(path: Path) -> dict[str, object]:
    import pandas as pd

    frame = pd.read_csv(path)
    status_columns = [f"layer{layer}_status" for layer in range(1, 5)]
    missing = [column for column in status_columns + ["final_score"] if column not in frame.columns]
    if missing:
        return {"ok": False, "error": f"results missing required columns: {', '.join(missing)}", "field": "results"}
    valid = frame[status_columns].eq("ok").all(axis=1) & frame["final_score"].notna()
    ranked = frame.loc[valid].sort_values("final_score", ascending=False, kind="mergesort")
    failed = frame.loc[~valid]
    return {
        "ok": True,
        "nRows": int(len(frame)),
        "nRanked": int(valid.sum()),
        "nFailed": int((~valid).sum()),
        "columns": [str(column) for column in frame.columns],
        "rows": _json_safe(frame.to_dict(orient="records")),
        "rankedRows": _json_safe(ranked.head(100).to_dict(orient="records")),
        "failedRows": _json_safe(failed.head(100).to_dict(orient="records")),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate = subparsers.add_parser("validate-molecule-set")
    validate.add_argument("--input", type=Path, required=True)
    validate.add_argument("--max-rows", type=int, default=100_000)
    summarize = subparsers.add_parser("summarize-results")
    summarize.add_argument("--input", type=Path, required=True)
    args = parser.parse_args()

    if args.command == "validate-molecule-set":
        result = validate_molecule_set(args.input.resolve(), args.max_rows)
    elif args.command == "summarize-results":
        result = summarize_results(args.input.resolve())
    else:  # pragma: no cover
        result = {"ok": False, "error": "unsupported command", "field": "command"}
    print(json.dumps(result, ensure_ascii=False))
    return 0 if result.get("ok") else 2


if __name__ == "__main__":
    raise SystemExit(main())
