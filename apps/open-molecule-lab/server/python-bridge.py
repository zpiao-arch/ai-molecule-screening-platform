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


def _read_results(path: Path):
    import pandas as pd

    return pd.read_csv(path, dtype={"id": "string", "smiles": "string"})


def _ranking_score_field(frame) -> str:
    if "final_score_dock" in frame.columns and frame["final_score_dock"].notna().any():
        return "final_score_dock"
    return "final_score"


def _validate_result_identity(frame, expected_input: Path | None) -> dict[str, object] | None:
    if "id" not in frame.columns or "smiles" not in frame.columns:
        return {"ok": False, "error": "results missing required columns: id, smiles", "field": "results"}
    result_ids = frame["id"].fillna("").astype(str)
    if result_ids.eq("").any():
        return {"ok": False, "error": "results contain blank ids", "field": "results"}
    duplicated = sorted(result_ids[result_ids.duplicated(keep=False)].unique().tolist())
    if duplicated:
        return {
            "ok": False,
            "error": f"results contain duplicate ids: {', '.join(duplicated[:10])}",
            "field": "results",
        }
    if expected_input is None:
        return None

    expected = _read_results(expected_input)
    if "id" not in expected.columns or "smiles" not in expected.columns:
        return {"ok": False, "error": "expected input missing id/smiles columns", "field": "expectedInput"}
    expected_ids = expected["id"].fillna("").astype(str)
    result_set = set(result_ids)
    expected_set = set(expected_ids)
    missing = sorted(expected_set - result_set)
    unexpected = sorted(result_set - expected_set)
    if len(frame) != len(expected) or missing or unexpected:
        return {
            "ok": False,
            "error": (
                f"result identity mismatch: expected={len(expected)} actual={len(frame)} "
                f"missing={missing[:10]} unexpected={unexpected[:10]}"
            ),
            "field": "results",
        }
    expected_smiles = dict(zip(expected_ids, expected["smiles"].fillna("").astype(str), strict=True))
    result_smiles = dict(zip(result_ids, frame["smiles"].fillna("").astype(str), strict=True))
    changed = sorted(mol_id for mol_id in expected_set if expected_smiles[mol_id] != result_smiles[mol_id])
    if changed:
        return {
            "ok": False,
            "error": f"results changed smiles for ids: {', '.join(changed[:10])}",
            "field": "results",
        }
    return None


def summarize_results(path: Path, expected_input: Path | None = None) -> dict[str, object]:
    frame = _read_results(path)
    identity_error = _validate_result_identity(frame, expected_input)
    if identity_error:
        return identity_error
    status_columns = [f"layer{layer}_status" for layer in range(1, 5)]
    missing = [column for column in status_columns + ["final_score"] if column not in frame.columns]
    if missing:
        return {"ok": False, "error": f"results missing required columns: {', '.join(missing)}", "field": "results"}
    ranking_score_field = _ranking_score_field(frame)
    valid = frame[status_columns].eq("ok").all(axis=1) & frame[ranking_score_field].notna()
    structure_docking_ok = (
        int(frame["structure_docking_status"].eq("ok").sum())
        if "structure_docking_status" in frame.columns
        else 0
    )
    return {
        "ok": True,
        "nRows": int(len(frame)),
        "nRanked": int(valid.sum()),
        "nFailed": int((~valid).sum()),
        "columns": [str(column) for column in frame.columns],
        "rankingScoreField": ranking_score_field,
        "structureDockingOk": structure_docking_ok,
    }


def page_results(path: Path, offset: int, limit: int, view: str) -> dict[str, object]:
    import pandas as pd

    if offset < 0 or limit < 1 or limit > 200:
        return {"ok": False, "error": "offset must be >= 0 and limit must be 1..200", "field": "pagination"}
    frame = _read_results(path)
    status_columns = [f"layer{layer}_status" for layer in range(1, 5)]
    missing = [column for column in status_columns + ["final_score"] if column not in frame.columns]
    if missing:
        return {"ok": False, "error": f"results missing required columns: {', '.join(missing)}", "field": "results"}
    ranking_score_field = _ranking_score_field(frame)
    valid = frame[status_columns].eq("ok").all(axis=1) & frame[ranking_score_field].notna()
    if view == "ranked":
        selected = frame.loc[valid].sort_values(ranking_score_field, ascending=False, kind="mergesort")
    elif view == "failed":
        selected = frame.loc[~valid]
    else:
        selected = frame
    page = selected.iloc[offset:offset + limit]
    return {
        "ok": True,
        "view": view,
        "offset": offset,
        "limit": limit,
        "total": int(len(selected)),
        "rankingScoreField": ranking_score_field,
        "rows": _json_safe(page.to_dict(orient="records")),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate = subparsers.add_parser("validate-molecule-set")
    validate.add_argument("--input", type=Path, required=True)
    validate.add_argument("--max-rows", type=int, default=100_000)
    summarize = subparsers.add_parser("summarize-results")
    summarize.add_argument("--input", type=Path, required=True)
    summarize.add_argument("--expected-input", type=Path)
    page = subparsers.add_parser("page-results")
    page.add_argument("--input", type=Path, required=True)
    page.add_argument("--offset", type=int, default=0)
    page.add_argument("--limit", type=int, default=100)
    page.add_argument("--view", choices=["ranked", "failed", "all"], default="ranked")
    args = parser.parse_args()

    if args.command == "validate-molecule-set":
        result = validate_molecule_set(args.input.resolve(), args.max_rows)
    elif args.command == "summarize-results":
        result = summarize_results(
            args.input.resolve(),
            args.expected_input.resolve() if args.expected_input else None,
        )
    elif args.command == "page-results":
        result = page_results(args.input.resolve(), args.offset, args.limit, args.view)
    else:  # pragma: no cover
        result = {"ok": False, "error": "unsupported command", "field": "command"}
    print(json.dumps(result, ensure_ascii=False))
    return 0 if result.get("ok") else 2


if __name__ == "__main__":
    raise SystemExit(main())
