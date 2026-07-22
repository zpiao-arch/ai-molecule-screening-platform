from __future__ import annotations

import argparse
import hashlib
import importlib
import importlib.metadata
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any


REQUIRED_IMPORTS = (
    "numpy",
    "scipy",
    "pandas",
    "pyarrow",
    "sklearn",
    "joblib",
    "rdkit",
    "torch",
    "unimol_tools",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def asset_record(path: Path) -> dict[str, Any]:
    exists = path.is_file()
    return {
        "path": str(path.resolve()),
        "exists": exists,
        "size": path.stat().st_size if exists else None,
        "sha256": sha256_file(path) if exists else None,
    }


def binary_record(path: Path) -> dict[str, Any]:
    record = asset_record(path)
    record["executable"] = bool(record["exists"] and os.access(path, os.X_OK))
    return record


def resolve_binary(name: str, root: Path) -> Path:
    direct = os.environ.get(f"{name.upper()}_BIN")
    if direct:
        configured = Path(direct).expanduser()
        candidate = configured / name if configured.is_dir() else configured
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    dock_bin_dir = os.environ.get("DOCK_BIN_DIR")
    if dock_bin_dir:
        candidate = Path(dock_bin_dir).expanduser() / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate.resolve()
    package_candidate = root / "bin" / name
    if package_candidate.is_file() and os.access(package_candidate, os.X_OK):
        return package_candidate.resolve()
    path_candidate = shutil.which(name)
    return Path(path_candidate).resolve() if path_candidate else package_candidate.resolve()


def _load_asset_integrity(scoring_dir: Path):
    candidates = (
        scoring_dir / "asset_integrity.py",
        Path(__file__).resolve().parents[2] / "scoring" / "asset_integrity.py",
    )
    for module_path in candidates:
        if not module_path.is_file():
            continue
        spec = importlib.util.spec_from_file_location("four_level_runtime_asset_integrity", module_path)
        if spec is None or spec.loader is None:
            continue
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    return None


def import_status(module_name: str) -> str:
    try:
        importlib.import_module(module_name)
        return "ok"
    except Exception as exc:
        return f"error:{type(exc).__name__}:{exc}"


def package_version(distribution: str) -> str | None:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return None


def run_backend_probe(scoring_dir: str | Path) -> dict[str, Any]:
    scoring_dir = Path(scoring_dir).resolve()
    path = str(scoring_dir)
    if path not in sys.path:
        sys.path.insert(0, path)
    from scoring import MoleculeScorer

    model_path = scoring_dir / "models" / "bindingdb_l2" / "l2_model_sklearn_1_7_2.joblib"
    scorer = MoleculeScorer(
        use_unimol=True,
        l2_method="bindingdb",
        default_target_text="CHEMBL2051 Neuraminidase Influenza A virus",
        l2_model_path=str(model_path),
        strict_backends=True,
    )
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    row = scorer.score_one(smiles, mol_id="runtime_probe")
    statuses = {
        "l1": row.get("layer1_status"),
        "l2": row.get("layer2_status"),
        "l3": row.get("layer3_status"),
        "l4": row.get("layer4_status"),
    }
    scores = {
        "l1": row.get("layer1_score"),
        "l2": row.get("docking_normalized"),
        "l3": row.get("admet_score"),
        "l4": row.get("unimol_score"),
        "final": row.get("final_score"),
    }
    ok = all(value == "ok" for value in statuses.values()) and all(
        isinstance(value, (int, float)) for value in scores.values()
    )
    return {
        "smiles": smiles,
        "statuses": statuses,
        "scores": scores,
        "overall": "ok" if ok else "failed",
    }


def inspect_runtime(
    scoring_dir: str | Path,
    *,
    probe_backends: bool = False,
) -> dict[str, Any]:
    scoring_dir = Path(scoring_dir).resolve()
    root = scoring_dir.parent
    model_dir = scoring_dir / "models"
    unimol_snapshot = (
        model_dir
        / "unimol"
        / "models--dptech--Uni-Mol-Models"
        / "snapshots"
        / "9f19c45c718192888a1c8a1c905f69f0755ea502"
    )
    imports = {name: import_status(name) for name in REQUIRED_IMPORTS}
    admet_paths = [model_dir / "admet" / f"{name}.pkl" for name in ("tox21", "bbbp", "clintox", "sider")]
    admet_assets = [asset_record(path) for path in admet_paths]
    asset_integrity = _load_asset_integrity(scoring_dir)
    manifest_probe = model_dir / "bindingdb_l2" / "l2_model_sklearn_1_7_2.joblib"
    manifest_path = asset_integrity.discover_manifest(manifest_probe) if asset_integrity else None
    if asset_integrity is None:
        manifest_result = {
            "ok": False,
            "manifest": None,
            "checked": 0,
            "missing": [],
            "mismatches": [],
            "error": "asset integrity module unavailable",
            "unsafe_override": False,
        }
    elif manifest_path is None:
        manifest_result = {
            "ok": False,
            "manifest": None,
            "checked": 0,
            "missing": [],
            "mismatches": [],
            "error": "asset manifest unavailable",
            "unsafe_override": asset_integrity.unsafe_override_enabled(),
        }
    else:
        manifest_root = (
            manifest_path.parent.parent if manifest_path.parent.name == "assets" else manifest_path.parent
        )
        manifest_result = asset_integrity.verify_manifest(manifest_root, manifest_path)

    result = {
        "python": list(sys.version_info[:3]),
        "executable": sys.executable,
        "packages": {
            "numpy": package_version("numpy"),
            "scipy": package_version("scipy"),
            "pandas": package_version("pandas"),
            "pyarrow": package_version("pyarrow"),
            "scikit-learn": package_version("scikit-learn"),
            "joblib": package_version("joblib"),
            "rdkit": package_version("rdkit"),
            "torch": package_version("torch"),
            "unimol-tools": package_version("unimol-tools"),
        },
        "imports": imports,
        "models": {
            "l2": asset_record(model_dir / "bindingdb_l2" / "l2_model.joblib"),
            "l2_compat_sklearn_1_7_2": asset_record(
                model_dir / "bindingdb_l2" / "l2_model_sklearn_1_7_2.joblib"
            ),
            "admet": {"count": sum(item["exists"] for item in admet_assets), "assets": admet_assets},
            "unimol_weights": asset_record(unimol_snapshot / "mol_pre_all_h_220816.pt"),
            "unimol_dict": asset_record(unimol_snapshot / "mol.dict.txt"),
            "unimol_references": asset_record(model_dir / "ref_embeddings.npz"),
        },
        "asset_manifest": manifest_result,
        "binaries": {
            "smina": binary_record(resolve_binary("smina", root)),
            "obabel": binary_record(resolve_binary("obabel", root)),
        },
    }
    required_ok = (
        result["python"][:2] == [3, 11]
        and all(status == "ok" for status in imports.values())
        and result["models"]["l2"]["exists"]
        and result["models"]["l2_compat_sklearn_1_7_2"]["exists"]
        and result["models"]["admet"]["count"] == 4
        and result["models"]["unimol_weights"]["exists"]
        and result["models"]["unimol_dict"]["exists"]
        and result["models"]["unimol_references"]["exists"]
        and result["asset_manifest"]["ok"]
        and not result["asset_manifest"].get("unsafe_override", False)
        and all(item["exists"] and item["executable"] for item in result["binaries"].values())
    )
    if probe_backends:
        try:
            result["backend_probe"] = run_backend_probe(scoring_dir)
        except Exception as exc:
            result["backend_probe"] = {
                "overall": "failed",
                "error": f"{type(exc).__name__}:{exc}",
            }
        required_ok = required_ok and result["backend_probe"].get("overall") == "ok"
    result["overall"] = "ok" if required_ok else "failed"
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the frozen four-level benchmark runtime")
    parser.add_argument(
        "--scoring-dir",
        default=Path(os.environ.get(
            "FOUR_LEVEL_SCORING_DIR",
            Path(__file__).resolve().parents[2] / "scoring",
        )),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    result = inspect_runtime(args.scoring_dir, probe_backends=args.strict)
    rendered = json.dumps(result, indent=2, ensure_ascii=False)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 1 if args.strict and result["overall"] != "ok" else 0


if __name__ == "__main__":
    raise SystemExit(main())
