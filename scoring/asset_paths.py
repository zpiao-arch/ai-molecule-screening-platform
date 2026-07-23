"""Resolve externally distributed model and receptor assets without symlinks."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


ASSET_ROOT_ENV = "FOUR_LEVEL_ASSET_ROOT"


@dataclass(frozen=True)
class AssetPaths:
    root: Path
    manifest: Path
    model_root: Path
    l2_model: Path
    l2_params: Path
    admet_model_dir: Path
    unimol_model_dir: Path
    receptor_root: Path


def resolve_asset_paths(asset_root: str | Path | None = None) -> AssetPaths | None:
    configured = asset_root if asset_root is not None else os.environ.get(ASSET_ROOT_ENV)
    if configured is None or not str(configured).strip():
        return None

    root = Path(configured).expanduser().resolve()
    model_root = (root / "scoring" / "models").resolve()
    l2_dir = model_root / "bindingdb_l2"
    compatible_model = l2_dir / "l2_model_sklearn_1_7_2.joblib"
    l2_model = compatible_model if compatible_model.is_file() else l2_dir / "l2_model.joblib"
    return AssetPaths(
        root=root,
        manifest=(root / "ASSET_MANIFEST.json").resolve(),
        model_root=model_root,
        l2_model=l2_model.resolve(),
        l2_params=(l2_dir / "l2_params.json").resolve(),
        admet_model_dir=(model_root / "admet").resolve(),
        unimol_model_dir=model_root,
        receptor_root=(root / "scoring" / "receptors").resolve(),
    )
