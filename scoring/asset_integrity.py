"""SHA-256 trust gate for external model and data assets."""

from __future__ import annotations

import hashlib
import json
import os
from functools import lru_cache
from pathlib import Path
from typing import Any


UNSAFE_ENV = "FOUR_LEVEL_ALLOW_UNVERIFIED_ASSETS"
MANIFEST_ENV = "FOUR_LEVEL_ASSET_MANIFEST"
_TRUE_VALUES = {"1", "true", "yes", "on"}


class AssetIntegrityError(RuntimeError):
    """Raised before an untrusted external asset can be deserialized."""


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def unsafe_override_enabled() -> bool:
    return os.environ.get(UNSAFE_ENV, "").strip().lower() in _TRUE_VALUES


def _manifest_root(manifest_path: Path) -> Path:
    return manifest_path.parent.parent if manifest_path.parent.name == "assets" else manifest_path.parent


def discover_manifest(asset_path: str | Path) -> Path | None:
    configured = os.environ.get(MANIFEST_ENV)
    if configured:
        return Path(configured).expanduser().resolve()

    path = Path(asset_path).expanduser().resolve()
    for parent in (path.parent, *path.parents):
        direct = parent / "ASSET_MANIFEST.json"
        if direct.is_file():
            return direct.resolve()
        source_manifest = parent / "assets" / "ASSET_MANIFEST.json"
        if source_manifest.is_file():
            return source_manifest.resolve()
    return None


@lru_cache(maxsize=256)
def _cached_sha256(
    resolved_path: str,
    size: int,
    mtime_ns: int,
    manifest_path: str,
    expected: str,
) -> str:
    del size, mtime_ns, manifest_path, expected
    return sha256_file(resolved_path)


def _load_manifest(manifest_path: Path) -> dict[str, Any]:
    if not manifest_path.is_file():
        raise AssetIntegrityError(f"asset manifest not found: {manifest_path}")
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AssetIntegrityError(f"asset manifest is unreadable: {manifest_path}: {exc}") from exc
    files = payload.get("files") if isinstance(payload, dict) else None
    if not isinstance(files, dict) or not files:
        raise AssetIntegrityError(f"asset manifest has no files mapping: {manifest_path}")
    return payload


def _unsafe_result(path: Path, manifest_path: Path | None, reason: str) -> dict[str, object]:
    return {
        "path": str(path),
        "manifest": str(manifest_path) if manifest_path is not None else None,
        "relative_path": None,
        "sha256": sha256_file(path) if path.is_file() else None,
        "verified": False,
        "unsafe_override": True,
        "reason": reason,
    }


def verify_asset(
    path: str | Path,
    *,
    manifest_path: str | Path | None = None,
) -> dict[str, object]:
    asset = Path(path).expanduser().resolve()
    if not asset.is_file():
        raise AssetIntegrityError(f"asset file not found: {asset}")

    manifest = (
        Path(manifest_path).expanduser().resolve()
        if manifest_path is not None
        else discover_manifest(asset)
    )
    try:
        if manifest is None:
            raise AssetIntegrityError(f"asset manifest not found for: {asset}")
        payload = _load_manifest(manifest)
        root = _manifest_root(manifest)
        try:
            relative = asset.relative_to(root).as_posix()
        except ValueError as exc:
            raise AssetIntegrityError(
                f"asset is outside manifest root and not listed: {asset}"
            ) from exc
        expected = payload["files"].get(relative)
        if not isinstance(expected, str):
            raise AssetIntegrityError(f"asset is not listed in manifest: {relative}")
        stat = asset.stat()
        actual = _cached_sha256(
            str(asset),
            stat.st_size,
            stat.st_mtime_ns,
            str(manifest),
            expected,
        )
        if actual != expected:
            raise AssetIntegrityError(
                f"asset sha256 mismatch: {relative}: expected {expected}, got {actual}"
            )
        return {
            "path": str(asset),
            "manifest": str(manifest),
            "relative_path": relative,
            "sha256": actual,
            "verified": True,
            "unsafe_override": False,
            "reason": "ok",
        }
    except AssetIntegrityError as exc:
        if unsafe_override_enabled():
            return _unsafe_result(asset, manifest, str(exc))
        raise


def verify_manifest(root: str | Path, manifest_path: str | Path) -> dict[str, object]:
    asset_root = Path(root).expanduser().resolve()
    manifest = Path(manifest_path).expanduser().resolve()
    try:
        payload = _load_manifest(manifest)
    except AssetIntegrityError as exc:
        return {
            "ok": False,
            "manifest": str(manifest),
            "root": str(asset_root),
            "checked": 0,
            "missing": [],
            "mismatches": [],
            "error": str(exc),
            "unsafe_override": unsafe_override_enabled(),
        }

    missing: list[str] = []
    mismatches: list[str] = []
    for relative, expected in payload["files"].items():
        path = asset_root / relative
        if not path.is_file():
            missing.append(relative)
        elif sha256_file(path) != expected:
            mismatches.append(relative)
    return {
        "ok": not missing and not mismatches,
        "manifest": str(manifest),
        "root": str(asset_root),
        "checked": len(payload["files"]),
        "missing": missing,
        "mismatches": mismatches,
        "error": None,
        "unsafe_override": unsafe_override_enabled(),
    }
