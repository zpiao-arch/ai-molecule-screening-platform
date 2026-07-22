from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import stat
import tarfile
import zipfile
from pathlib import Path, PurePosixPath


FIXED_TIMESTAMP = (2026, 7, 21, 0, 0, 0)
EXCLUDED_DIRECTORY_NAMES = {
    ".git",
    ".pytest_cache",
    ".venv",
    "__pycache__",
    "build",
    "dist",
    "logs",
    "node_modules",
    "outputs",
    "runs",
}
EXCLUDED_SUFFIXES = {
    ".joblib",
    ".log",
    ".npz",
    ".pdbqt",
    ".pickle",
    ".pkl",
    ".pt",
    ".pyc",
    ".pyo",
    ".tsbuildinfo",
    ".zip",
}
ALLOWED_EXTERNALIZED_SUFFIX_PATHS = {
    PurePosixPath("validation/frozen_run/per_target_metrics.parquet"),
}


def _is_excluded(relative: Path) -> bool:
    posix = PurePosixPath(relative.as_posix())
    if any(part in EXCLUDED_DIRECTORY_NAMES for part in posix.parts):
        return True
    if any(part.endswith(".egg-info") for part in posix.parts):
        return True
    if posix.parts and posix.parts[0] == "data":
        return True
    if posix.parts[:2] in (("scoring", "models"), ("scoring", "receptors")):
        return True
    if posix.parts and posix.parts[0] == "legacy" and posix.suffix in {".csv", ".json", ".parquet"}:
        return True
    if posix.name == ".DS_Store" or posix.name.endswith(".tar.gz"):
        return True
    if posix.suffix == ".parquet" and posix not in ALLOWED_EXTERNALIZED_SUFFIX_PATHS:
        return True
    return posix.suffix in EXCLUDED_SUFFIXES


def _is_asset_excluded(relative: Path) -> bool:
    posix = PurePosixPath(relative.as_posix())
    if any(part in EXCLUDED_DIRECTORY_NAMES for part in posix.parts):
        return True
    if posix.name == ".DS_Store" or posix.name.startswith("._"):
        return True
    return posix.suffix in {".log", ".pyc", ".pyo"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_source_zip(source_dir: str | Path, output: str | Path) -> dict[str, object]:
    source = Path(source_dir).resolve()
    destination = Path(output).resolve()
    if not source.is_dir():
        raise ValueError(f"source directory does not exist: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)

    files = sorted(
        path
        for path in source.rglob("*")
        if path.is_file() and not _is_excluded(path.relative_to(source))
    )
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with zipfile.ZipFile(
        temporary,
        mode="w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=9,
    ) as archive:
        for path in files:
            relative = path.relative_to(source)
            archive_name = (PurePosixPath(source.name) / PurePosixPath(relative.as_posix())).as_posix()
            info = zipfile.ZipInfo(archive_name, date_time=FIXED_TIMESTAMP)
            info.create_system = 3
            mode = stat.S_IMODE(path.stat().st_mode)
            info.external_attr = (stat.S_IFREG | mode) << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    temporary.replace(destination)
    return {
        "output": str(destination),
        "sha256": _sha256(destination),
        "compressed_bytes": destination.stat().st_size,
        "source_files": len(files),
        "top_level_directory": source.name,
    }


def build_asset_tar(asset_dir: str | Path, output: str | Path) -> dict[str, object]:
    source = Path(asset_dir).resolve()
    destination = Path(output).resolve()
    if not source.is_dir():
        raise ValueError(f"asset directory does not exist: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)

    files = sorted(
        path
        for path in source.rglob("*")
        if path.is_file() and not _is_asset_excluded(path.relative_to(source))
    )
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    with temporary.open("wb") as raw_handle:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_handle,
            compresslevel=9,
            mtime=0,
        ) as gzip_handle:
            with tarfile.open(fileobj=gzip_handle, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for path in files:
                    relative = path.relative_to(source)
                    archive_name = (
                        PurePosixPath(source.name) / PurePosixPath(relative.as_posix())
                    ).as_posix()
                    info = tarfile.TarInfo(archive_name)
                    info.size = path.stat().st_size
                    info.mode = stat.S_IMODE(path.stat().st_mode)
                    info.mtime = 0
                    info.uid = 0
                    info.gid = 0
                    info.uname = ""
                    info.gname = ""
                    with path.open("rb") as handle:
                        archive.addfile(info, handle)
    temporary.replace(destination)
    return {
        "output": str(destination),
        "sha256": _sha256(destination),
        "compressed_bytes": destination.stat().st_size,
        "asset_files": len(files),
        "top_level_directory": source.name,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build deterministic source and asset archives")
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--asset-dir", type=Path)
    parser.add_argument("--asset-output", type=Path)
    args = parser.parse_args(argv)
    if not args.output and not args.asset_output:
        parser.error("至少需要 --output 或 --asset-output")
    if bool(args.asset_dir) != bool(args.asset_output):
        parser.error("--asset-dir 与 --asset-output 必须同时提供")
    result = {}
    if args.output:
        result["source"] = build_source_zip(args.source_dir, args.output)
    if args.asset_output:
        result["assets"] = build_asset_tar(args.asset_dir, args.asset_output)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
