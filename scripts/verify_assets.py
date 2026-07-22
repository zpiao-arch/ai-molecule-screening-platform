from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify(root: Path, manifest_path: Path) -> dict[str, object]:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    missing = []
    mismatches = []
    for relative, expected in payload.get("files", {}).items():
        path = root / relative
        if not path.is_file():
            missing.append(relative)
        elif sha256(path) != expected:
            mismatches.append(relative)
    return {
        "ok": not missing and not mismatches,
        "missing": missing,
        "mismatches": mismatches,
        "checked": len(payload.get("files", {})),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify external assets for the four-level CLI")
    parser.add_argument("--asset-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=Path("assets/ASSET_MANIFEST.json"))
    args = parser.parse_args()
    result = verify(args.asset_root.resolve(), args.manifest.resolve())
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
