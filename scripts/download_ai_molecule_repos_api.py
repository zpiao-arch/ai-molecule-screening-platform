#!/usr/bin/env python3
"""Mirror selected AI molecular-design repositories via the GitHub API.

This intentionally skips large model/data artifacts so the local copy remains
usable for code review and integration when normal git clone is unreliable.
"""

from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path


DEST = Path("./data/external_repos")

REPOS = [
    (
        "MolecularAI",
        "REINVENT4",
        "main",
        "REINVENT4",
        ["", "configs", "reinvent", "contrib/notebooks", "contrib/reinvent-doc/example_cfgs"],
    ),
    ("MolecularAI", "DockStream", "master", "DockStream", ["", "dockstream", "examples", "interfaces"]),
    ("MolecularAI", "GraphINVENT", "master", "GraphINVENT", ["", "graphinvent", "environments", "tools", "tutorials"]),
    ("CDDLeiden", "DrugEx", "master", "DrugEx", ["", "drugex", "docs", "tutorial/CLI", "tutorial/advanced"]),
]

SKIP_PARTS = {
    ".git",
    "__pycache__",
    ".env",
    ".ipynb_checkpoints",
    "outputs",
    "models",
    "checkpoints",
}

SKIP_SUFFIXES = {
    ".as",
    ".ckpt",
    ".gz",
    ".h5",
    ".hdf5",
    ".joblib",
    ".oeb",
    ".pkl",
    ".pickle",
    ".pt",
    ".pth",
    ".tar",
    ".tgz",
    ".zip",
}

TEXT_OR_SMALL_SUFFIXES = {
    "",
    ".cfg",
    ".conf",
    ".csv",
    ".dockerfile",
    ".ini",
    ".ipynb",
    ".json",
    ".license",
    ".md",
    ".py",
    ".rst",
    ".sh",
    ".smi",
    ".smiles",
    ".svg",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}

MAX_FILE_SIZE = 1_500_000
MAX_NOTEBOOK_SIZE = 800_000


def curl_bytes(url: str) -> bytes:
    last_error = ""
    for attempt in range(6):
        proc = subprocess.run(
            [
                "curl",
                "-fL",
                "--connect-timeout",
                "20",
                "--max-time",
                "120",
                "-H",
                "Accept: application/vnd.github+json",
                "-H",
                "User-Agent: codex-ai-molecule-repo-downloader",
                url,
            ],
            capture_output=True,
            check=False,
        )
        if proc.returncode == 0:
            return proc.stdout
        last_error = proc.stderr.decode("utf-8", errors="replace")
        time.sleep(2**attempt)
    raise RuntimeError(f"curl failed for {url}: {last_error}")


def request_json(url: str) -> dict:
    return json.loads(curl_bytes(url).decode("utf-8"))


def should_keep(path: str, size: int) -> bool:
    parts = set(Path(path).parts)
    if parts & SKIP_PARTS:
        return False

    lower = path.lower()
    suffix = Path(lower).suffix
    if suffix in SKIP_SUFFIXES:
        return False
    if suffix == ".ipynb" and size > MAX_NOTEBOOK_SIZE:
        return False
    if size > MAX_FILE_SIZE:
        return False
    if suffix in TEXT_OR_SMALL_SUFFIXES:
        return True

    # Keep setup-like files without extensions.
    name = Path(lower).name
    return name in {
        "dockerfile",
        "license",
        "makefile",
        "manifest.in",
        "requirements",
    }


def raw_url(owner: str, repo: str, branch: str, path: str) -> str:
    quoted = "/".join(part.replace("#", "%23").replace("?", "%3F") for part in Path(path).parts)
    return f"https://raw.githubusercontent.com/{owner}/{repo}/{branch}/{quoted}"


def contents_url(owner: str, repo: str, branch: str, path: str) -> str:
    if path:
        quoted = "/".join(part.replace("#", "%23").replace("?", "%3F") for part in Path(path).parts)
        return f"https://api.github.com/repos/{owner}/{repo}/contents/{quoted}?ref={branch}"
    return f"https://api.github.com/repos/{owner}/{repo}/contents?ref={branch}"


def write_readme(dest: Path, owner: str, repo: str, branch: str, kept: int, skipped: int) -> None:
    marker = dest / "_DOWNLOAD_NOTE.md"
    marker.write_text(
        "\n".join(
            [
                f"# {repo} local source mirror",
                "",
                f"Source: https://github.com/{owner}/{repo}",
                f"Branch: `{branch}`",
                "",
                "Downloaded through the GitHub API because direct `git clone` was unstable in this environment.",
                "Large model weights, binary datasets, archives, and generated output were intentionally skipped.",
                "",
                f"Files downloaded: {kept}",
                f"Files skipped: {skipped}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def iter_contents(owner: str, repo: str, branch: str, root: str) -> list[dict]:
    """Depth-first directory walk using the contents API.

    This avoids GitHub's huge recursive tree response for large repos.
    """

    stack = [root]
    files: list[dict] = []
    seen_dirs: set[str] = set()
    while stack:
        path = stack.pop()
        if path in seen_dirs:
            continue
        seen_dirs.add(path)
        try:
            payload = request_json(contents_url(owner, repo, branch, path))
        except RuntimeError as exc:
            print(f"[{repo}] warning: could not list {path or '<root>'}: {exc}")
            continue
        items = payload if isinstance(payload, list) else [payload]
        for item in items:
            item_type = item.get("type")
            item_path = item.get("path", "")
            if item_type == "dir":
                parts = set(Path(item_path).parts)
                if not (parts & SKIP_PARTS):
                    stack.append(item_path)
            elif item_type == "file":
                files.append(item)
    return files


def mirror_repo(owner: str, repo: str, branch: str, dirname: str, roots: list[str]) -> None:
    dest = DEST / dirname
    dest.mkdir(parents=True, exist_ok=True)
    by_path: dict[str, dict] = {}
    for root in roots:
        print(f"[{repo}] listing {root or '<root>'}", flush=True)
        for item in iter_contents(owner, repo, branch, root):
            by_path[item["path"]] = item
    entries = list(by_path.values())

    kept = 0
    skipped = 0
    print(f"[{repo}] candidate files: {len(entries)}", flush=True)
    for entry in entries:
        path = entry["path"]
        size = int(entry.get("size") or 0)
        if not should_keep(path, size):
            skipped += 1
            continue
        target = dest / path
        if target.exists() and target.stat().st_size == size:
            kept += 1
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        download_url = entry.get("download_url") or raw_url(owner, repo, branch, path)
        data = curl_bytes(download_url)
        target.write_bytes(data)
        kept += 1
        if kept % 50 == 0:
            print(f"[{repo}] downloaded {kept} files", flush=True)

    write_readme(dest, owner, repo, branch, kept, skipped)
    print(f"[{repo}] done: downloaded={kept}, skipped={skipped}, dest={dest}", flush=True)


def main() -> None:
    DEST.mkdir(parents=True, exist_ok=True)
    for owner, repo, branch, dirname, roots in REPOS:
        mirror_repo(owner, repo, branch, dirname, roots)


if __name__ == "__main__":
    main()
