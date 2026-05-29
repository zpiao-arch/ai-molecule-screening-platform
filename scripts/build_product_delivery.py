#!/usr/bin/env python3
"""Build a standard product delivery zip for one project round."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import zipfile
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "webapp"))

import server  # noqa: E402


def run_tests() -> dict:
    cmd = [sys.executable, "-m", "pytest", "webapp/tests", "-q"]
    proc = subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=False)
    return {
        "command": " ".join(cmd),
        "returncode": proc.returncode,
        "stdout_tail": proc.stdout[-4000:],
        "stderr_tail": proc.stderr[-4000:],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build standard product delivery package.")
    parser.add_argument("--project", required=True)
    parser.add_argument("--round", type=int, default=1)
    parser.add_argument("--projects-root", default="")
    parser.add_argument("--title", default="")
    parser.add_argument("--skip-tests", action="store_true")
    args = parser.parse_args()

    if args.projects_root:
        server.PROJECTS_ROOT = Path(args.projects_root)
    project_dir = server.PROJECTS_ROOT / args.project
    if not project_dir.exists():
        print(f"project not found: {project_dir}", file=sys.stderr)
        return 2

    safe_round = max(1, int(args.round))
    title = args.title or f"{args.project} Standard Product Delivery"
    test_result = {"skipped": True} if args.skip_tests else run_tests()
    if not args.skip_tests and test_result.get("returncode") != 0:
        print(json.dumps(test_result, ensure_ascii=False, indent=2), file=sys.stderr)
        return 3

    full_export = server.stage8_full_export_bundle(args.project, project_dir, safe_round, title)
    bootstrap = server.product_bootstrap_health_payload(args.project, safe_round)
    action_plan = server.stage8_action_plan_payload(args.project, project_dir, safe_round)
    generator_adapters = server.generator_adapters_payload()

    export_dir = project_dir / "exports"
    export_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = export_dir / f"round_{safe_round}_standard_delivery_manifest.json"
    zip_path = export_dir / f"round_{safe_round}_standard_delivery.zip"

    manifest = {
        "schema_version": "0.1",
        "command": "standard-delivery-build",
        "project": args.project,
        "round": safe_round,
        "title": title,
        "generated_at": datetime.now().isoformat(),
        "test_result": test_result,
        "included_files": ["full_export", "bootstrap_health", "stage8_action_plan", "generator_adapters", "standard_manifest"],
        "files": {
            "zip": str(zip_path),
            "manifest": str(manifest_path),
            "full_export_zip": str(full_export.get("files", {}).get("zip", "")),
        },
        "boundary": [
            "Standard delivery package for computational screening demo and handoff.",
            "No potency, efficacy, toxicity, safety, dosing, clinical benefit, or therapeutic claim is created.",
        ],
    }
    server.write_json(manifest_path, manifest)

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.write(manifest_path, arcname=manifest_path.relative_to(project_dir).as_posix())
        for key, obj in [
            ("bootstrap_health.json", bootstrap),
            ("stage8_action_plan.json", action_plan),
            ("generator_adapters.json", generator_adapters),
            ("full_export_summary.json", full_export),
        ]:
            zf.writestr(key, json.dumps(obj, ensure_ascii=False, indent=2))
        full_zip = Path(str(full_export.get("files", {}).get("zip", "")))
        if full_zip.exists():
            zf.write(full_zip, arcname=f"exports/{full_zip.name}")
        zf.writestr(
            "README.md",
            "\n".join(
                [
                    f"# {title}",
                    "",
                    f"- Project: `{args.project}`",
                    f"- Round: `{safe_round}`",
                    "- This is a computational screening product delivery package.",
                    "- Open `bootstrap_health.json` first, then `stage8_action_plan.json`.",
                    "- No efficacy, safety, toxicity, dosing, clinical benefit, or therapeutic claim is made.",
                ]
            )
            + "\n",
        )

    print(f"delivery_manifest: {manifest_path}")
    print(f"delivery_zip: {zip_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
