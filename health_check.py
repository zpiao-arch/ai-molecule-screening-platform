#!/usr/bin/env python3
"""Summarize the Product Ops system health endpoint for local delivery checks."""

from __future__ import annotations

import json
import os
import sys
import urllib.error
import urllib.request


def main() -> int:
    port = os.environ.get("PORT", "8765")
    url = os.environ.get("SYSTEM_HEALTH_URL", f"http://localhost:{port}/api/system/health")
    try:
        with urllib.request.urlopen(url, timeout=10) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except urllib.error.URLError as exc:
        print(f"system health unreachable: {url}")
        print(f"error: {exc}")
        print("try: python3 webapp/server.py")
        return 1

    overall_status = payload.get("overall_status", "unknown")
    print(f"overall_status: {overall_status}")
    print(f"url: {payload.get('url', '-')}")
    print(f"default_project: {payload.get('default_project', '-')}")

    sections = payload.get("sections", {}) or {}
    for key, section in sections.items():
        label = section.get("label", key)
        status = section.get("status", "-")
        print(f"\n[{label}] {status}")
        for check in section.get("checks", []) or []:
            check_status = check.get("status", "-")
            evidence = check.get("evidence", "-")
            remedy = check.get("remedy", "")
            line = f"- {check.get('label', check.get('check_id', '-'))}: {check_status} ({evidence})"
            if remedy:
                line += f" -> {remedy}"
            print(line)

    print("\nrecommended_commands:")
    for command in payload.get("recommended_commands", []) or []:
        print(f"- {command}")

    if payload.get("boundary"):
        print("\nboundary:")
        for item in payload.get("boundary", []) or []:
            print(f"- {item}")

    return 0 if overall_status in {"ready", "ready_with_warnings"} else 2


if __name__ == "__main__":
    sys.exit(main())
