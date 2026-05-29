#!/usr/bin/env python3
"""Run a browser-driven Stage 1-8 closed-loop acceptance test.

This script intentionally uses the web UI for the operational steps. It keeps
API calls only for independent evidence capture after the UI has created the
artifacts.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from playwright.sync_api import sync_playwright


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_ROOT = ROOT / "outputs" / "frontend_e2e_20260526"


def dump_json(out_dir: Path, name: str, data: object) -> None:
    (out_dir / name).write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


def api(page, path: str, method: str = "GET", body: object | None = None):
    script = """async ({path, method, body}) => {
      const opts = { method, headers: {'Content-Type': 'application/json'} };
      if (body !== null && body !== undefined) opts.body = JSON.stringify(body);
      const res = await fetch(path, opts);
      const text = await res.text();
      let data = null;
      try { data = text ? JSON.parse(text) : null; } catch (e) { data = {raw: text}; }
      return { ok: res.ok, status: res.status, data };
    }"""
    result = page.evaluate(script, {"path": path, "method": method, "body": body})
    if not result["ok"]:
        raise AssertionError(f"API {method} {path} failed {result['status']}: {result['data']}")
    return result["data"]


def wait_button_done(page, selector: str, timeout: int = 180_000) -> None:
    page.wait_for_function(
        "sel => { const b = document.querySelector(sel); return b && !b.disabled; }",
        arg=selector,
        timeout=timeout,
    )


def select_project(page, selector: str, project: str) -> None:
    page.wait_for_function(
        """({selector, project}) => {
          const s = document.querySelector(selector);
          return s && Array.from(s.options).some(o => o.value === project);
        }""",
        arg={"selector": selector, "project": project},
        timeout=30_000,
    )
    page.select_option(selector, project)


def assert_tokens(page, selector: str, tokens: list[str]) -> str:
    locator = page.locator(selector).first
    text = locator.inner_text(timeout=10_000) if page.locator(selector).count() else ""
    if not any(token in text for token in tokens):
        text = page.locator("body").inner_text(timeout=10_000)
    if not any(token in text for token in tokens):
        raise AssertionError(f"{selector} missing tokens {tokens}. Text: {text[:500]}")
    return text


def run_frontend_closed_loop(base_url: str, project: str, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    round_no = "1"
    errors: list[dict] = []

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        context = browser.new_context(viewport={"width": 1440, "height": 1200}, ignore_https_errors=True)
        page = context.new_page()
        page.on("console", lambda msg: errors.append({"type": msg.type, "text": msg.text}) if msg.type == "error" else None)
        page.on("pageerror", lambda exc: errors.append({"type": "pageerror", "text": str(exc)}))

        page.goto(base_url.rstrip("/") + f"/?e2e={project}", wait_until="domcontentloaded", timeout=30_000)
        page.wait_for_selector('#sidebar-nav .nav-item[data-page="projects"]', timeout=30_000)
        page.wait_for_timeout(1_000)
        page.screenshot(path=str(out_dir / "00_home.png"), full_page=True)

        page.click('[data-page="projects"]')
        page.wait_for_selector("#page-projects.active", timeout=10_000)
        page.click('button:has-text("新建项目")')
        page.fill("#new-proj-name", project)
        page.click("#new-project-modal button.btn-primary")
        page.wait_for_function(
            "name => Array.from(document.querySelectorAll('#projects-list *')).some(el => el.textContent.includes(name))",
            arg=project,
            timeout=30_000,
        )
        page.screenshot(path=str(out_dir / "01_project_created.png"), full_page=True)

        page.click('[data-page="brief"]')
        page.wait_for_selector("#page-brief.active", timeout=10_000)
        select_project(page, "#bf-project-select", project)
        page.fill("#st1-disease", "甲流")
        page.fill("#st1-source", "https://www.rcsb.org/structure/3TI6")
        page.fill("#st1-target-hint", "neuraminidase oseltamivir pocket")
        page.click("#st1-select-btn")
        wait_button_done(page, "#st1-select-btn", timeout=60_000)
        assert_tokens(page, "#stage1-target-results", ["influenza_a_h1n1_na"])
        page.click("#st1-intake-btn")
        wait_button_done(page, "#st1-intake-btn", timeout=60_000)
        page.click("#st1-target-pack-btn")
        wait_button_done(page, "#st1-target-pack-btn", timeout=60_000)
        page.click("#st1-target-pack-validate-btn")
        wait_button_done(page, "#st1-target-pack-validate-btn", timeout=60_000)
        page.screenshot(path=str(out_dir / "02_stage1_complete.png"), full_page=True)
        stage1 = api(page, f"/api/projects/{project}/stage1")
        dump_json(out_dir, "stage1.json", stage1)
        if not all(stage1.get(k) for k in ["has_target_selection", "has_brief", "has_prompt", "has_target_pack"]):
            raise AssertionError(f"Stage 1 incomplete: {stage1}")

        page.click('[data-page="stage2"]')
        page.wait_for_selector("#page-stage2.active", timeout=10_000)
        select_project(page, "#st2-project-select", project)
        page.fill("#st2-disease", "甲流")
        page.fill("#st2-target", "influenza_a_h1n1_na")
        page.click("#st2-run-btn")
        wait_button_done(page, "#st2-run-btn", timeout=90_000)
        assert_tokens(page, "#stage2-evidence-results", ["influenza_a_h1n1_na"])
        page.screenshot(path=str(out_dir / "03_stage2_complete.png"), full_page=True)
        stage2 = api(page, f"/api/projects/{project}/stage2")
        dump_json(out_dir, "stage2.json", stage2)
        if not stage2.get("has_matrix"):
            raise AssertionError(f"Stage 2 matrix missing: {stage2}")

        page.click('[data-page="stage8"]')
        page.wait_for_selector("#page-stage8.active", timeout=10_000)
        select_project(page, "#st8-project-select", project)
        page.fill("#st8-round", round_no)
        page.select_option("#st8-candidate-source-mode", "proxy")
        page.fill("#st8-candidate-n", "12")
        page.fill("#st8-candidate-top", "4")
        page.click("#st8-candidate-btn")
        wait_button_done(page, "#st8-candidate-btn", timeout=120_000)
        assert_tokens(page, "#stage8-candidate-entry", ["通过过滤"])
        page.screenshot(path=str(out_dir / "04_stage3_candidate_entry.png"), full_page=True)

        page.click('[data-page="stage4"]')
        page.wait_for_selector("#page-stage4.active", timeout=10_000)
        select_project(page, "#st4-project-select", project)
        page.fill("#st4-round", round_no)
        page.fill("#st4-target", "influenza_a_h1n1_na")
        page.fill("#st4-top", "4")
        page.fill("#st4-decoys", "4")
        page.fill("#st4-max-conformers", "1")
        page.fill("#st4-pdb-id", "3TI6")
        page.fill("#st4-pocket-center", "-28.914,14.334,20.794")
        page.fill("#st4-pocket-size", "23.585,20.45,24.18")
        page.fill("#st4-pocket-source", "co_crystal_ligand_frontend_e2e")
        page.check("#st4-fetch-receptor")
        page.check("#st4-run-docking")
        page.select_option("#st4-docking-backend", "vina")
        page.fill("#st4-docking-timeout", "600")
        page.click("#st4-run-btn")
        wait_button_done(page, "#st4-run-btn", timeout=300_000)
        assert_tokens(page, "#stage4-status", ["Stage 4"])
        page.wait_for_function(
            "() => { const t=document.querySelector('#stage4-structure-status')?.textContent || ''; return t.includes('已加载') || t.includes('结构投影'); }",
            timeout=30_000,
        )
        page.screenshot(path=str(out_dir / "05_stage4_real_docking.png"), full_page=True)
        stage4 = api(page, f"/api/projects/{project}/stage4?round=1")
        dump_json(out_dir, "stage4.json", stage4)
        docking_plan = stage4.get("docking_plan") or {}
        if docking_plan.get("status") != "completed" or not stage4.get("docking_results"):
            raise AssertionError(f"Stage 4 docking incomplete: {docking_plan}")

        page.click("#st45-run-btn")
        wait_button_done(page, "#st45-run-btn", timeout=300_000)
        page.click("#st46-run-btn")
        wait_button_done(page, "#st46-run-btn", timeout=120_000)
        page.screenshot(path=str(out_dir / "06_stage45_stage46.png"), full_page=True)
        stage45 = api(page, f"/api/projects/{project}/stage45?round=1")
        stage46 = api(page, f"/api/projects/{project}/stage46?round=1")
        readiness = api(page, f"/api/projects/{project}/scientific-readiness?round=1")
        dump_json(out_dir, "stage45.json", stage45)
        dump_json(out_dir, "stage46.json", stage46)
        dump_json(out_dir, "scientific_readiness.json", readiness)
        docking = (stage45.get("validation") or {}).get("docking") or {}
        if docking.get("status") != "completed" or int(docking.get("scored_count") or 0) <= 0:
            raise AssertionError(f"Stage 4.5 docking incomplete: {docking}")
        if not stage46.get("has_benchmark"):
            raise AssertionError(f"Stage 4.6 benchmark missing: {stage46}")

        page.click('[data-page="stage8"]')
        page.wait_for_selector("#page-stage8.active", timeout=10_000)
        select_project(page, "#st8-project-select", project)
        page.fill("#st8-round", round_no)
        for selector in [
            "#st8-demo-package-btn",
            "#st8-preflight-btn",
            "#st8-report-btn",
            "#st8-evidence-pack-btn",
            "#st8-full-export-btn",
            "#st8-demo-doctor-btn",
        ]:
            page.click(selector)
            wait_button_done(page, selector, timeout=120_000)
        page.screenshot(path=str(out_dir / "07_stage8_exports_doctor.png"), full_page=True)
        command = api(page, f"/api/projects/{project}/stage8/command-center?round=1")
        preflight = api(page, f"/api/projects/{project}/stage8/preflight?round=1")
        doctor = api(page, f"/api/projects/{project}/demo-doctor?round=1")
        dump_json(out_dir, "stage8_command_center.json", command)
        dump_json(out_dir, "stage8_preflight.json", preflight)
        dump_json(out_dir, "demo_doctor.json", doctor)

        rail = {str(row.get("stage")): row.get("status") for row in command.get("stage_rail", [])}
        if any(rail.get(str(i)) != "ready" for i in range(1, 8)):
            raise AssertionError(f"Stage rail not ready: {rail}")
        if preflight.get("status") != "ready":
            raise AssertionError(f"Stage 8 preflight not ready: {preflight.get('status')}")
        if doctor.get("overall_status") not in {"ready_for_demo", "ready_with_warnings"}:
            raise AssertionError(f"Demo Doctor not ready: {doctor.get('overall_status')}")

        visual_probe = {
            "status_text": page.locator("#stage4-structure-status").inner_text() if page.locator("#stage4-structure-status").count() else "",
            "pose_options": 0,
            "canvas_count": 0,
            "svg_count": 0,
        }
        page.click('[data-page="stage4"]')
        select_project(page, "#st4-project-select", project)
        page.wait_for_function(
            "() => { const t=document.querySelector('#stage4-structure-status')?.textContent || ''; return t.includes('已加载') || t.includes('结构投影'); }",
            timeout=30_000,
        )
        page.wait_for_timeout(1_000)
        visual_probe = {
            "status_text": page.locator("#stage4-structure-status").inner_text(),
            "meta_text": page.locator("#stage4-viewer-meta").inner_text(),
            "pose_options": page.locator("#stage4-pose-select option").count(),
            "canvas_count": page.locator("#stage4-structure-viewer canvas").count(),
            "svg_count": page.locator("#stage4-structure-viewer svg").count(),
        }
        dump_json(out_dir, "stage4_visual_probe.json", visual_probe)
        page.screenshot(path=str(out_dir / "stage4_visual_probe.png"), full_page=True)
        if visual_probe["pose_options"] <= 0 or (visual_probe["canvas_count"] <= 0 and visual_probe["svg_count"] <= 0):
            raise AssertionError(f"Stage 4 visual probe failed: {visual_probe}")

        final = {
            "project": project,
            "stage1_flags": {k: stage1.get(k) for k in ["has_target_selection", "has_brief", "has_prompt", "has_target_pack"]},
            "stage4_docking_status": docking_plan.get("status"),
            "stage4_docking_results": len(stage4.get("docking_results") or []),
            "stage45_docking": docking,
            "stage46_auc": ((stage46.get("benchmark") or {}).get("metrics") or {}).get("roc_auc"),
            "scientific_readiness": readiness.get("readiness_level"),
            "claim_level": readiness.get("claim_level"),
            "stage8_preflight": preflight.get("status"),
            "stage8_rail": rail,
            "demo_doctor": doctor.get("overall_status"),
            "stage4_visual_probe": visual_probe,
            "console_errors": errors,
            "screenshots_dir": str(out_dir),
        }
        dump_json(out_dir, "final_frontend_e2e_summary.json", final)
        browser.close()

    if errors:
        raise AssertionError(f"Browser console errors: {errors}")
    return final


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the browser-driven AI molecule closed-loop acceptance test.")
    parser.add_argument("--base-url", default="http://localhost:8765/", help="Running FastAPI web URL.")
    parser.add_argument("--project", default=f"frontend_full_loop_{time.strftime('%Y%m%d_%H%M%S')}", help="Fresh project name to create.")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT, help="Directory for screenshots and JSON evidence.")
    args = parser.parse_args()

    out_dir = args.output_root / args.project
    final = run_frontend_closed_loop(args.base_url, args.project, out_dir)
    print(json.dumps(final, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
