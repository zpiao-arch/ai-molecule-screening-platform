import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

from fastapi.testclient import TestClient

import server

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_stage6_stage7_api as stage67_fixture  # noqa: E402


class Stage8ApiTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=Path.cwd())
        self.projects_root = Path(self.tmp.name) / "projects"
        self.projects_root.mkdir(parents=True, exist_ok=True)
        self.original_projects_root = server.PROJECTS_ROOT
        server.PROJECTS_ROOT = self.projects_root
        self.client = TestClient(server.app)

    def tearDown(self):
        server.PROJECTS_ROOT = self.original_projects_root
        self.tmp.cleanup()

    def create_project(self, name="stage8_demo"):
        return stage67_fixture.Stage6Stage7ApiTests.create_project(self, name)

    def test_stage8_command_center_aggregates_stage1_to_stage7(self):
        name, _project = self.create_project()
        self.client.post(f"/api/projects/{name}/stage6/validate", json={"round": 1, "top": 5})
        self.client.post(f"/api/projects/{name}/stage7/package", json={"round": 1, "title": "Stage 8 package"})

        response = self.client.get(f"/api/projects/{name}/stage8/command-center?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 8)
        self.assertEqual(payload["project"], name)
        self.assertEqual(payload["round"], 1)
        self.assertTrue(payload["has_stage5"])
        self.assertTrue(payload["has_stage6"])
        self.assertTrue(payload["has_stage7"])
        self.assertEqual(payload["target"]["id"], "influenza_a_h1n1_na")
        for key in ["target_evidence", "candidate_intake", "rdkit_validation", "docking", "feedback"]:
            self.assertIn(key, payload["readiness"])
        for key in ["raw_candidates", "ranked_candidates", "advanced"]:
            self.assertIn(key, payload["metrics"])
        self.assertEqual([row["stage"] for row in payload["stage_rail"]], [1, 2, 3, 4, 5, 6, 7])
        self.assertGreaterEqual(payload["quality_gate_summary"]["total"], 1)
        self.assertEqual(payload["hit_triage"][0]["id"], "mol_001")
        deliverable_keys = {row["key"] for row in payload["deliverables"]}
        self.assertIn("executive_summary", deliverable_keys)
        self.assertIn("stage6_quality_gates", deliverable_keys)
        self.assertTrue(payload["next_actions"])
        boundary = " ".join(payload["boundary"]).lower()
        self.assertIn("computational", boundary)
        self.assertIn("efficacy", boundary)

    def test_stage8_command_center_surfaces_missing_downstream_assets(self):
        name, _project = self.create_project("stage8_missing_demo")

        response = self.client.get(f"/api/projects/{name}/stage8/command-center?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 8)
        self.assertFalse(payload["has_stage6"])
        self.assertFalse(payload["has_stage7"])
        statuses = {row["stage"]: row["status"] for row in payload["stage_rail"]}
        self.assertEqual(statuses[6], "missing")
        self.assertEqual(statuses[7], "missing")
        self.assertTrue(any("Stage 6" in action for action in payload["next_actions"]))
        self.assertTrue(payload["files"]["stage5_dashboard_data"].endswith("stage5/dashboard_data.json"))

    def test_stage8_demo_runner_builds_missing_closed_loop_artifacts(self):
        response = self.client.post("/api/projects", json={"name": "runner_demo"})
        self.assertEqual(response.status_code, 200, response.text)

        response = self.client.post(
            "/api/projects/runner_demo/stage8/demo-runner",
            json={
                "round": 1,
                "disease": "甲流",
                "target": "influenza_a_h1n1_na",
                "source_kind": "url",
                "source": "https://www.rcsb.org/structure/3TI6",
                "target_hint": "neuraminidase oseltamivir pocket",
                "candidates": 18,
                "top": 5,
                "run_docking": False,
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage8-demo-runner")
        self.assertEqual(payload["project"], "runner_demo")
        self.assertEqual(payload["round"], 1)
        self.assertEqual(payload["status"], "completed")
        steps = {step["step_id"]: step for step in payload["steps"]}
        for step_id in [
            "target_selection",
            "target_intake",
            "target_pack",
            "stage2_evidence",
            "candidate_generation",
            "proxy_scoring",
            "ranking",
            "feedback",
            "stage4_assets",
            "stage45_controls",
            "stage46_benchmark",
            "stage5_dashboard",
            "stage6_validation",
            "stage7_delivery",
            "stage8_command_center",
        ]:
            self.assertIn(step_id, steps)
            self.assertIn(steps[step_id]["status"], {"completed", "ready", "completed_with_warnings"})
        self.assertTrue(payload["artifacts"]["target_pack_json"].endswith("target_pack.json"))
        self.assertTrue(payload["artifacts"]["ranked_candidates"].endswith("round_1_ranked.csv"))
        self.assertGreater(payload["summary"]["candidate_count"], 0)
        self.assertTrue(payload["stage8"]["has_stage6"])
        self.assertTrue(payload["stage8"]["has_stage7"])
        target_selection_path = self.projects_root / "runner_demo" / "targets" / "target_selection.csv"
        target_report_path = self.projects_root / "runner_demo" / "reports" / "target_selection.md"
        self.assertTrue(target_selection_path.exists())
        self.assertTrue(target_report_path.exists())
        target_rows = server.csv_to_dicts(target_selection_path)
        self.assertGreaterEqual(len(target_rows), 1)
        self.assertEqual(target_rows[0]["target_id"], "influenza_a_h1n1_na")
        self.assertTrue((self.projects_root / "runner_demo" / "target_pack.json").exists())
        self.assertTrue((self.projects_root / "runner_demo" / "ranked" / "round_1_ranked.csv").exists())
        self.assertTrue((self.projects_root / "runner_demo" / "stage7" / "round_1_delivery_manifest.json").exists())

    def test_stage8_demo_runner_real_docking_args_fetch_receptor_and_forward_pocket_box(self):
        req = server.Stage8DemoRunnerRequest(
            round=1,
            disease="甲流",
            target="influenza_a_h1n1_na",
            candidates=8,
            top=3,
            decoys=2,
            run_docking=True,
            docking_backend="vina",
        )
        target_pack = {
            "pdb_id": "3TI6",
            "pocket": {
                "center": [-28.914, 14.334, 20.794],
                "size": [23.585, 20.45, 24.18],
                "source": "co_crystal_ligand",
            },
        }

        args = server.stage8_stage4_real_args(
            self.projects_root / "arg_demo",
            1,
            req,
            "influenza_a_h1n1_na",
            target_pack,
            {"pdb_id": "3TI6"},
            "vina",
        )

        self.assertTrue(args.fetch_receptor)
        self.assertTrue(args.run_docking)
        self.assertEqual(args.pdb_id, "3TI6")
        self.assertEqual(args.pocket_center, "-28.914,14.334,20.794")
        self.assertEqual(args.pocket_size, "23.585,20.45,24.18")
        self.assertEqual(args.pocket_source, "co_crystal_ligand")

    def test_stage8_report_preflight_repair_catalog_and_candidate_entrypoints(self):
        response = self.client.post("/api/projects", json={"name": "ops_demo"})
        self.assertEqual(response.status_code, 200, response.text)

        preflight = self.client.get("/api/projects/ops_demo/stage8/preflight?round=1")
        self.assertEqual(preflight.status_code, 200, preflight.text)
        preflight_payload = preflight.json()
        self.assertEqual(preflight_payload["stage"], 8)
        self.assertTrue(any(item["status"] == "missing" for item in preflight_payload["checks"]))
        self.assertTrue(preflight_payload["can_auto_repair"])

        repair = self.client.post(
            "/api/projects/ops_demo/stage8/repair",
            json={"round": 1, "target": "influenza_a_h1n1_na", "candidates": 12, "top": 4},
        )
        self.assertEqual(repair.status_code, 200, repair.text)
        self.assertEqual(repair.json()["command"], "stage8-repair")
        self.assertTrue(repair.json()["runner"]["summary"]["candidate_count"])

        report = self.client.post("/api/projects/ops_demo/stage8/report", json={"round": 1})
        self.assertEqual(report.status_code, 200, report.text)
        report_payload = report.json()
        self.assertEqual(report_payload["command"], "stage8-report")
        self.assertIn("target_rationale", report_payload["sections"])
        self.assertIn("computational", " ".join(report_payload["boundary"]).lower())
        self.assertTrue((self.projects_root / "ops_demo" / "reports" / "stage8_round_1_closed_loop_report.md").exists())

        catalog = self.client.get("/api/target-catalog?query=甲流")
        self.assertEqual(catalog.status_code, 200, catalog.text)
        self.assertGreaterEqual(catalog.json()["count"], 1)
        self.assertTrue(any(row["id"] == "influenza_a_h1n1_na" for row in catalog.json()["targets"]))

        custom = self.client.post(
            "/api/projects/ops_demo/target-catalog/custom",
            json={
                "target_id": "demo_custom_target",
                "display_name": "Demo custom target",
                "disease": "demo disease",
                "pdb_id": "1ABC",
                "reference_ligand": "demo ligand",
                "known_drugs": "control_a; control_b",
            },
        )
        self.assertEqual(custom.status_code, 200, custom.text)
        self.assertEqual(custom.json()["target"]["id"], "demo_custom_target")
        self.assertTrue((self.projects_root / "ops_demo" / "targets" / "custom_target_catalog.json").exists())

        candidate = self.client.post(
            "/api/projects/ops_demo/stage3/candidates",
            json={
                "round": 2,
                "source_mode": "text",
                "source_text": "CCO\nCCN\nc1ccccc1O",
                "n": 6,
                "top": 2,
                "use_openai": False,
            },
        )
        self.assertEqual(candidate.status_code, 200, candidate.text)
        candidate_payload = candidate.json()
        self.assertEqual(candidate_payload["command"], "stage3-candidates")
        self.assertEqual(candidate_payload["assets"]["raw_candidates"], 3)
        self.assertGreaterEqual(candidate_payload["assets"]["passed_filter"], 1)
        self.assertTrue((self.projects_root / "ops_demo" / "stage3" / "round_2_stage3_assets.json").exists())

    def test_stage8_acceptance_demo_builds_clean_project_and_report(self):
        response = self.client.post(
            "/api/stage8/acceptance-demo",
            json={
                "project_name": "flu_acceptance_demo",
                "force": True,
                "round": 1,
                "candidates": 12,
                "top": 4,
                "run_docking": False,
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage8-acceptance-demo")
        self.assertEqual(payload["project"], "flu_acceptance_demo")
        self.assertEqual(payload["status"], "passed")
        self.assertEqual(payload["target"]["id"], "influenza_a_h1n1_na")
        self.assertEqual(payload["preflight"]["status"], "ready")
        checks = {row["check_id"]: row for row in payload["checks"]}
        for check_id in [
            "target_catalog",
            "target_selection",
            "target_pack",
            "candidate_intake",
            "stage4_assets",
            "stage5_dashboard",
            "stage6_validation",
            "stage7_delivery",
            "stage8_preflight",
            "closed_loop_report",
            "download_materials",
        ]:
            self.assertEqual(checks[check_id]["status"], "passed", check_id)
        target_selection_path = self.projects_root / "flu_acceptance_demo" / "targets" / "target_selection.csv"
        self.assertTrue(target_selection_path.exists())
        self.assertGreaterEqual(len(server.csv_to_dicts(target_selection_path)), 1)
        self.assertTrue((self.projects_root / "flu_acceptance_demo" / "acceptance" / "stage8_acceptance_report.json").exists())
        self.assertTrue((self.projects_root / "flu_acceptance_demo" / "acceptance" / "stage8_acceptance_report.md").exists())
        self.assertTrue((self.projects_root / "flu_acceptance_demo" / "reports" / "stage8_round_1_closed_loop_report.md").exists())

    def test_stage8_demo_guide_exposes_ordered_wizard_steps_and_next_action(self):
        name, _project = self.create_project("guide_demo")

        response = self.client.get(f"/api/projects/{name}/stage8/demo-guide?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage8-demo-guide")
        self.assertEqual(payload["product_mode"], "guided_demo")
        step_ids = [row["step_id"] for row in payload["steps"]]
        for step_id in [
            "project_setup",
            "target_pack",
            "target_evidence",
            "candidate_intake",
            "stage4_validation",
            "stage5_dashboard",
            "stage6_validation",
            "stage7_delivery",
            "stage8_report",
            "evidence_pack",
        ]:
            self.assertIn(step_id, step_ids)
        self.assertEqual(payload["next_primary_action"]["step_id"], "target_pack")
        self.assertEqual(payload["next_primary_action"]["status"], "missing")
        self.assertIn("api", payload["next_primary_action"])
        self.assertGreater(payload["completion"]["total"], payload["completion"]["ready"])
        self.assertTrue(any("computational" in item.lower() for item in payload["boundary"]))

    def test_stage8_review_mode_summarizes_demo_story_for_reviewers(self):
        name, _project = self.create_project("review_demo")
        self.client.post(f"/api/projects/{name}/stage6/validate", json={"round": 1, "top": 5})
        self.client.post(f"/api/projects/{name}/stage7/package", json={"round": 1, "title": "Review package"})
        self.client.post(f"/api/projects/{name}/stage8/evidence-pack", json={"round": 1, "title": "Review evidence pack"})

        response = self.client.get(f"/api/projects/{name}/stage8/review-mode?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 8)
        self.assertEqual(payload["command"], "stage8-review-mode")
        self.assertEqual(payload["product_mode"], "review_mode")
        self.assertEqual(payload["project"], name)
        self.assertEqual(payload["round"], 1)
        self.assertIn("computational screening", payload["positioning"].lower())
        story_ids = [row["section_id"] for row in payload["storyline"]]
        self.assertEqual(
            story_ids,
            ["project_goal", "target_evidence", "candidate_screening", "computational_validation", "delivery_evidence"],
        )
        self.assertEqual(payload["target"]["id"], "influenza_a_h1n1_na")
        self.assertGreaterEqual(payload["evidence_strength"]["ready_steps"], 5)
        self.assertIn(payload["review_status"], {"ready_for_review", "ready_with_warnings", "needs_repair"})
        self.assertIn("stage8_report", payload["primary_actions"])
        self.assertIn("evidence_pack_zip", payload["download_links"])
        self.assertTrue(payload["talk_track"])
        self.assertTrue(any("no efficacy" in item.lower() for item in payload["claim_boundary"]))

    def test_stage8_evidence_pack_exports_zip_manifest_report_and_download_links(self):
        name, project = self.create_project("evidence_pack_demo")
        self.client.post(f"/api/projects/{name}/stage6/validate", json={"round": 1, "top": 5})
        self.client.post(f"/api/projects/{name}/stage7/package", json={"round": 1, "title": "Evidence package"})

        response = self.client.post(
            f"/api/projects/{name}/stage8/evidence-pack",
            json={"round": 1, "title": "Evidence Pack"},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage8-evidence-pack")
        self.assertEqual(payload["project"], name)
        self.assertIn(payload["package_status"], {"ready", "ready_with_warnings"})
        self.assertTrue(payload["files"]["zip"].endswith("round_1_evidence_pack.zip"))
        self.assertTrue(payload["files"]["manifest"].endswith("round_1_evidence_manifest.json"))
        for key in ["target_pack", "closed_loop_report", "stage8_command_center", "ranked_candidates", "stage7_manifest", "evidence_manifest"]:
            self.assertIn(key, payload["included_files"])
            self.assertNotIn(key, payload["missing_files"])
        self.assertIn("evidence_pack_zip", payload["download_links"])
        self.assertIn(f"/api/projects/{name}/artifact?path=", payload["download_links"]["evidence_pack_zip"])

        download = self.client.get(payload["download_links"]["evidence_pack_zip"])
        self.assertEqual(download.status_code, 200, download.text)
        zip_path = project / "exports" / "round_1_evidence_pack.zip"
        self.assertTrue(zip_path.exists())
        with zipfile.ZipFile(zip_path) as zf:
            names = set(zf.namelist())
        self.assertIn("stage8_command_center.json", names)
        self.assertIn("reports/stage8_round_1_closed_loop_report.md", names)
        self.assertIn("exports/round_1_evidence_manifest.json", names)

    def test_stage3_candidate_entrypoint_reports_api_key_safety_without_leaking_secret(self):
        response = self.client.post("/api/projects", json={"name": "key_safety_demo"})
        self.assertEqual(response.status_code, 200, response.text)

        secret = "test-api-key-secret-value"
        response = self.client.post(
            "/api/projects/key_safety_demo/stage3/candidates",
            json={
                "round": 1,
                "source_mode": "text",
                "source_text": "CCO\nCCN",
                "n": 4,
                "top": 2,
                "use_openai": False,
                "api_key": secret,
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        self.assertNotIn(secret, response.text)
        payload = response.json()
        self.assertFalse(payload["api_key_security"]["api_key_persisted"])
        self.assertTrue(payload["api_key_security"]["api_key_input_received"])
        self.assertIn("not written", payload["api_key_security"]["policy"])

        project_dir = self.projects_root / "key_safety_demo"
        for path in project_dir.rglob("*"):
            if path.is_file():
                try:
                    content = path.read_text(encoding="utf-8")
                except UnicodeDecodeError:
                    continue
                self.assertNotIn(secret, content, str(path))


if __name__ == "__main__":
    unittest.main()
