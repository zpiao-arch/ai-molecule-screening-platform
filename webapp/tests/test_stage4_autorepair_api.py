import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class Stage4AutoRepairApiTests(unittest.TestCase):
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

    def create_project(self, name="stage4_autorepair_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        project = self.projects_root / name
        server.write_json(
            project / "briefs" / "target_brief.json",
            {
                "target_catalog_id": "influenza_a_h1n1_na",
                "target_name": "Influenza A(H1N1) neuraminidase",
                "disease_context": "甲流",
                "binding_site": {
                    "reference_ligand": "oseltamivir",
                    "key_residues": ["Arg118", "Glu276"],
                },
            },
        )
        return name, project

    def test_stage4_repair_generates_missing_candidates_for_requested_round(self):
        name, project = self.create_project()
        response = self.client.get(f"/api/projects/{name}/stage4?round=2")
        self.assertEqual(response.status_code, 200, response.text)
        before = response.json()
        candidate_gate = next(item for item in before["stage4_preflight"] if item["step_id"] == "candidate_input")
        self.assertEqual(candidate_gate["status"], "missing")

        response = self.client.post(f"/api/projects/{name}/stage4/repair", json={"round": 2})

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage4-repair")
        self.assertTrue(payload["repaired"])
        repairs = {item["repair_id"]: item for item in payload["repairs"]}
        self.assertEqual(repairs["stage3_candidates"]["status"], "completed")
        self.assertTrue((project / "candidates" / "round_2_candidates.csv").exists())
        self.assertTrue((project / "stage3" / "round_2_stage3_assets.json").exists())
        after_gate = next(item for item in payload["stage4_preflight"] if item["step_id"] == "candidate_input")
        self.assertEqual(after_gate["status"], "ready")
        joined_next_actions = "\n".join(payload.get("next_actions", []))
        self.assertNotIn("可先补候选", joined_next_actions)
        self.assertIn("运行 Stage 4", joined_next_actions)

    def test_stage4_doctor_marks_candidate_gap_as_auto_repairable(self):
        name, _project = self.create_project()

        response = self.client.get(f"/api/projects/{name}/doctor?round=2")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        issue = next(item for item in payload["issues"] if item["step_id"] == "candidate_input")
        self.assertTrue(issue["auto_repair"])
        self.assertEqual(issue["repair_endpoint"], f"/api/projects/{name}/stage4/repair")
        self.assertIn("Stage 4修复", issue["next_action"])

    def test_stage8_repair_forwards_run_docking_policy(self):
        name, _project = self.create_project("stage8_repair_policy_demo")

        response = self.client.post(
            f"/api/projects/{name}/stage8/repair",
            json={"round": 1, "run_docking": True, "candidates": 6, "top": 3, "decoys": 2},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertIn("runner", payload)
        stage4_step = next(step for step in payload["runner"]["steps"] if step["step_id"] == "stage4_assets")
        self.assertIn("docking_plan", stage4_step["files"])
        plan = server.read_json(Path(stage4_step["files"]["docking_plan"]))
        self.assertTrue(plan.get("run_docking_requested"))


if __name__ == "__main__":
    unittest.main()
