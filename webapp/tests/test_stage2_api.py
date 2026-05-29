import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class Stage2ApiTests(unittest.TestCase):
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

    def create_project(self, name="stage2_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        return name, self.projects_root / name

    def test_stage2_evidence_writes_matrix_assets_and_report(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage2/evidence",
            json={"disease": "甲流", "target": "influenza_a_h1n1_na", "top": 1, "offline": True},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 2)
        self.assertEqual(payload["command"], "evidence-stage2")
        self.assertEqual(payload["disease"], "甲流")
        self.assertEqual(len(payload["matrix"]), 1)
        self.assertEqual(payload["matrix"][0]["target_id"], "influenza_a_h1n1_na")
        self.assertIn("targets", payload["assets"])
        self.assertIn("Stage 2 Target Evidence Source Report", payload["report"])
        self.assertTrue((project / "evidence" / "stage2_target_sources.csv").exists())
        self.assertTrue((project / "evidence" / "stage2_closed_loop_assets.json").exists())
        self.assertTrue((project / "reports" / "stage2_evidence_report.md").exists())

    def test_stage2_status_reports_existing_matrix_assets_and_report(self):
        name, _project = self.create_project()
        self.client.post(
            f"/api/projects/{name}/stage2/evidence",
            json={"disease": "甲流", "target": "influenza_a_h1n1_na", "top": 1, "offline": True},
        )

        response = self.client.get(f"/api/projects/{name}/stage2")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertTrue(payload["has_matrix"])
        self.assertTrue(payload["has_assets"])
        self.assertTrue(payload["has_report"])
        self.assertEqual(payload["matrix"][0]["target_id"], "influenza_a_h1n1_na")
        self.assertEqual(payload["assets"]["disease"], "甲流")


if __name__ == "__main__":
    unittest.main()
