import json
import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class Stage1ApiTests(unittest.TestCase):
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

    def create_project(self, name="stage1_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        return name, self.projects_root / name

    def test_stage1_target_select_writes_ranking_and_returns_rows(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage1/target-select",
            json={"disease": "甲流", "top": 2},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 1)
        self.assertEqual(payload["command"], "target-select")
        self.assertEqual(payload["disease"], "甲流")
        self.assertEqual(len(payload["targets"]), 2)
        self.assertEqual(payload["targets"][0]["rank"], "1")
        self.assertTrue(payload["targets"][0]["target_id"])
        self.assertIn("files", payload)
        self.assertTrue((project / "targets" / "target_selection.csv").exists())
        self.assertTrue((project / "reports" / "target_selection.md").exists())

    def test_stage1_brief_from_target_writes_brief_prompt_and_updates_config(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage1/brief-from-target",
            json={
                "disease": "甲流",
                "target": "influenza_a_h1n1_na",
                "free_text": "围绕甲流神经氨酸酶口袋生成虚拟筛选候选分子。",
                "force": True,
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 1)
        self.assertEqual(payload["command"], "brief-from-target")
        self.assertEqual(payload["brief"]["target_catalog_id"], "influenza_a_h1n1_na")
        self.assertIn("Generator Prompt", payload["prompt"])
        self.assertTrue((project / "briefs" / "target_brief.json").exists())
        self.assertTrue((project / "prompts" / "generator_prompt.md").exists())
        config = json.loads((project / "config.json").read_text(encoding="utf-8"))
        self.assertEqual(config["target"]["target_catalog_id"], "influenza_a_h1n1_na")

    def test_stage1_status_reports_existing_selection_brief_and_prompt(self):
        name, _project = self.create_project()
        self.client.post(
            f"/api/projects/{name}/stage1/target-select",
            json={"disease": "甲流", "top": 1},
        )
        self.client.post(
            f"/api/projects/{name}/stage1/brief-from-target",
            json={"disease": "甲流", "target": "influenza_a_h1n1_na", "force": True},
        )

        response = self.client.get(f"/api/projects/{name}/stage1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertTrue(payload["has_target_selection"])
        self.assertTrue(payload["has_brief"])
        self.assertTrue(payload["has_prompt"])
        self.assertEqual(len(payload["target_selection"]), 1)
        self.assertEqual(payload["brief"]["target_catalog_id"], "influenza_a_h1n1_na")


if __name__ == "__main__":
    unittest.main()
