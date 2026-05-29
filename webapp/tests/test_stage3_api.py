import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class Stage3ApiTests(unittest.TestCase):
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

    def create_project(self, name="stage3_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        return name, self.projects_root / name

    def test_stage3_status_reports_existing_candidate_assets(self):
        name, project = self.create_project()
        response = self.client.post(
            f"/api/projects/{name}/stage3/candidates",
            json={
                "round": 1,
                "source_mode": "text",
                "source_text": "CCO\nCCN\nc1ccccc1O",
                "n": 6,
                "top": 2,
                "use_openai": False,
            },
        )
        self.assertEqual(response.status_code, 200, response.text)

        response = self.client.get(f"/api/projects/{name}/stage3?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 3)
        self.assertEqual(payload["project"], name)
        self.assertEqual(payload["round"], 1)
        self.assertTrue(payload["has_assets"])
        self.assertGreaterEqual(len(payload["raw_candidates"]), 1)
        self.assertGreaterEqual(len(payload["candidates"]), 1)
        self.assertIn("stage3_assets", payload["files"])
        self.assertIn("raw_candidates", payload["files"])
        self.assertFalse(payload["api_key_security"]["api_key_persisted"])
        self.assertIn("computational", payload["boundary"].lower())
        self.assertTrue((project / "stage3" / "round_1_stage3_assets.json").exists())


if __name__ == "__main__":
    unittest.main()
