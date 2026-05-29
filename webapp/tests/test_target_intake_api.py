import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class TargetIntakeApiTests(unittest.TestCase):
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

    def create_project(self, name="target_intake_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        return name, self.projects_root / name

    def test_target_intake_normalizes_url_text_and_prompt_inputs(self):
        name, _project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/target-intake",
            json={
                "source_kind": "url",
                "source": "https://www.rcsb.org/structure/3TI6",
                "disease": "甲流",
                "target_hint": "neuraminidase",
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 1)
        self.assertEqual(payload["source_kind"], "url")
        self.assertTrue(payload["normalized_target"]["target_name"])
        self.assertTrue(payload["brief"]["target_name"])
        self.assertTrue(payload["prompt"])
        self.assertTrue(payload["files"]["target_brief"].endswith("target_brief.json"))
        self.assertTrue(payload["files"]["generator_prompt"].endswith("generator_prompt.md"))


if __name__ == "__main__":
    unittest.main()
