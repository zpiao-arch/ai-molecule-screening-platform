import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class ProjectRegistryApiTests(unittest.TestCase):
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

    def create_project(self, name="registry_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        return name, self.projects_root / name

    def test_project_list_reports_candidate_receptor_and_docking_status(self):
        name, project = self.create_project()
        server.write_csv(
            project / "candidates" / "round_2_candidates.csv",
            [{"round": "2", "id": "mol_001", "smiles": "CCO", "parent": "", "source": "unit"}],
            ["round", "id", "smiles", "parent", "source"],
        )
        server.write_json(
            project / "stage4" / "round_2_receptor_package.json",
            {
                "target_id": "influenza_a_h1n1_na",
                "pdb_id": "3TI6",
                "local_receptor_pdb": str(project / "stage4" / "receptors" / "3TI6.pdb"),
                "local_receptor_pdbqt": str(project / "stage4" / "receptors" / "3TI6.pdbqt"),
            },
        )
        (project / "stage4" / "receptors").mkdir(parents=True, exist_ok=True)
        (project / "stage4" / "receptors" / "3TI6.pdb").write_text("ATOM\n", encoding="utf-8")
        (project / "stage4" / "receptors" / "3TI6.pdbqt").write_text("ATOM\n", encoding="utf-8")
        server.write_json(
            project / "stage4" / "round_2_docking_plan.json",
            {"status": "completed", "selected_backend": "vina"},
        )
        server.write_csv(
            project / "stage4" / "round_2_docking_scores_template.csv",
            [{"id": "mol_001", "smiles": "CCO", "docking_score": "-7.2", "pose_pass": "true", "backend": "vina"}],
            ["id", "smiles", "docking_score", "pose_pass", "backend"],
        )

        response = self.client.get("/api/projects")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        project_row = next(row for row in payload if row["name"] == name)
        self.assertEqual(project_row["latest_round"], 2)
        status = project_row["asset_status"]
        self.assertEqual(status["candidate"]["status"], "ready")
        self.assertEqual(status["candidate"]["count"], 1)
        self.assertEqual(status["receptor"]["status"], "ready")
        self.assertEqual(status["receptor"]["pdb_id"], "3TI6")
        self.assertEqual(status["docking"]["status"], "completed")
        self.assertEqual(status["docking"]["score_count"], 1)
        self.assertEqual(status["docking"]["backend"], "vina")


if __name__ == "__main__":
    unittest.main()
