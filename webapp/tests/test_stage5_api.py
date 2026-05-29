import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class Stage5ApiTests(unittest.TestCase):
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

    def create_project(self, name="stage5_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        project = self.projects_root / name
        server.write_json(
            project / "briefs" / "target_brief.json",
            {
                "target_catalog_id": "influenza_a_h1n1_na",
                "target_name": "Influenza A(H1N1) neuraminidase",
                "disease_context": "甲流",
                "binding_site": {"reference_ligand": "oseltamivir", "key_residues": ["Arg118", "Glu276"]},
            },
        )
        server.write_csv(
            project / "evidence" / "stage2_target_sources.csv",
            [
                {
                    "target_id": "influenza_a_h1n1_na",
                    "target": "Influenza A(H1N1) neuraminidase",
                    "evidence_score": "0.9967",
                    "readiness": "ready_for_mvp_closed_loop",
                    "primary_pdb": "3TI6",
                    "positive_controls": "oseltamivir; zanamivir; peramivir",
                    "pubmed_articles": "8",
                }
            ],
            ["target_id", "target", "evidence_score", "readiness", "primary_pdb", "positive_controls", "pubmed_articles"],
        )
        server.write_csv(
            project / "candidates" / "round_1_candidates.csv",
            [{"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"}],
            ["round", "id", "smiles", "parent", "source"],
        )
        server.write_csv(
            project / "scores" / "round_1_scores.csv",
            [{"id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "total_proxy": "0.72", "score_source": "proxy"}],
            ["id", "smiles", "total_proxy", "score_source"],
        )
        server.write_csv(
            project / "ranked" / "round_1_ranked.csv",
            [
                {
                    "rank": "1",
                    "id": "mol_001",
                    "smiles": "CCOC(=O)c1ccccc1",
                    "total_proxy": "0.72",
                    "docking_proxy": "0.5",
                    "pose_proxy": "0.8",
                    "decision": "advance",
                }
            ],
            ["rank", "id", "smiles", "total_proxy", "docking_proxy", "pose_proxy", "decision"],
        )
        server.write_json(
            project / "feedback" / "round_1_feedback.json",
            {"round": 1, "selected_count": 1, "selected_smiles": ["CCOC(=O)c1ccccc1"]},
        )
        server.write_csv(
            project / "seeds" / "round_1_seeds.csv",
            [{"id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "note": "advance"}],
            ["id", "smiles", "note"],
        )
        return name, project

    def test_stage5_dashboard_generates_assets_and_returns_data(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage5/dashboard",
            json={"round": 1, "title": "Stage 5 API Demo"},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 5)
        self.assertEqual(payload["command"], "stage5-dashboard")
        self.assertTrue(payload["has_dashboard"])
        self.assertEqual(payload["dashboard"]["stage"], 5)
        self.assertEqual(payload["dashboard"]["round"], 1)
        self.assertEqual(payload["dashboard"]["title"], "Stage 5 API Demo")
        self.assertEqual(payload["dashboard"]["metrics"]["raw_candidates"], 1)
        self.assertEqual(payload["dashboard"]["metrics"]["advanced"], 1)
        self.assertTrue((project / "stage5" / "dashboard_data.json").exists())
        self.assertTrue((project / "stage5" / "index.html").exists())
        self.assertTrue((project / "reports" / "stage5_dashboard_report.md").exists())

    def test_stage5_status_reads_existing_dashboard_and_report(self):
        name, _project = self.create_project()
        self.client.post(
            f"/api/projects/{name}/stage5/dashboard",
            json={"round": 1, "title": "Stage 5 API Demo"},
        )

        response = self.client.get(f"/api/projects/{name}/stage5?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 5)
        self.assertTrue(payload["has_dashboard"])
        self.assertTrue(payload["has_report"])
        self.assertEqual(payload["dashboard"]["target"]["id"], "influenza_a_h1n1_na")
        self.assertIn("Stage 5 Product Dashboard", payload["report"])

    def test_stage5_status_recomputes_requested_round_instead_of_serving_stale_dashboard(self):
        name, project = self.create_project()
        stage5_dir = project / "stage5"
        stage5_dir.mkdir(parents=True, exist_ok=True)
        server.write_json(
            stage5_dir / "dashboard_data.json",
            {
                "stage": 5,
                "round": 99,
                "readiness": {
                    "target_evidence": "ready_for_mvp_closed_loop",
                    "candidate_intake": "missing",
                    "rdkit_validation": "missing",
                    "control_panel": "missing",
                    "decoy_panel": "missing",
                    "docking": "missing",
                    "feedback": "missing",
                },
                "metrics": {},
                "target": {},
                "tables": {},
                "boundary": [],
            },
        )
        server.write_text(project / "reports" / "stage5_dashboard_report.md", "# stale report\n")

        response = self.client.get(f"/api/projects/{name}/stage5?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        dashboard = payload["dashboard"]
        self.assertEqual(dashboard["round"], 1)
        self.assertEqual(dashboard["readiness"]["target_evidence"], "ready_for_mvp_closed_loop")
        self.assertEqual(dashboard["readiness"]["candidate_intake"], "ready")
        self.assertEqual(dashboard["readiness"]["feedback"], "ready")


if __name__ == "__main__":
    unittest.main()
