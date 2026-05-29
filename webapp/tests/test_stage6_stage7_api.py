import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class Stage6Stage7ApiTests(unittest.TestCase):
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

    def create_project(self, name="stage6_stage7_demo"):
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
            [
                {"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"},
                {"round": 1, "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "parent": "", "source": "unit"},
            ],
            ["round", "id", "smiles", "parent", "source"],
        )
        server.write_csv(
            project / "scores" / "round_1_scores.csv",
            [
                {"id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "total_proxy": "0.72", "score_source": "proxy"},
                {"id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "total_proxy": "0.61", "score_source": "proxy"},
            ],
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
                    "score_source": "proxy",
                    "decision": "advance",
                },
                {
                    "rank": "2",
                    "id": "mol_002",
                    "smiles": "CCS(=O)(=O)Nc1ccccc1",
                    "total_proxy": "0.61",
                    "docking_proxy": "0.4",
                    "pose_proxy": "0.7",
                    "score_source": "proxy",
                    "decision": "hold",
                },
            ],
            ["rank", "id", "smiles", "total_proxy", "docking_proxy", "pose_proxy", "score_source", "decision"],
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
        server.write_json(
            project / "stage4" / "round_1_stage4_assets.json",
            {"stage": 4, "round": 1, "target_id": "influenza_a_h1n1_na", "valid_count": 2, "invalid_count": 0},
        )
        server.write_json(
            project / "stage4" / "round_1_validation_metrics.json",
            {
                "readiness": {
                    "rdkit_descriptors": "ready",
                    "control_panel": "ready",
                    "decoy_panel": "ready",
                    "docking": "skipped",
                },
                "panel_counts": {"candidate": 2, "positive_control": 2, "reference_control": 1, "decoy": 2},
            },
        )
        server.write_json(
            project / "stage4" / "round_1_docking_plan.json",
            {"status": "skipped", "selected_backend": "none", "expected_scores_csv": str(project / "stage4" / "round_1_docking_scores_template.csv")},
        )
        server.write_csv(
            project / "stage4" / "round_1_real_descriptors.csv",
            [
                {"id": "mol_001", "canonical_smiles": "CCOC(=O)c1ccccc1", "valid": "1", "qed": "0.71", "mw": "150.18", "lipinski_violations": "0"},
                {"id": "mol_002", "canonical_smiles": "CCS(=O)(=O)Nc1ccccc1", "valid": "1", "qed": "0.62", "mw": "199.27", "lipinski_violations": "0"},
            ],
            ["id", "canonical_smiles", "valid", "qed", "mw", "lipinski_violations"],
        )
        server.write_csv(
            project / "stage4" / "round_1_benchmark_panel.csv",
            [
                {"panel_type": "candidate", "id": "mol_001", "canonical_smiles": "CCOC(=O)c1ccccc1", "qed": "0.71", "mw": "150.18", "nearest_control": "oseltamivir", "max_similarity_to_controls": "0.31"},
                {"panel_type": "candidate", "id": "mol_002", "canonical_smiles": "CCS(=O)(=O)Nc1ccccc1", "qed": "0.62", "mw": "199.27", "nearest_control": "zanamivir", "max_similarity_to_controls": "0.24"},
                {"panel_type": "positive_control", "id": "oseltamivir", "canonical_smiles": "CCC", "qed": "0.75", "mw": "312.4", "nearest_control": "", "max_similarity_to_controls": ""},
                {"panel_type": "reference_control", "id": "zanamivir", "canonical_smiles": "CCCC", "qed": "0.68", "mw": "332.3", "nearest_control": "", "max_similarity_to_controls": ""},
                {"panel_type": "decoy", "id": "decoy_001", "canonical_smiles": "c1ccccc1", "qed": "0.45", "mw": "78.1", "nearest_control": "", "max_similarity_to_controls": ""},
                {"panel_type": "decoy", "id": "decoy_002", "canonical_smiles": "CCN", "qed": "0.48", "mw": "45.1", "nearest_control": "", "max_similarity_to_controls": ""},
            ],
            ["panel_type", "id", "canonical_smiles", "qed", "mw", "nearest_control", "max_similarity_to_controls"],
        )
        server.write_text(project / "reports" / "stage3_round_1_report.md", "# Stage 3 report\n")
        server.write_text(project / "reports" / "stage4_round_1_report.md", "# Stage 4 report\n")
        return name, project

    def test_stage6_validate_generates_assets_and_returns_operations_payload(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage6/validate",
            json={"round": 1, "top": 5},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 6)
        self.assertEqual(payload["command"], "stage6-validate")
        self.assertTrue(payload["has_assets"])
        self.assertEqual(payload["assets"]["stage"], 6)
        self.assertIn("target_evidence", {row["gate_id"] for row in payload["quality_gates"]})
        self.assertEqual(payload["hit_triage"][0]["id"], "mol_001")
        self.assertTrue(payload["assay_queue"])
        self.assertTrue(payload["risk_register"])
        self.assertIn("Stage 6 Validation Operations", payload["report"])
        self.assertTrue((project / "stage6" / "round_1_validation_assets.json").exists())
        self.assertTrue((project / "stage6" / "round_1_quality_gates.csv").exists())
        self.assertTrue((project / "reports" / "stage6_round_1_validation_report.md").exists())

    def test_stage6_status_reads_existing_validation_assets(self):
        name, _project = self.create_project()
        self.client.post(f"/api/projects/{name}/stage6/validate", json={"round": 1, "top": 5})

        response = self.client.get(f"/api/projects/{name}/stage6?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertTrue(payload["has_assets"])
        self.assertTrue(payload["has_quality_gates"])
        self.assertEqual(payload["quality_gates"][0]["gate_id"], "target_evidence")
        self.assertEqual(payload["hit_triage"][0]["id"], "mol_001")
        self.assertIn("computational", " ".join(payload["boundary"]).lower())

    def test_stage7_package_generates_manifest_and_delivery_previews(self):
        name, project = self.create_project()
        self.client.post(f"/api/projects/{name}/stage6/validate", json={"round": 1, "top": 5})

        response = self.client.post(
            f"/api/projects/{name}/stage7/package",
            json={"round": 1, "title": "Flu NA delivery package"},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 7)
        self.assertEqual(payload["command"], "stage7-package")
        self.assertTrue(payload["has_manifest"])
        self.assertEqual(payload["manifest"]["stage"], 7)
        self.assertIn("executive_summary", payload["deliverables"])
        self.assertIn("stage6_quality_gates", payload["deliverables"])
        self.assertIn("stage8_frontend_spec", payload["deliverables"])
        self.assertIn("Stage 7 Delivery Package", payload["report"])
        self.assertIn("Flu NA delivery package", payload["executive_summary"])
        self.assertIn("Reproducibility Runbook", payload["reproducibility"])
        self.assertTrue(payload["checklist"])
        self.assertTrue((project / "stage7" / "round_1_delivery_manifest.json").exists())
        self.assertTrue((project / "stage7" / "stage8_frontend_product_spec.md").exists())
        self.assertTrue((project / "reports" / "stage7_round_1_delivery_report.md").exists())

    def test_stage7_status_reads_existing_delivery_package(self):
        name, _project = self.create_project()
        self.client.post(f"/api/projects/{name}/stage7/package", json={"round": 1, "title": "Status Package"})

        response = self.client.get(f"/api/projects/{name}/stage7?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertTrue(payload["has_manifest"])
        self.assertIn("executive_summary", payload["deliverables"])
        self.assertIn("Status Package", payload["executive_summary"])
        self.assertIn("Stage 8 Frontend Product Specification", payload["stage8_spec"])
        self.assertIn("computational", " ".join(payload["boundary"]).lower())


if __name__ == "__main__":
    unittest.main()
