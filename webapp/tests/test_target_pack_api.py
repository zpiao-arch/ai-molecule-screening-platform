import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class TargetPackApiTests(unittest.TestCase):
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

    def create_project(self, name="target_pack_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        project = self.projects_root / name
        server.write_json(
            project / "briefs" / "target_brief.json",
            {
                "schema_version": "0.1",
                "target_catalog_id": "influenza_a_h1n1_na",
                "target_name": "Influenza A(H1N1) neuraminidase",
                "disease_context": "甲流",
                "target_rationale": "Validated NA target with public co-crystal structures.",
                "free_text_requirement": "围绕3TI6共晶口袋生成候选分子。",
                "protein": {"name": "Influenza A(H1N1) neuraminidase", "pdb_id": "3TI6"},
                "binding_site": {
                    "description": "Neuraminidase active site",
                    "reference_ligand": "oseltamivir",
                    "key_residues": ["Arg118", "Glu276"],
                    "center": [],
                    "size": [],
                    "source": "co_crystal",
                    "box_strategy": "Use co-crystal ligand centroid.",
                },
                "design_intent": {"style": "analog_design"},
                "generation_constraints": {"max_heavy_atoms": 55, "max_molecular_weight": 550.0},
                "validation_plan": {"known_positive_controls": ["oseltamivir", "zanamivir"]},
                "compliance_boundary": {"scope": "computational screening"},
                "source_urls": ["https://www.rcsb.org/structure/3TI6"],
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
                    "reference_ligand": "oseltamivir",
                    "positive_controls": "oseltamivir; zanamivir; peramivir",
                    "reference_controls": "laninamivir",
                    "assay_path": "Neuraminidase inhibition assay",
                    "binding_site_source": "co_crystal",
                    "box_strategy": "Use co-crystallized oseltamivir centroid from 3TI6.",
                }
            ],
            [
                "target_id",
                "target",
                "evidence_score",
                "readiness",
                "primary_pdb",
                "reference_ligand",
                "positive_controls",
                "reference_controls",
                "assay_path",
                "binding_site_source",
                "box_strategy",
            ],
        )
        server.write_json(
            project / "config.json",
            {
                "target": {
                    "name": "Influenza A(H1N1) neuraminidase",
                    "target_catalog_id": "influenza_a_h1n1_na",
                    "pdb_id": "3TI6",
                    "reference_ligand_sdf": "oseltamivir",
                    "pocket": {"center": [0.0, 0.0, 0.0], "size": [20.0, 20.0, 20.0], "source": "co_crystal"},
                }
            },
        )
        server.write_csv(
            project / "candidates" / "round_1_candidates.csv",
            [{"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"}],
            ["round", "id", "smiles", "parent", "source"],
        )
        return name, project

    def test_target_pack_endpoint_builds_pack_and_lists_artifacts(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/target-pack",
            json={"round": 1, "target": "influenza_a_h1n1_na", "pdb_id": "3TI6", "reference_ligand": "oseltamivir"},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 1)
        self.assertIn(payload["status"], {"ready", "ready_with_warnings"})
        self.assertEqual(payload["target_pack"]["target_id"], "influenza_a_h1n1_na")
        self.assertTrue(payload["target_pack"]["pocket"]["center"])
        self.assertTrue(payload["target_pack"]["controls"]["positive"])
        self.assertTrue(payload["files"]["target_pack_json"].endswith("target_pack.json"))
        self.assertTrue((project / "target_pack.json").exists())
        self.assertTrue((project / "reports" / "target_pack_report.md").exists())

        status = self.client.get(f"/api/projects/{name}/target-pack?round=1")
        self.assertEqual(status.status_code, 200, status.text)
        self.assertTrue(status.json()["has_target_pack"])


if __name__ == "__main__":
    unittest.main()
