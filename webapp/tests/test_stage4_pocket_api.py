import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class Stage4PocketApiTests(unittest.TestCase):
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

    def create_project(self, name="stage4_pocket_demo"):
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
                    "assay_path": "Neuraminidase inhibition assay",
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
                "assay_path",
            ],
        )
        server.write_csv(
            project / "candidates" / "round_1_candidates.csv",
            [{"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"}],
            ["round", "id", "smiles", "parent", "source"],
        )
        receptor_dir = project / "stage4" / "receptors"
        receptor_dir.mkdir(parents=True, exist_ok=True)
        (receptor_dir / "3TI6.pdb").write_text(
            "\n".join(
                [
                    "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N",
                    "HETATM    2  C1  G39 A 801       8.000   9.000  10.000  1.00 10.00           C",
                    "HETATM    3  C2  G39 A 801       9.000  10.000  11.000  1.00 10.00           C",
                    "HETATM    4  O1  G39 A 801       7.000   9.500  10.500  1.00 10.00           O",
                    "END",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return name, project

    def test_stage4_pocket_pack_uses_co_crystal_reference(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage4/pocket-pack",
            json={"round": 1, "pdb_id": "3TI6", "reference_ligand": "oseltamivir"},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 4)
        self.assertIn(payload["pocket"]["status"], {"ready", "needs_review"})
        self.assertTrue(payload["pocket"]["center"])
        self.assertIn(payload["pocket"]["source"], {"co_crystal", "curated", "manual"})
        self.assertTrue(payload["files"]["pocket_pack"].endswith("pocket_pack.json"))
        self.assertTrue((project / "stage4" / "pocket_pack.json").exists())


if __name__ == "__main__":
    unittest.main()
