import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class ExportPackageApiTests(unittest.TestCase):
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

    def create_project(self, name="export_bundle_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        project = self.projects_root / name
        server.write_json(
            project / "briefs" / "target_brief.json",
            {
                "target_catalog_id": "influenza_a_h1n1_na",
                "target_name": "Influenza A(H1N1) neuraminidase",
                "disease_context": "甲流",
                "binding_site": {"reference_ligand": "oseltamivir"},
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
                }
            ],
            ["target_id", "target", "evidence_score", "readiness", "primary_pdb", "reference_ligand", "positive_controls"],
        )
        server.write_csv(
            project / "candidates" / "round_1_candidates.csv",
            [{"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"}],
            ["round", "id", "smiles", "parent", "source"],
        )
        server.write_json(
            project / "stage7" / "round_1_delivery_manifest.json",
            {
                "delivery_status": "ready_for_demo_package",
                "deliverables": {
                    "stage7_manifest": {
                        "exists": True,
                        "path": str(project / "stage7" / "round_1_delivery_manifest.json"),
                        "relative_path": "stage7/round_1_delivery_manifest.json",
                        "description": "Manifest",
                    }
                },
            },
        )
        return name, project

    def test_export_package_bundle_lists_core_files(self):
        name, project = self.create_project()

        response = self.client.post(f"/api/projects/{name}/export-package", json={"round": 1})

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 8)
        self.assertIn(payload["package_status"], {"ready", "ready_with_warnings"})
        self.assertTrue(payload["files"]["zip"].endswith(".zip"))
        self.assertIn("target_pack_json", payload["included_files"])
        self.assertIn("stage7_manifest", payload["included_files"])
        self.assertTrue((project / "exports").exists())


if __name__ == "__main__":
    unittest.main()
