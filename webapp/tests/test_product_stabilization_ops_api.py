import tempfile
import time
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server

import test_stage6_stage7_api as stage67_fixture


class ProductStabilizationOpsApiTests(unittest.TestCase):
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

    def create_project(self, name="stabilization_demo"):
        return stage67_fixture.Stage6Stage7ApiTests.create_project(self, name)

    def add_scientific_controls(self, project: Path):
        server.write_csv(
            project / "stage4_5" / "round_1_control_docking_scores.csv",
            [
                {
                    "panel_type": "positive_control",
                    "id": "oseltamivir",
                    "docking_score": "-8.1",
                    "pose_pass": "true",
                    "backend": "vina",
                },
                {
                    "panel_type": "candidate",
                    "id": "mol_001",
                    "docking_score": "-7.2",
                    "pose_pass": "true",
                    "backend": "vina",
                },
                {
                    "panel_type": "decoy",
                    "id": "decoy_001",
                    "docking_score": "-4.3",
                    "pose_pass": "false",
                    "backend": "vina",
                },
            ],
            ["panel_type", "id", "docking_score", "pose_pass", "backend"],
        )
        server.write_json(
            project / "stage4_5" / "round_1_control_validation.json",
            {
                "counts": {"candidate": 1, "known_control": 1, "decoy": 1},
                "docking": {"status": "completed", "backend": "vina", "pose_pass_count": 2},
                "candidate_vs_controls": {
                    "best_candidate_delta_kcal_mol": 0.9,
                    "candidate_beats_or_matches_best_control": False,
                },
                "redocking": {"rmsd_status": "not_available"},
            },
        )
        server.write_json(
            project / "stage4_6" / "round_1_retrospective_benchmark.json",
            {
                "metrics": {
                    "roc_auc": 1.0,
                    "top_k": {
                        "top_1": {
                            "k": 1,
                            "control_hit_count": 1,
                            "control_hit_rate": 1.0,
                            "enrichment_factor": 3.0,
                        }
                    },
                    "best_candidate": {"id": "mol_001", "rank": 2},
                },
                "counts": {"positives": 1, "negatives": 1, "candidates": 1},
            },
        )

    def test_system_health_exposes_product_readiness_checks(self):
        self.create_project("flu_na_real_demo")

        response = self.client.get("/api/system/health")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertIn(payload["overall_status"], {"ready", "ready_with_warnings", "needs_setup"})
        self.assertIn("runtime", payload["sections"])
        self.assertIn("frontend", payload["sections"])
        self.assertIn("computational_tools", payload["sections"])
        self.assertIn("delivery", payload["sections"])
        self.assertTrue(payload["recommended_commands"])
        flattened = [item for section in payload["sections"].values() for item in section["checks"]]
        self.assertTrue(flattened)
        for item in flattened:
            self.assertIn("status", item)
            self.assertIn("evidence", item)
            self.assertIn("remedy", item)
        self.assertIn("Computational", " ".join(payload["boundary"]))

    def test_system_health_uses_current_stage4_capability_schema(self):
        self.create_project("flu_na_real_demo")
        original_capabilities = server.stage4_capabilities
        server.stage4_capabilities = lambda: {
            "modules": {
                "rdkit": {"name": "rdkit", "status": "available", "version": "2022.09.5"},
                "openbabel": {"name": "openbabel", "status": "available", "version": "3.1.0"},
                "posebusters": {"name": "posebusters", "status": "available", "version": "0.6.5"},
                "meeko": {"name": "meeko", "status": "available", "version": "0.7.1"},
                "gemmi": {"name": "gemmi", "status": "available", "version": "0.7.5"},
            },
            "executables": {
                "vina": {"name": "vina", "status": "found", "path": "/tmp/vina"},
                "gnina": {"name": "gnina", "status": "not_found", "path": ""},
                "obabel": {"name": "obabel", "status": "found", "path": "/tmp/obabel"},
                "mk_prepare_ligand.py": {"name": "mk_prepare_ligand.py", "status": "found", "path": "/tmp/mk_prepare_ligand.py"},
                "mk_prepare_receptor.py": {"name": "mk_prepare_receptor.py", "status": "found", "path": "/tmp/mk_prepare_receptor.py"},
                "bust": {"name": "bust", "status": "found", "path": "/tmp/bust"},
            },
            "docking_backend": {
                "status": "available",
                "available_backends": ["vina"],
                "message": "Vina is available.",
            },
        }
        try:
            response = self.client.get("/api/system/health")
        finally:
            server.stage4_capabilities = original_capabilities

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["sections"]["computational_tools"]["status"], "ready")
        checks = {
            item["check_id"]: item
            for item in payload["sections"]["computational_tools"]["checks"]
        }
        self.assertEqual(checks["rdkit"]["status"], "ready")
        self.assertEqual(checks["openbabel"]["status"], "ready")
        self.assertEqual(checks["posebusters"]["status"], "ready")
        self.assertEqual(checks["meeko"]["status"], "ready")
        self.assertEqual(checks["vina"]["status"], "ready")
        self.assertEqual(checks["docking_backend"]["status"], "ready")
        self.assertEqual(checks["gnina"]["status"], "ready")
        self.assertIn("optional", checks["gnina"]["evidence"].lower())
        self.assertEqual(payload["overall_status"], "ready")

    def test_target_pack_validate_reports_schema_and_missing_coordinates(self):
        name, _project = self.create_project()
        self.client.post(f"/api/projects/{name}/target-pack", json={"round": 1})

        response = self.client.get(f"/api/projects/{name}/target-pack/validate?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 1)
        self.assertEqual(payload["command"], "target-pack-validate")
        self.assertIn(payload["status"], {"ready", "ready_with_warnings", "invalid"})
        self.assertIn("target_id", payload["checks"])
        self.assertIn("pdb_id", payload["checks"])
        self.assertIn("controls", payload["checks"])
        self.assertIn("pocket_box", payload["checks"])
        self.assertIn("evidence_matrix", payload["checks"])
        self.assertTrue(payload["next_actions"])
        self.assertTrue(Path(payload["files"]["validation_json"]).exists())

    def test_scientific_readiness_summarizes_controls_decoys_and_claim_level(self):
        name, project = self.create_project()
        self.add_scientific_controls(project)

        response = self.client.get(f"/api/projects/{name}/scientific-readiness?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "scientific-readiness")
        self.assertIn(payload["readiness_level"], {"strong_computational", "moderate_computational", "planning_only"})
        self.assertEqual(payload["claim_level"], "computational_screening_only")
        self.assertGreaterEqual(payload["metrics"]["stage45_scored"], 3)
        self.assertTrue(payload["controls"]["has_positive_control"])
        self.assertTrue(payload["decoys"]["has_decoy"])
        self.assertIn("control_decoy_separation", payload["checks"])
        self.assertTrue(Path(payload["files"]["readiness_json"]).exists())

    def test_job_queue_runs_full_export_and_exposes_status(self):
        name, _project = self.create_project()

        response = self.client.post(
            "/api/jobs",
            json={"task": "stage8_full_export", "project": name, "round": 1, "payload": {"title": "Queued Export"}},
        )

        self.assertEqual(response.status_code, 200, response.text)
        job = response.json()
        self.assertEqual(job["task"], "stage8_full_export")
        self.assertIn(job["status"], {"queued", "running", "completed"})

        job_id = job["job_id"]
        final = job
        for _ in range(30):
            status_response = self.client.get(f"/api/jobs/{job_id}")
            self.assertEqual(status_response.status_code, 200, status_response.text)
            final = status_response.json()
            if final["status"] in {"completed", "failed"}:
                break
            time.sleep(0.05)

        self.assertEqual(final["status"], "completed", final)
        self.assertEqual(final["result"]["command"], "stage8-full-export")
        self.assertTrue(Path(final["result"]["files"]["zip"]).exists())

        list_response = self.client.get("/api/jobs")
        self.assertEqual(list_response.status_code, 200, list_response.text)
        self.assertTrue(any(row["job_id"] == job_id for row in list_response.json()["jobs"]))

    def test_stage8_full_export_adds_health_validation_and_scientific_readiness(self):
        name, project = self.create_project()
        self.add_scientific_controls(project)

        response = self.client.post(
            f"/api/projects/{name}/stage8/full-export",
            json={"round": 1, "title": "Full Product Export"},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage8-full-export")
        self.assertIn(payload["package_status"], {"ready", "ready_with_warnings"})
        self.assertIn("system_health", payload["included_files"])
        self.assertIn("target_pack_validation", payload["included_files"])
        self.assertIn("scientific_readiness", payload["included_files"])
        self.assertTrue(Path(payload["files"]["zip"]).exists())
        self.assertTrue(Path(payload["files"]["manifest"]).exists())
        self.assertIn("download_links", payload)
        self.assertIn("full_export_zip", payload["download_links"])
        self.assertIn("full_export_manifest", payload["download_links"])
        manifest = server.read_json(Path(payload["files"]["manifest"]))
        self.assertIn("scientific_readiness", manifest["included_files"])
        self.assertIn("Computational", " ".join(manifest["boundary"]))
        self.assertTrue((project / "exports" / "round_1_scientific_readiness.json").exists())


if __name__ == "__main__":
    unittest.main()
