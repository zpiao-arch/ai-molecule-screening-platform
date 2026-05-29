import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server

import test_stage6_stage7_api as stage67_fixture


class ProductNextStepApiTests(unittest.TestCase):
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

    def create_project(self, name="next_step_demo"):
        return stage67_fixture.Stage6Stage7ApiTests.create_project(self, name)

    def test_product_bootstrap_health_exposes_startup_and_demo_checks(self):
        name, _project = self.create_project("flu_na_real_demo")

        response = self.client.get(f"/api/product/bootstrap-health?project={name}&round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "product-bootstrap-health")
        self.assertEqual(payload["project"], name)
        self.assertIn(payload["startup_status"], {"ready", "ready_with_warnings", "needs_setup"})
        self.assertIn("startup", payload["sections"])
        self.assertIn("demo", payload["sections"])
        self.assertIn("delivery", payload["sections"])
        self.assertTrue(payload["quick_start"])
        self.assertTrue(any("./start_web.sh" in item for item in payload["quick_start"]))
        self.assertIn("http://localhost:8765/", payload["urls"]["app"])
        self.assertTrue(any("Computational" in item for item in payload["boundary"]))

    def test_stage4_operator_guide_summarizes_docking_qc_and_templates(self):
        name, _project = self.create_project("stage4_operator_demo")

        response = self.client.get(f"/api/projects/{name}/stage4/operator-guide?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage4-operator-guide")
        self.assertEqual(payload["project"], name)
        self.assertIn("docking_status", payload)
        self.assertIn("score_csv_template", payload)
        self.assertTrue(payload["workflow_steps"])
        self.assertTrue(any("external" in step["action"].lower() for step in payload["workflow_steps"]))
        self.assertIn("posebusters", {item["tool"] for item in payload["tool_readiness"]})
        self.assertTrue(any("Docking" in item for item in payload["boundary"]))

    def test_stage8_action_plan_surfaces_primary_action_and_autofix_api(self):
        response = self.client.post("/api/projects", json={"name": "action_plan_demo"})
        self.assertEqual(response.status_code, 200, response.text)

        response = self.client.get("/api/projects/action_plan_demo/stage8/action-plan?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage8-action-plan")
        self.assertEqual(payload["project"], "action_plan_demo")
        self.assertIn(payload["overall_status"], {"ready", "ready_with_warnings", "needs_repair"})
        self.assertIn("primary_action", payload)
        self.assertIn("autofix", payload["primary_action"])
        self.assertTrue(payload["missing_or_warn"])
        self.assertTrue(any(item["api"].endswith("/stage8/repair") for item in payload["actions"]))

    def test_target_catalog_includes_cross_disease_expansion_targets(self):
        response = self.client.get("/api/target-catalog")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        ids = {row["id"] for row in payload["targets"]}
        self.assertIn("egfr_tyrosine_kinase", ids)
        self.assertIn("bace1_amyloid_beta_secretase", ids)
        self.assertIn("hiv1_protease", ids)
        egfr = next(row for row in payload["targets"] if row["id"] == "egfr_tyrosine_kinase")
        self.assertTrue(egfr["primary_pdb"])
        self.assertTrue(egfr["known_drugs"])

    def test_generator_adapters_report_available_external_generation_options(self):
        response = self.client.get("/api/generator-adapters")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "generator-adapters")
        adapter_ids = {row["id"] for row in payload["adapters"]}
        for adapter_id in ["proxy_smiles", "openai_prompt", "reinvent4", "drugex", "diffsbdd", "colabfold"]:
            self.assertIn(adapter_id, adapter_ids)
        self.assertTrue(any(row["status"] in {"available", "local_repo", "planned"} for row in payload["adapters"]))
        self.assertIn("api_key_policy", payload)
        self.assertIn("not persisted", payload["api_key_policy"].lower())

    def test_product_deep_integration_collects_six_next_value_tracks(self):
        name, _project = self.create_project("deep_integration_demo")

        response = self.client.get(f"/api/projects/{name}/product/deep-integration?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "product-deep-integration")
        self.assertEqual(payload["project"], name)
        self.assertEqual(payload["round"], 1)
        self.assertIn(payload["overall_status"], {"ready", "ready_with_warnings", "needs_work"})
        self.assertIn("summary", payload)
        self.assertIn("canonical_urls", payload)
        self.assertEqual(len(payload["focus_areas"]), 6)
        area_ids = {row["area_id"] for row in payload["focus_areas"]}
        for area_id in [
            "real_demo_hardening",
            "generator_adapters",
            "stage4_scientific_trust",
            "stage8_command_center",
            "delivery_reproducibility",
            "scientific_report_layer",
        ]:
            self.assertIn(area_id, area_ids)
        for row in payload["focus_areas"]:
            self.assertIn(row["status"], {"ready", "warn", "missing"})
            self.assertTrue(row["label"])
            self.assertTrue(row["evidence"])
            self.assertTrue(row["route"])
            self.assertTrue(row["api"])
            self.assertTrue(row["actions"])
        self.assertIn("/api/generator-adapters", payload["canonical_urls"]["generator_adapters"])
        self.assertTrue(any("computational" in item.lower() for item in payload["boundary"]))

    def test_standard_delivery_build_script_creates_manifest_and_zip(self):
        name, _project = self.create_project("delivery_script_demo")

        result = subprocess.run(
            [
                sys.executable,
                "scripts/build_product_delivery.py",
                "--project",
                name,
                "--round",
                "1",
                "--projects-root",
                str(self.projects_root),
                "--skip-tests",
            ],
            cwd=server.ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        self.assertIn("delivery_zip:", result.stdout)
        project = self.projects_root / name
        manifest = project / "exports" / "round_1_standard_delivery_manifest.json"
        self.assertTrue(manifest.exists())
        payload = server.read_json(manifest)
        self.assertEqual(payload["command"], "standard-delivery-build")
        self.assertIn("full_export", payload["included_files"])
        self.assertTrue(Path(payload["files"]["zip"]).exists())


if __name__ == "__main__":
    unittest.main()
