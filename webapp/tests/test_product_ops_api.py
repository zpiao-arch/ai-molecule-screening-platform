import sys
import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_stage6_stage7_api as stage67_fixture  # noqa: E402


class ProductOpsApiTests(unittest.TestCase):
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

    def create_project(self, name="product_ops_demo"):
        return stage67_fixture.Stage6Stage7ApiTests.create_project(self, name)

    def test_health_reports_runtime_and_first_available_default_project(self):
        self.create_project("flu_na_real_demo")
        self.create_project("aaa_current_project")

        response = self.client.get("/api/health")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["default_project"], "aaa_current_project")
        self.assertIn("http://localhost:8765", payload["url"])
        self.assertTrue(payload["checks"]["index_html"])
        self.assertTrue(payload["checks"]["projects_root"])

    def test_static_vendor_assets_are_served_locally(self):
        response = self.client.get("/static/vendor/3Dmol-min.js")

        self.assertEqual(response.status_code, 200, response.text[:200])
        self.assertIn("javascript", response.headers.get("content-type", ""))
        self.assertGreater(len(response.content), 100_000)

    def test_index_response_disables_browser_cache_for_frontend_updates(self):
        response = self.client.get("/")

        self.assertEqual(response.status_code, 200, response.text[:200])
        self.assertIn("no-store", response.headers.get("cache-control", ""))
        self.assertIn("结构投影兜底视图", response.text)
        self.assertNotIn("缺少 receptor 或 pose URL", response.text)

    def test_stage8_demo_package_generates_demo_assets_downloads_and_command_center(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage8/demo-package",
            json={"round": 1, "top": 5, "title": "One click demo"},
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 8)
        self.assertEqual(payload["command"], "stage8-demo-package")
        self.assertTrue(payload["generated"]["stage5"])
        self.assertTrue(payload["generated"]["stage6"])
        self.assertTrue(payload["generated"]["stage7"])
        self.assertTrue(payload["has_stage5"])
        self.assertTrue(payload["has_stage6"])
        self.assertTrue(payload["has_stage7"])
        self.assertTrue((project / "stage5" / "index.html").exists())
        self.assertTrue((project / "stage7" / "round_1_delivery_manifest.json").exists())
        self.assertIn("stage7_manifest", payload["download_links"])
        self.assertIn("ranked_candidates", payload["download_links"])

        download = self.client.get(payload["download_links"]["stage7_manifest"])
        self.assertEqual(download.status_code, 200, download.text)
        self.assertIn("ready_for_demo_package", download.text)

    def test_demo_doctor_aggregates_stage8_static_assets_and_system_checks(self):
        name, _project = self.create_project()
        self.client.post(
            f"/api/projects/{name}/stage8/demo-package",
            json={"round": 1, "top": 5, "title": "One click demo"},
        )

        response = self.client.get(f"/api/projects/{name}/demo-doctor?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 8)
        self.assertEqual(payload["project"], name)
        self.assertIn(payload["overall_status"], {"ready_for_demo", "ready_with_warnings"})
        self.assertIn("checks", payload)
        self.assertIn("stage_pipeline", payload["checks"])
        self.assertIn("static_3dmol", payload["checks"])
        self.assertEqual(payload["checks"]["static_3dmol"]["status"], "ready")
        self.assertIn("system_tools", payload["checks"])
        self.assertIn("stage8", payload["checks"])
        self.assertTrue(payload["stage8"]["has_stage6"])
        self.assertTrue(payload["stage8"]["has_stage7"])
        self.assertTrue(payload["next_actions"])

    def test_artifact_download_rejects_paths_outside_project(self):
        name, _project = self.create_project("artifact_guard_demo")

        response = self.client.get(f"/api/projects/{name}/artifact?path=../../server.py")

        self.assertEqual(response.status_code, 400, response.text)

    def test_doctor_reports_local_eval_tool_bin_executables(self):
        vina = server.DEFAULT_EVAL_TOOLS / "bin" / "vina"
        self.assertTrue(vina.exists(), "fixture expects the local eval tool bin vina already installed")

        response = self.client.get("/api/doctor")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertIn("vina", payload["executables"])
        self.assertEqual(payload["executables"]["vina"]["status"], "found")
        vina_path = Path(payload["executables"]["vina"]["path"])
        self.assertTrue(vina_path.exists())
        self.assertTrue(vina_path.is_file())
        self.assertIn(str(server.DEFAULT_EVAL_TOOLS), str(vina_path))
        self.assertFalse(payload["proxy_only"])
        self.assertTrue(payload["capabilities"]["docking_backend"]["available_backends"])


class ProductOpsScriptContractTests(unittest.TestCase):
    def test_start_script_and_web_readme_document_health_check(self):
        root = Path(__file__).resolve().parents[2]
        script = root / "start_web.sh"
        readme = root / "WEB_DEMO_README.md"

        self.assertTrue(script.exists(), "start_web.sh should provide one-command startup")
        self.assertTrue(readme.exists(), "WEB_DEMO_README.md should explain demo startup")
        script_text = script.read_text(encoding="utf-8")
        readme_text = readme.read_text(encoding="utf-8")
        health_check = root / "health_check.py"
        self.assertTrue(health_check.exists(), "health_check.py should summarize /api/system/health")
        health_text = health_check.read_text(encoding="utf-8")
        for expected in [
            "webapp/server.py",
            "/api/health",
            "/api/system/health",
            "localhost:8765",
            "health_check.py",
            "Existing service is not responding",
            "Existing service became ready",
            "Use a different port",
        ]:
            self.assertIn(expected, script_text)
        for expected in ["/api/system/health", "overall_status", "recommended_commands"]:
            self.assertIn(expected, health_text)
        for expected in ["start_web.sh", "/api/health", "/api/system/health", "/api/jobs", "完整交付导出", "产品指挥台"]:
            self.assertIn(expected, readme_text)


if __name__ == "__main__":
    unittest.main()
