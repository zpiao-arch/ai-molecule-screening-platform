import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class Stage8FrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_sidebar_exposes_stage8_product_command_center(self):
        for expected in [
            'data-page="stage8"',
            'id="page-stage8"',
            "Stage 8 · 产品指挥台",
            "产品指挥台",
        ]:
            self.assertIn(expected, self.html)

    def test_stage8_page_exposes_command_center_regions(self):
        for expected in [
            'id="st8-project-select"',
            'id="st8-round"',
            'id="st8-demo-mode"',
            'id="st8-preflight-btn"',
            'id="st8-repair-btn"',
            'id="st8-acceptance-btn"',
            'id="st8-report-btn"',
            'id="st8-catalog-btn"',
            'id="st8-candidate-btn"',
            'id="stage8-acceptance-status"',
            'id="stage8-acceptance-report"',
            'id="st8-demo-package-btn"',
            'id="st8-demo-doctor-btn"',
            'id="st8-demo-runner-btn"',
            'id="st8-demo-guide-btn"',
            'id="st8-evidence-pack-btn"',
            'id="st8-review-mode-btn"',
            'id="stage8-review-status"',
            'id="stage8-review-mode"',
            'id="stage8-demo-guide-status"',
            'id="stage8-demo-guide"',
            'id="stage8-evidence-pack-status"',
            'id="stage8-evidence-pack"',
            'id="stage8-preflight-status"',
            'id="stage8-preflight"',
            'id="stage8-report-status"',
            'id="stage8-report-preview"',
            'id="stage8-target-catalog"',
            'id="stage8-candidate-entry"',
            'id="stage8-demo-runner-status"',
            'id="stage8-demo-runner"',
            'id="stage8-status"',
            'id="stage8-stage-rail"',
            'id="stage8-readiness"',
            'id="stage8-metrics"',
            'id="stage8-target-summary"',
            'id="stage8-next-actions"',
            'id="stage8-funnel"',
            'id="stage8-quality-gates"',
            'id="stage8-hit-triage"',
            'id="stage8-deliverables"',
            'id="stage8-downloads"',
            'id="stage8-demo-doctor"',
            'id="stage8-demo-doctor-status"',
            'id="stage8-boundary"',
        ]:
            self.assertIn(expected, self.html)
        self.assertNotIn('id="st8-demo-mode" checked', self.html)

    def test_stage8_frontend_calls_command_center_api_and_renders_payload(self):
        for expected in [
            "loadStage8",
            "loadStage8CommandCenter()",
            "renderStage8CommandCenter(",
            "renderStage8Rail(",
            "renderStage8Deliverables(",
            "renderStage8Downloads(",
            "renderStage8DemoDoctor(",
            "runStage8DemoPackage()",
            "runStage8DemoDoctor()",
            "runStage8DemoRunner()",
            "loadStage8DemoGuide()",
            "runStage8EvidencePack()",
            "loadStage8ReviewMode()",
            "renderStage8ReviewMode(",
            "renderStage8DemoGuide(",
            "renderStage8EvidencePack(",
            "loadStage8Preflight()",
            "runStage8Repair()",
            "runStage8AcceptanceDemo()",
            "generateStage8Report()",
            "loadStage8TargetCatalog()",
            "runStage8CandidateEntry()",
            "renderStage8Acceptance(",
            "renderStage8DemoRunner(",
            "renderStage8Preflight(",
            "renderStage8Report(",
            "renderStage8TargetCatalog(",
            "renderStage8CandidateEntry(",
            "selectStage8DefaultProject(",
            "openStage8Route(",
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/command-center?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/preflight?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/repair"',
            '"/api/stage8/acceptance-demo"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/report"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/demo-package"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/demo-runner"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/demo-guide?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/evidence-pack"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/review-mode?round="+encodeURIComponent(roundNo)',
            '"/api/target-catalog?query="+encodeURIComponent(q)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage3/candidates"',
            '"/api/projects/"+encodeURIComponent(n)+"/demo-doctor?round="+encodeURIComponent(roundNo)',
        ]:
            self.assertIn(expected, self.html)

    def test_stage8_frontend_exposes_deep_integration_panel(self):
        for expected in [
            'id="stage8-deep-integration-status"',
            'id="stage8-deep-integration"',
            "六项深化融合",
            "真实 Demo 加硬",
            "多生成器适配",
            "Stage 4 科学可信度",
            "Stage 8 总控",
            "交付复现",
            "科学报告解释层",
            "loadStage8DeepIntegration()",
            "renderStage8DeepIntegration(",
            '"/api/projects/"+encodeURIComponent(n)+"/product/deep-integration?round="+encodeURIComponent(roundNo)',
        ]:
            self.assertIn(expected, self.html)

    def test_stage8_frontend_explains_openai_api_key_safety(self):
        for expected in [
            "API Key 只随本次请求发送",
            "不会写入项目文件",
            "stage8-api-key-safety",
            "api_key_security",
            "api_key_persisted",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
