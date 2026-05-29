import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class Stage6Stage7FrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_sidebar_exposes_stage6_and_stage7_pages(self):
        for expected in [
            'data-page="stage6"',
            'data-page="stage7"',
            'id="page-stage6"',
            'id="page-stage7"',
            "验证运营",
            "交付包",
        ]:
            self.assertIn(expected, self.html)

    def test_stage6_page_exposes_controls_tables_and_report_regions(self):
        for expected in [
            "Stage 6 · 验证运营",
            'id="st6-project-select"',
            'id="st6-round"',
            'id="st6-top"',
            'id="stage6-status"',
            'id="stage6-summary"',
            'id="stage6-gates-table"',
            'id="stage6-triage-table"',
            'id="stage6-queue-table"',
            'id="stage6-risk-table"',
            'id="stage6-runbook-preview"',
            'id="stage6-report-preview"',
            "runStage6Validation()",
            "loadStage6Status()",
            "renderStage6Status(",
        ]:
            self.assertIn(expected, self.html)

    def test_stage7_page_exposes_controls_deliverables_and_previews(self):
        for expected in [
            "Stage 7 · 交付包",
            'id="st7-project-select"',
            'id="st7-round"',
            'id="st7-title"',
            'id="stage7-status"',
            'id="stage7-deliverables-table"',
            'id="stage7-checklist-table"',
            'id="stage7-summary-preview"',
            'id="stage7-runbook-preview"',
            'id="stage7-stage8-preview"',
            'id="stage7-report-preview"',
            "generateStage7Package()",
            "loadStage7Status()",
            "renderStage7Status(",
        ]:
            self.assertIn(expected, self.html)

    def test_stage6_stage7_frontend_calls_api_routes(self):
        for expected in [
            '"/api/projects/"+encodeURIComponent(n)+"/stage6?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage6/validate"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage7?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage7/package"',
            "loadStage6",
            "loadStage7",
            "renderStage6Table(",
            "renderStage7Deliverables(",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
