import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class Stage5FrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_dashboard_exposes_stage5_controls_and_regions(self):
        for expected in [
            'id="stage5-dashboard-workspace"',
            'id="st5-project-select"',
            'id="st5-round"',
            'id="st5-title"',
            'id="stage5-metrics"',
            'id="stage5-readiness"',
            'id="stage5-target-summary"',
            'id="stage5-ranked-table"',
            'id="stage5-boundary"',
            "刷新Stage 5",
            "生成Stage 5看板",
        ]:
            self.assertIn(expected, self.html)

    def test_stage5_frontend_calls_stage5_api_routes(self):
        for expected in [
            '"/api/projects/"+encodeURIComponent(n)+"/stage5?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage5/dashboard"',
            "loadStage5Dashboard()",
            "generateStage5Dashboard()",
            "syncStage5RoundFromProject(",
            "latestRoundForProject(",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
