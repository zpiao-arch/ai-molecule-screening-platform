import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class DashboardFrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_dashboard_is_project_status_home(self):
        for expected in [
            'id="dashboard-overview-shell"',
            'id="dashboard-status-hero"',
            "项目状态首页",
            "项目概况、Stage 5 看板和近期处理记录",
            'id="dashboard-stats"',
            "dashboard-metrics-grid",
        ]:
            self.assertIn(expected, self.html)

    def test_dashboard_keeps_stage5_workspace_with_denser_toolbar(self):
        for expected in [
            'id="stage5-dashboard-workspace"',
            'id="dashboard-stage5-toolbar"',
            'id="stage5-status"',
            'id="stage5-metrics"',
            'id="stage5-readiness"',
            'id="stage5-target-summary"',
            'id="stage5-boundary"',
            'id="stage5-ranked-table"',
            "刷新Stage 5",
            "生成Stage 5看板",
        ]:
            self.assertIn(expected, self.html)

    def test_dashboard_recent_projects_render_status_rows(self):
        for expected in [
            'id="dashboard-recent"',
            'id="dashboard-recent-list"',
            "project-status-row",
            "status-dot",
            "renderDashboardStats(",
            "renderDashboardRecentProjects(",
            "overviewMetricCard(",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
