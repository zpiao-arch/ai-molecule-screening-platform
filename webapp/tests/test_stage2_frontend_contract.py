import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class Stage2FrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_sidebar_exposes_stage2_evidence_page(self):
        for expected in [
            'data-page="stage2"',
            'id="page-stage2"',
            "靶点证据",
            "Stage 2 · 证据矩阵",
        ]:
            self.assertIn(expected, self.html)

    def test_stage2_page_exposes_evidence_controls_and_outputs(self):
        for expected in [
            'id="st2-project-select"',
            'id="st2-disease"',
            'id="st2-target"',
            'id="st2-top"',
            'id="st2-offline"',
            'id="stage2-evidence-results"',
            'id="stage2-evidence-detail"',
            'id="stage2-report-preview"',
            "runStage2Evidence()",
            "loadStage2Status()",
        ]:
            self.assertIn(expected, self.html)

    def test_stage2_frontend_calls_stage2_api_routes(self):
        for expected in [
            '"/api/projects/"+encodeURIComponent(n)+"/stage2"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage2/evidence"',
        ]:
            self.assertIn(expected, self.html)

    def test_stage2_evidence_matrix_exposes_strength_markers(self):
        for expected in [
            'id="stage2-evidence-strength-legend"',
            "证据强弱",
            "强证据",
            "中等证据",
            "待补证据",
            "stage2EvidenceStrength(",
            "renderStage2EvidenceStrengthBadge(",
            '"evidence_strength"',
            "evidence-strength-badge",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
