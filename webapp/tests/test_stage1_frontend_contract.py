import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class Stage1FrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_brief_page_exposes_stage1_target_selection_controls(self):
        for expected in [
            'id="stage1-target-workspace"',
            'id="st1-disease"',
            'id="st1-target"',
            'id="st1-top"',
            'id="st1-source-kind"',
            'id="st1-source"',
            'id="st1-target-hint"',
            'id="st1-target-pack-btn"',
            'id="stage1-target-pack"',
            'id="stage1-target-results"',
            'runStage1TargetSelect()',
            'runStage1TargetIntake()',
            'buildStage1TargetPack()',
            'loadStage1Status()',
        ]:
            self.assertIn(expected, self.html)

    def test_stage1_frontend_calls_stage1_api_routes(self):
        for expected in [
            '"/api/projects/"+encodeURIComponent(n)+"/stage1"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage1/target-select"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage1/brief-from-target"',
            '"/api/projects/"+encodeURIComponent(n)+"/target-intake"',
            '"/api/projects/"+encodeURIComponent(n)+"/target-pack"',
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
