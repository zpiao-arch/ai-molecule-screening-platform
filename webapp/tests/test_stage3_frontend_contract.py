import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class Stage3FrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_sidebar_exposes_stage3_candidate_page(self):
        for expected in [
            'data-page="stage3"',
            'id="page-stage3"',
            "Stage 3 · 候选生成与初筛",
            "候选生成",
        ]:
            self.assertIn(expected, self.html)

    def test_stage3_page_exposes_candidate_generation_controls(self):
        for expected in [
            'id="st3-project-select"',
            'id="st3-round"',
            'id="st3-source-mode"',
            'id="st3-n"',
            'id="st3-top"',
            'id="st3-source"',
            'id="st3-use-openai"',
            'id="st3-openai-model"',
            'id="st3-api-key"',
            'id="stage3-api-key-safety"',
            'id="st3-prompt"',
            'id="st3-no-score"',
            'id="st3-run-btn"',
            'id="stage3-status"',
            'id="stage3-metrics"',
            'id="stage3-context-summary"',
            'id="stage3-raw-table"',
            'id="stage3-candidates-table"',
            'id="stage3-ranked-table"',
            'id="stage3-feedback-panel"',
            'id="stage3-prompt-preview"',
            'id="stage3-report-preview"',
        ]:
            self.assertIn(expected, self.html)

    def test_stage3_frontend_calls_stage3_api_and_renders_payload(self):
        for expected in [
            "loadStage3",
            "currentStage3Project()",
            "syncStage3RoundFromProject(",
            "runStage3Candidates()",
            "loadStage3Status()",
            "renderStage3Status(",
            "renderStage3CandidateEntry(",
            '"/api/projects/"+encodeURIComponent(n)+"/stage3?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage3/candidates"',
            "use_openai:",
            "openai_model:",
            "api_key:",
            "source_mode:",
            "source_text:",
            "prompt:",
            "no_score:",
            "raw_candidates",
            "stage3_assets",
            "prompt_file",
            "uses_feedback_context",
        ]:
            self.assertIn(expected, self.html)

    def test_stage3_page_explains_feedback_context_and_boundary(self):
        for expected in [
            "上一轮反馈",
            "Prompt / Context",
            "读取 Stage 1 + Stage 2",
            "读取上一轮闭环反馈",
            "计算筛选结果，不代表实验活性",
            "API Key 只随本次请求发送",
            "不会写入项目文件",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
