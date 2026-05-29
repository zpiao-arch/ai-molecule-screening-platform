import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class TargetBriefFrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_target_brief_page_is_grouped_as_brief_form(self):
        for expected in [
            'id="target-brief-form-shell"',
            'id="brief-section-disease"',
            'id="brief-section-virus"',
            'id="brief-section-protein"',
            'id="brief-section-pocket"',
            'id="brief-section-reference-drug"',
            'id="brief-section-design-intent"',
            "疾病与适应症",
            "病毒/病原体",
            "蛋白与结构",
            "口袋与关键残基",
            "参考药/控药",
            "设计意图",
        ]:
            self.assertIn(expected, self.html)

    def test_target_brief_redesign_preserves_existing_stage1_controls(self):
        for expected in [
            'id="stage1-target-workspace"',
            'id="stage1-target-results"',
            'id="stage1-target-pack"',
            'id="stage1-target-pack-validation"',
            'id="stage1-prompt-preview"',
            'id="stage1-prompt-text"',
            'id="st1-disease"',
            'id="st1-target"',
            'id="st1-top"',
            'id="st1-source-kind"',
            'id="st1-target-hint"',
            'id="st1-source"',
            'id="bf-target-name"',
            'id="bf-disease"',
            'id="bf-protein"',
            'id="bf-pdb-id"',
            'id="bf-reference-ligand"',
            'id="bf-key-residues"',
            'id="bf-pocket"',
            'id="bf-must-have"',
            'id="bf-avoid"',
            'id="bf-force"',
            'id="bf-submit-btn"',
            "runStage1TargetSelect()",
            "runStage1TargetIntake()",
            "generateStage1BriefFromTarget()",
            "buildStage1TargetPack()",
            "runStage1TargetPackValidate()",
            "submitBrief()",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
