import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class Stage5To12FrontendPolishContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_stage4_exposes_internal_anchor_navigation_without_hiding_tools(self):
        for expected in [
            'id="stage4-anchor-nav"',
            'href="#stage4-section-run"',
            'href="#stage4-section-structure"',
            'href="#stage4-section-controls"',
            'href="#stage4-section-artifacts"',
            'id="stage4-section-run"',
            'id="stage4-section-structure"',
            'id="stage4-section-controls"',
            'id="stage4-section-artifacts"',
            "运行与依赖",
            "受体/口袋/3D",
            "对照校准",
            "产物与表格",
        ]:
            self.assertIn(expected, self.html)

    def test_stage6_is_framed_as_validation_plan_page(self):
        for expected in [
            'id="stage6-validation-plan-shell"',
            'id="stage6-validation-plan-steps"',
            "验证计划页",
            "实验验证候选",
            "质量门",
            "风险复核",
        ]:
            self.assertIn(expected, self.html)

    def test_stage7_prioritizes_delivery_generation_and_download_materials(self):
        for expected in [
            'id="stage7-delivery-actions"',
            'id="stage7-material-lane"',
            "交付材料生成页",
            "生成交付包",
            "下载材料",
            "复现说明",
        ]:
            self.assertIn(expected, self.html)

    def test_stage8_highlights_primary_product_actions(self):
        for expected in [
            'id="stage8-primary-actions"',
            'id="stage8-product-focus"',
            "一键演示",
            "一键导出",
            "缺口修复",
            "聚合状态",
        ]:
            self.assertIn(expected, self.html)

    def test_workflow_results_config_system_get_product_framing(self):
        for expected in [
            'id="workflow-closed-loop-map"',
            "闭环运行",
            "Stage 1 到 Stage 8",
            'id="results-candidate-pool-guide"',
            "最终候选池",
            "证据解释",
            'id="config-basic-form-guide"',
            "普通表单优先",
            "高级 JSON",
            'id="system-dependency-diagnostics"',
            "依赖诊断",
            "RDKit",
            "Vina",
            "OpenBabel",
            "Meeko",
            "PoseBusters",
            "3Dmol",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
