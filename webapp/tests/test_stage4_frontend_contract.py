import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"
STATIC = INDEX.parent


class Stage4FrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_sidebar_exposes_stage4_real_library_page(self):
        for expected in [
            'data-page="stage4"',
            'id="page-stage4"',
            "真实库校验",
            "Stage 4 · 真实化学校验",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_page_exposes_controls_and_outputs(self):
        for expected in [
            'id="st4-project-select"',
            'id="st4-round"',
            'id="st4-target"',
            'id="st4-preset-select"',
            'id="st4-top"',
            'id="st4-decoys"',
            'id="st4-max-conformers"',
            'id="st4-receptor-pdb"',
            'id="st4-docking-backend"',
            'id="st4-pocket-center"',
            'id="st4-pocket-size"',
            'id="st4-pocket-source"',
            'id="st4-fetch-receptor"',
            'id="st4-run-docking"',
            'id="st4-docking-timeout"',
            'id="st4-external-scores"',
            'id="st4-rescore"',
            "runStage4SmokeTest()",
            "repairStage4Project()",
            'id="stage4-status"',
            'id="stage4-capability-status"',
            'id="stage4-metrics"',
            'id="stage4-readiness"',
            'id="stage4-preflight"',
            'id="stage4-project-doctor"',
            'id="stage4-library-status"',
            'id="stage4-receptor-summary"',
            'id="stage4-wizard-guide"',
            'id="stage4-docking-plan"',
            'id="stage4-structure-viewer"',
            'id="stage4-structure-status"',
            'id="stage4-pose-select"',
            'id="stage4-viewer-meta"',
            'id="stage4-docking-results"',
            'id="stage4-artifact-links"',
            'id="stage4-image-strip"',
            'id="stage4-descriptors-table"',
            'id="stage4-similarity-table"',
            'id="stage4-diverse-table"',
            'id="stage4-benchmark-table"',
            'id="stage45-status"',
            'id="stage45-metrics"',
            'id="stage45-control-summary"',
            'id="stage45-scores-table"',
            'id="stage45-report-preview"',
            'id="stage46-status"',
            'id="stage46-metrics"',
            'id="stage46-topk"',
            'id="stage46-ranking-table"',
            'id="stage46-report-preview"',
            'id="stage4-report-preview"',
            "runStage45Validation()",
            "runStage46Benchmark()",
            "runStage4Real()",
            "loadStage4Status()",
            "syncStage4RoundFromProject(",
            "latestRoundForProject(",
            "applyStage4Preset()",
        ]:
            self.assertIn(expected, self.html)
        self.assertNotIn('id="st4-round" value="2"', self.html)

    def test_stage4_frontend_calls_stage4_api_routes(self):
        for expected in [
            '"/api/projects/"+encodeURIComponent(n)+"/stage4?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage4/real"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage4/repair"',
            '"/api/projects/"+encodeURIComponent(n)+"/doctor"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage45?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage45/validate"',
            '"/api/projects/"+encodeURIComponent(n)+"/stage46?round="+encodeURIComponent(roundNo)',
            '"/api/projects/"+encodeURIComponent(n)+"/stage46/benchmark"',
            '"/api/stage4/capabilities"',
            '"/api/stage4/presets"',
            '"/api/stage4/smoke-test"',
            "loadStage4Capabilities(",
            "loadStage4Presets(",
            "renderStage4Capabilities(",
            "fetch_receptor:",
            "receptor_pdb:",
            "run_docking:",
            "docking_timeout:",
            "external_scores:",
            "pocket_center:",
            "pocket_size:",
            "pocket_source:",
            "renderStage4Status(",
            "renderStage4Preflight(",
            "renderStage4ProjectDoctor(",
            "renderStage4Wizard(",
            "renderStage4DockingResults(",
            "renderStage4StructureViewer(",
            "loadStage4StructurePose(",
            "initStage4Viewer(",
            "3Dmol",
            "visualization_assets",
            "renderStage4Artifacts(",
            "renderStage4Images(",
            "renderStage45Status(",
            "renderStage46Status(",
            "syncStage4RoundFromProject(",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_3dmol_viewer_uses_local_vendor_asset(self):
        self.assertIn('src="/static/vendor/3Dmol-min.js"', self.html)
        self.assertNotIn("cdn.jsdelivr.net/npm/3dmol", self.html)
        vendor = STATIC / "vendor" / "3Dmol-min.js"
        self.assertTrue(vendor.exists())
        self.assertGreater(vendor.stat().st_size, 100_000)

    def test_stage4_viewer_keeps_missing_pose_message_on_reset(self):
        for expected in [
            "stage4ViewerMissingPoseMessage(",
            'assets.status==="missing_pose"',
            'a.status==="missing_pose"',
            "缺少 docking pose，正在显示受体结构",
            "已加载 receptor-only 结构视图",
            "loadStage4ReceptorOnlyView(",
            "需要真实 Vina/GNINA 生成 pose",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_viewer_can_render_receptor_without_pose(self):
        for expected in [
            "loadStage4ReceptorOnlyView(",
            "已加载 receptor-only 结构视图",
            "缺少 docking pose，已显示受体结构。",
            "fetch(rec.url)",
            "addStage4DockingBox(stage4Viewer,a.docking_box||{})",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_viewer_can_overlay_reference_ligand_without_pose(self):
        for expected in [
            "stage4ViewerReferenceLigandMessage(",
            "a.reference_ligand",
            "fetch(ref.url)",
            "referenceModel.setStyle",
            "已加载受体与共晶参考配体",
            "不是候选 docking pose",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_viewer_has_visible_projection_fallback(self):
        for expected in [
            "renderStage4StructureFallback(",
            "stage4ParsePdbAtoms(",
            "stage4DrawProjectionFallback(",
            "stage4Render3DmolOrFallback(",
            "stage4-3dmol-layer",
            "结构投影兜底视图",
            "proteinAtoms",
            "ligandAtoms",
        ]:
            self.assertIn(expected, self.html)
        self.assertNotIn('overlay.style.display="none"', self.html)

    def test_stage4_projection_fallback_focuses_on_ligand_pocket(self):
        for expected in [
            "stage4LigandCenter(",
            "stage4FilterPocketAtoms(",
            "stage4InferPocketRadius(",
            "pocketAtoms",
            "口袋局部投影",
            "stage4-pocket-shell",
            "stage4-ligand-focus",
            "参考配体",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_viewer_promotes_3dmol_over_projection_when_available(self):
        for expected in [
            "stage4-viewer-3d-ready",
            "stage4Clear3DReady(",
            "stage4Mark3DReady(",
            "stage4ReferenceSelection(",
            "口袋参考视图",
            "zoomTo(referenceSelection)",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_hash_route_opens_stage4_page(self):
        for expected in [
            "activatePageFromHash(",
            "stage4-section-run",
            "stage4-section-structure",
            "window.addEventListener(\"hashchange\"",
            'activatePage("stage4"',
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_exposes_interactive_guide(self):
        for expected in [
            'id="stage4-interactive-guide"',
            'id="stage4-guide-steps"',
            'id="stage4-guide-next-action"',
            "Stage 4 交互式引导",
            "运行前检查",
            "受体/口袋",
            "Docking/分数",
            "对照校准",
            "产物查看",
            "下一步",
            "stage4GuideSteps(",
            "renderStage4InteractiveGuide(",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_frontend_surfaces_auto_repair_hints(self):
        for expected in [
            "auto_repair",
            "repair_hint",
            "repair_endpoint",
            "可修复",
            "showStage4RepairSummary(",
            "Stage 4资产修复完成",
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_refresh_preserves_manual_round(self):
        self.assertIn("syncStage4RoundFromProject(currentStage4Project(), true)", self.html)
        self.assertIn("function loadStage4Status(){var n=currentStage4Project();var roundNo=stage4Round();", self.html)
        self.assertNotIn("function loadStage4Status(){var n=currentStage4Project();var roundNo=stage4Round();if(!n){setStage4Status(\"请选择项目\",\"warn\");renderStage45Status(null);renderStage46Status(null);renderStage4ProjectDoctor(null);return;}syncStage4RoundFromProject(n);", self.html)


if __name__ == "__main__":
    unittest.main()
