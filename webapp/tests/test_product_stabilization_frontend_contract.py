import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class ProductStabilizationFrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_stage1_exposes_target_pack_validation(self):
        for expected in [
            'id="st1-target-pack-validate-btn"',
            'id="stage1-target-pack-validation"',
            "runStage1TargetPackValidate()",
            "renderStage1TargetPackValidation(",
            '"/api/projects/"+encodeURIComponent(n)+"/target-pack/validate?round="+encodeURIComponent(latestRoundForProject(n))',
        ]:
            self.assertIn(expected, self.html)

    def test_stage4_exposes_scientific_readiness(self):
        for expected in [
            'id="st4-scientific-readiness-btn"',
            'id="stage4-scientific-readiness-status"',
            'id="stage4-scientific-readiness"',
            "loadScientificReadiness()",
            "renderScientificReadiness(",
            '"/api/projects/"+encodeURIComponent(n)+"/scientific-readiness?round="+encodeURIComponent(roundNo)',
        ]:
            self.assertIn(expected, self.html)

    def test_stage8_exposes_full_export_and_job_queue(self):
        for expected in [
            'id="st8-full-export-btn"',
            'id="st8-full-export-job-btn"',
            'id="stage8-full-export-status"',
            'id="stage8-full-export"',
            'id="stage8-jobs"',
            "runStage8FullExport()",
            "runStage8FullExportJob()",
            "loadProductJobs()",
            "function renderStage8FullExport(",
            "function renderProductJobs(",
            '"/api/projects/"+encodeURIComponent(n)+"/stage8/full-export"',
            '"/api/jobs"',
            '"/api/jobs/"+encodeURIComponent(jobId)',
        ]:
            self.assertIn(expected, self.html)

    def test_system_page_uses_product_health_endpoint(self):
        for expected in [
            'id="sys-health-status"',
            'id="sys-health-sections"',
            'id="sys-jobs"',
            "loadSystemHealth()",
            "function loadSystemHealth()",
            "function renderSystemHealth(",
            "function renderSystemHealthSections(",
            '"/api/system/health"',
            '"/api/jobs"',
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
