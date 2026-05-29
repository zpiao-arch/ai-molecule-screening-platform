import unittest
from pathlib import Path


INDEX = Path(__file__).resolve().parents[1] / "static" / "index.html"


class ProjectRegistryFrontendContractTests(unittest.TestCase):
    def setUp(self):
        self.html = INDEX.read_text(encoding="utf-8")

    def test_project_registry_exposes_status_aware_list(self):
        for expected in [
            'id="page-projects"',
            'id="projects-list"',
            "项目状态",
            "候选",
            "受体",
            "Docking",
            "project-registry-list",
            "project-asset-pill",
        ]:
            self.assertIn(expected, self.html)

    def test_project_registry_renders_asset_status_from_api(self):
        for expected in [
            "renderProjectRegistry(",
            "projectAssetPill(",
            "projectAssetStatus(",
            "asset_status",
            "candidate",
            "receptor",
            "docking",
        ]:
            self.assertIn(expected, self.html)


if __name__ == "__main__":
    unittest.main()
