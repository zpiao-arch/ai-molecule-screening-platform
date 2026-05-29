import argparse
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "ai_mol_loop.py"
SPEC = importlib.util.spec_from_file_location("ai_mol_loop_module", MODULE_PATH)
assert SPEC and SPEC.loader
ai_mol_loop = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ai_mol_loop)


def rows_from_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


class Stage5DashboardTests(unittest.TestCase):
    def test_stage5_dashboard_writes_product_assets_from_prior_stage_outputs(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            ai_mol_loop.ensure_project_dirs(project)
            ai_mol_loop.write_json(project / "config.json", ai_mol_loop.default_config())
            ai_mol_loop.write_json(
                project / "briefs" / "target_brief.json",
                {
                    "target_catalog_id": "influenza_a_h1n1_na",
                    "target_name": "Influenza A neuraminidase",
                    "disease_context": "甲流",
                    "binding_site": {"reference_ligand": "oseltamivir", "key_residues": ["Arg118", "Glu276"]},
                },
            )
            ai_mol_loop.write_csv(
                project / "evidence" / "stage2_target_sources.csv",
                [
                    {
                        "target_id": "influenza_a_h1n1_na",
                        "target": "Influenza A neuraminidase",
                        "evidence_score": "0.86",
                        "readiness": "ready_for_mvp_closed_loop",
                        "primary_pdb": "3TI6",
                        "positive_controls": "oseltamivir; zanamivir",
                        "pubmed_articles": "5",
                    }
                ],
                ["target_id", "target", "evidence_score", "readiness", "primary_pdb", "positive_controls", "pubmed_articles"],
            )
            ai_mol_loop.write_json(
                project / "evidence" / "stage2_closed_loop_assets.json",
                {
                    "targets": [
                        {
                            "target_id": "influenza_a_h1n1_na",
                            "target": "Influenza A neuraminidase",
                            "evidence_quality": {"readiness": "ready_for_mvp_closed_loop", "score": 0.86},
                        }
                    ]
                },
            )
            ai_mol_loop.write_json(
                project / "stage3" / "round_1_stage3_assets.json",
                {"round": 1, "raw_candidates": 3, "passed_filter": 2, "failed_filter": 1, "sources": ["unit"]},
            )
            ai_mol_loop.write_csv(
                project / "candidates" / "round_1_candidates.csv",
                [
                    {"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"},
                    {"round": 1, "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "parent": "", "source": "unit"},
                ],
                ["round", "id", "smiles", "parent", "source"],
            )
            ai_mol_loop.write_csv(
                project / "filtered" / "round_1_filtered.csv",
                [
                    {"id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "passed_filter": "1"},
                    {"id": "bad_001", "smiles": "not_a_smiles", "passed_filter": "0"},
                ],
                ["id", "smiles", "passed_filter"],
            )
            ai_mol_loop.write_csv(
                project / "ranked" / "round_1_ranked.csv",
                [
                    {"rank": "1", "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "total_proxy": "0.72", "decision": "advance"},
                    {"rank": "2", "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "total_proxy": "0.61", "decision": "hold"},
                ],
                ["rank", "id", "smiles", "total_proxy", "decision"],
            )
            ai_mol_loop.write_csv(
                project / "scores" / "round_1_scores.csv",
                [
                    {"id": "mol_001", "score_source": "proxy+stage4_rdkit"},
                    {"id": "mol_002", "score_source": "proxy+stage4_rdkit"},
                ],
                ["id", "score_source"],
            )
            ai_mol_loop.write_json(
                project / "feedback" / "round_1_feedback.json",
                {"round": 1, "selected_count": 1, "selected_smiles": ["CCOC(=O)c1ccccc1"]},
            )
            ai_mol_loop.write_csv(
                project / "seeds" / "round_1_seeds.csv",
                [{"id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "note": "advance"}],
                ["id", "smiles", "note"],
            )
            ai_mol_loop.write_csv(
                project / "stage4" / "round_1_real_descriptors.csv",
                [
                    {"id": "mol_001", "canonical_smiles": "CCOC(=O)c1ccccc1", "valid": "1", "qed": "0.65", "mw": "150.0", "lipinski_violations": "0"},
                    {"id": "mol_002", "canonical_smiles": "CCS(=O)(=O)Nc1ccccc1", "valid": "1", "qed": "0.58", "mw": "199.0", "lipinski_violations": "0"},
                ],
                ["id", "canonical_smiles", "valid", "qed", "mw", "lipinski_violations"],
            )
            ai_mol_loop.write_csv(
                project / "stage4" / "round_1_benchmark_panel.csv",
                [
                    {"panel_type": "candidate", "id": "mol_001", "canonical_smiles": "CCOC(=O)c1ccccc1", "qed": "0.65", "nearest_control": "oseltamivir"},
                    {"panel_type": "positive_control", "id": "oseltamivir", "canonical_smiles": "CC", "qed": "0.55", "nearest_control": "oseltamivir"},
                    {"panel_type": "decoy", "id": "decoy_001", "canonical_smiles": "CCO", "qed": "0.42", "nearest_control": "oseltamivir"},
                ],
                ["panel_type", "id", "canonical_smiles", "qed", "nearest_control"],
            )
            ai_mol_loop.write_json(
                project / "stage4" / "round_1_docking_plan.json",
                {"status": "skipped", "selected_backend": "", "message": "No backend installed."},
            )
            ai_mol_loop.write_json(
                project / "stage4" / "round_1_validation_metrics.json",
                {
                    "panel_counts": {"candidate": 2, "positive_control": 1, "decoy": 1},
                    "readiness": {"rdkit_descriptors": "ready", "control_panel": "ready", "decoy_panel": "ready", "docking": "skipped"},
                    "docking_status": "skipped",
                },
            )
            image_dir = project / "stage4" / "round_1_2d"
            image_dir.mkdir(parents=True)
            (image_dir / "candidate_mol_001.png").write_bytes(b"not-a-real-png-but-path-exists")
            ai_mol_loop.write_json(
                project / "stage4" / "round_1_stage4_assets.json",
                {
                    "stage": 4,
                    "round": 1,
                    "target_id": "influenza_a_h1n1_na",
                    "candidate_count": 2,
                    "valid_count": 2,
                    "controls_count": 1,
                    "decoy_count": 1,
                    "rendered_2d_images": [{"panel_type": "candidate", "id": "mol_001", "path": str(image_dir / "candidate_mol_001.png")}],
                    "files": {"real_descriptors": str(project / "stage4" / "round_1_real_descriptors.csv")},
                    "docking_plan": {"status": "skipped"},
                },
            )

            ai_mol_loop.stage5_dashboard(argparse.Namespace(project=str(project), round=1, title="Unit Dashboard"))

            data_path = project / "stage5" / "dashboard_data.json"
            html_path = project / "stage5" / "index.html"
            css_path = project / "stage5" / "styles.css"
            js_path = project / "stage5" / "app.js"
            report_path = project / "reports" / "stage5_dashboard_report.md"

            self.assertTrue(data_path.exists())
            self.assertTrue(html_path.exists())
            self.assertTrue(css_path.exists())
            self.assertTrue(js_path.exists())
            self.assertTrue(report_path.exists())

            data = json.loads(data_path.read_text(encoding="utf-8"))
            self.assertEqual(data["stage"], 5)
            self.assertEqual(data["round"], 1)
            self.assertEqual(data["project"]["name"], "project")
            self.assertEqual(data["target"]["id"], "influenza_a_h1n1_na")
            self.assertEqual(data["readiness"]["rdkit_validation"], "ready")
            self.assertEqual(data["readiness"]["docking"], "skipped")
            self.assertEqual(data["metrics"]["raw_candidates"], 3)
            self.assertEqual(data["metrics"]["filtered_candidates"], 2)
            self.assertEqual(data["metrics"]["valid_rdkit"], 2)
            self.assertEqual(data["metrics"]["advanced"], 1)
            self.assertIn("ranked_top", data["tables"])
            self.assertIn("benchmark_panel", data["tables"])
            self.assertIn("validation_metrics", data["tables"])
            self.assertEqual(data["molecule_images"][0]["relative_path"], "../stage4/round_1_2d/candidate_mol_001.png")

            html = html_path.read_text(encoding="utf-8")
            self.assertIn("dashboard_data.json", html)
            self.assertIn("styles.css", html)
            self.assertIn("app.js", html)
            self.assertIn("dashboard-data", html)
            self.assertIn("Unit Dashboard", html)

            report = report_path.read_text(encoding="utf-8")
            self.assertIn("Stage 5 Product Dashboard", report)
            self.assertIn("Computational screening", report)


if __name__ == "__main__":
    unittest.main()
