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


def make_productized_project(project: Path) -> None:
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
                "assay_path": "NA docking -> pose check -> neuraminidase inhibition assay",
            }
        ],
        ["target_id", "target", "evidence_score", "readiness", "primary_pdb", "positive_controls", "assay_path"],
    )
    ai_mol_loop.write_json(
        project / "stage3" / "round_2_stage3_assets.json",
        {"round": 2, "raw_candidates": 4, "passed_filter": 3, "failed_filter": 1, "sources": ["unit"]},
    )
    ai_mol_loop.write_csv(
        project / "candidates" / "round_2_candidates.csv",
        [
            {"round": 2, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"},
            {"round": 2, "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "parent": "", "source": "unit"},
            {"round": 2, "id": "mol_003", "smiles": "COc1ccccc1O", "parent": "", "source": "unit"},
        ],
        ["round", "id", "smiles", "parent", "source"],
    )
    ai_mol_loop.write_csv(
        project / "scores" / "round_2_scores.csv",
        [
            {"id": "mol_001", "score_source": "proxy+stage4_rdkit", "total_proxy": "0.74"},
            {"id": "mol_002", "score_source": "proxy+stage4_rdkit", "total_proxy": "0.68"},
            {"id": "mol_003", "score_source": "proxy+stage4_rdkit", "total_proxy": "0.59"},
        ],
        ["id", "score_source", "total_proxy"],
    )
    ai_mol_loop.write_csv(
        project / "ranked" / "round_2_ranked.csv",
        [
            {"rank": "1", "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "total_proxy": "0.74", "qed_proxy": "0.65", "decision": "advance", "score_source": "proxy+stage4_rdkit"},
            {"rank": "2", "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "total_proxy": "0.68", "qed_proxy": "0.58", "decision": "advance", "score_source": "proxy+stage4_rdkit"},
            {"rank": "3", "id": "mol_003", "smiles": "COc1ccccc1O", "total_proxy": "0.59", "qed_proxy": "0.52", "decision": "hold", "score_source": "proxy+stage4_rdkit"},
        ],
        ["rank", "id", "smiles", "total_proxy", "qed_proxy", "decision", "score_source"],
    )
    ai_mol_loop.write_json(
        project / "feedback" / "round_2_feedback.json",
        {"round": 2, "selected_count": 2, "selected_smiles": ["CCOC(=O)c1ccccc1", "CCS(=O)(=O)Nc1ccccc1"]},
    )
    ai_mol_loop.write_csv(
        project / "seeds" / "round_2_seeds.csv",
        [
            {"id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "note": "advance"},
            {"id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "note": "advance"},
        ],
        ["id", "smiles", "note"],
    )
    ai_mol_loop.write_csv(
        project / "stage4" / "round_2_real_descriptors.csv",
        [
            {"id": "mol_001", "canonical_smiles": "CCOC(=O)c1ccccc1", "valid": "1", "qed": "0.65", "mw": "150.0", "lipinski_violations": "0"},
            {"id": "mol_002", "canonical_smiles": "CCS(=O)(=O)Nc1ccccc1", "valid": "1", "qed": "0.58", "mw": "199.0", "lipinski_violations": "0"},
            {"id": "mol_003", "canonical_smiles": "COc1ccccc1O", "valid": "1", "qed": "0.52", "mw": "124.0", "lipinski_violations": "0"},
        ],
        ["id", "canonical_smiles", "valid", "qed", "mw", "lipinski_violations"],
    )
    ai_mol_loop.write_csv(
        project / "stage4" / "round_2_benchmark_panel.csv",
        [
            {"panel_type": "candidate", "id": "mol_001", "canonical_smiles": "CCOC(=O)c1ccccc1", "qed": "0.65", "nearest_control": "oseltamivir"},
            {"panel_type": "positive_control", "id": "oseltamivir", "canonical_smiles": "CC", "qed": "0.55", "nearest_control": "oseltamivir"},
            {"panel_type": "decoy", "id": "decoy_001", "canonical_smiles": "CCO", "qed": "0.42", "nearest_control": "oseltamivir"},
        ],
        ["panel_type", "id", "canonical_smiles", "qed", "nearest_control"],
    )
    ai_mol_loop.write_json(
        project / "stage4" / "round_2_docking_plan.json",
        {"status": "skipped", "selected_backend": "", "expected_scores_csv": str(project / "stage4" / "round_2_docking_scores_template.csv")},
    )
    ai_mol_loop.write_json(
        project / "stage4" / "round_2_validation_metrics.json",
        {
            "panel_counts": {"candidate": 3, "positive_control": 1, "decoy": 1},
            "readiness": {"rdkit_descriptors": "ready", "control_panel": "ready", "decoy_panel": "ready", "docking": "skipped"},
            "docking_status": "skipped",
        },
    )
    ai_mol_loop.write_json(
        project / "stage4" / "round_2_stage4_assets.json",
        {
            "stage": 4,
            "round": 2,
            "target_id": "influenza_a_h1n1_na",
            "candidate_count": 3,
            "valid_count": 3,
            "controls_count": 1,
            "decoy_count": 1,
            "docking_plan": {"status": "skipped"},
        },
    )
    ai_mol_loop.stage5_dashboard(argparse.Namespace(project=str(project), round=2, title="Unit Stage 5"))


class Stage6Stage7ProductizationTests(unittest.TestCase):
    def test_stage6_validate_writes_quality_gates_triage_and_validation_plan(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            make_productized_project(project)

            ai_mol_loop.stage6_validate(argparse.Namespace(project=str(project), round=2, top=2))

            assets_path = project / "stage6" / "round_2_validation_assets.json"
            gates_path = project / "stage6" / "round_2_quality_gates.csv"
            triage_path = project / "stage6" / "round_2_hit_triage.csv"
            queue_path = project / "stage6" / "round_2_assay_queue.csv"
            report_path = project / "reports" / "stage6_round_2_validation_report.md"

            self.assertTrue(assets_path.exists())
            self.assertTrue(gates_path.exists())
            self.assertTrue(triage_path.exists())
            self.assertTrue(queue_path.exists())
            self.assertTrue(report_path.exists())

            assets = json.loads(assets_path.read_text(encoding="utf-8"))
            self.assertEqual(assets["stage"], 6)
            self.assertEqual(assets["round"], 2)
            self.assertIn(assets["overall_status"], {"ready_for_computational_demo", "computational_demo_ready_real_docking_missing", "blocked"})
            self.assertEqual(assets["docking_status"], "skipped")
            self.assertIn("quality_gate_summary", assets)
            self.assertIn("next_actions", assets)

            gates = rows_from_csv(gates_path)
            self.assertTrue(any(row["gate_id"] == "target_evidence" and row["status"] == "pass" for row in gates))
            self.assertTrue(any(row["gate_id"] == "real_docking" and row["status"] == "warn" for row in gates))

            triage = rows_from_csv(triage_path)
            self.assertEqual(len(triage), 2)
            self.assertTrue(all(row["validation_tier"] == "tier_1_structure_validation" for row in triage))
            self.assertTrue(all("docking" in row["next_action"] for row in triage))

            queue = rows_from_csv(queue_path)
            self.assertTrue(any(row["queue_type"] == "computational_docking" for row in queue))
            self.assertTrue(any(row["queue_type"] == "wet_lab_assay_planning" for row in queue))

            report = report_path.read_text(encoding="utf-8")
            self.assertIn("Stage 6 Validation Operations", report)
            self.assertIn("not efficacy proof", report)

    def test_stage7_package_writes_delivery_manifest_repro_and_stage8_frontend_spec(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            make_productized_project(project)
            ai_mol_loop.stage6_validate(argparse.Namespace(project=str(project), round=2, top=2))

            ai_mol_loop.stage7_package(argparse.Namespace(project=str(project), round=2, title="Unit Delivery"))

            manifest_path = project / "stage7" / "round_2_delivery_manifest.json"
            summary_path = project / "stage7" / "round_2_executive_summary.md"
            repro_path = project / "stage7" / "round_2_reproducibility.md"
            checklist_path = project / "stage7" / "round_2_investor_demo_checklist.csv"
            stage8_spec_path = project / "stage7" / "stage8_frontend_product_spec.md"
            report_path = project / "reports" / "stage7_round_2_delivery_report.md"

            self.assertTrue(manifest_path.exists())
            self.assertTrue(summary_path.exists())
            self.assertTrue(repro_path.exists())
            self.assertTrue(checklist_path.exists())
            self.assertTrue(stage8_spec_path.exists())
            self.assertTrue(report_path.exists())

            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["stage"], 7)
            self.assertEqual(manifest["round"], 2)
            self.assertIn("stage5_dashboard", manifest["deliverables"])
            self.assertIn("stage6_validation_assets", manifest["deliverables"])
            self.assertIn("stage8_frontend_spec", manifest["deliverables"])
            self.assertIn("claims_boundary", manifest)
            self.assertTrue(any(item["exists"] for item in manifest["deliverables"].values()))

            summary = summary_path.read_text(encoding="utf-8")
            self.assertIn("Unit Delivery", summary)
            self.assertIn("computational screening", summary)

            repro = repro_path.read_text(encoding="utf-8")
            self.assertIn("stage5-dashboard", repro)
            self.assertIn("stage6-validate", repro)
            self.assertIn("stage7-package", repro)

            checklist = rows_from_csv(checklist_path)
            self.assertTrue(any(row["section"] == "demo_story" for row in checklist))
            self.assertTrue(any(row["section"] == "risk_boundary" for row in checklist))

            stage8_spec = stage8_spec_path.read_text(encoding="utf-8")
            self.assertIn("Stage 8 Frontend Product Specification", stage8_spec)
            self.assertIn("Project Command Center", stage8_spec)
            self.assertIn("Target Evidence Workspace", stage8_spec)
            self.assertIn("Candidate Funnel", stage8_spec)
            self.assertIn("Validation Operations", stage8_spec)


if __name__ == "__main__":
    unittest.main()
