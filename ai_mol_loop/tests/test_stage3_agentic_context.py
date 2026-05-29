import argparse
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


class Stage3AgenticContextTests(unittest.TestCase):
    def build_prior_round_project(self, project: Path) -> None:
        ai_mol_loop.ensure_project_dirs(project)
        config = ai_mol_loop.default_config()
        config["generator_prompt_file"] = str(project / "prompts" / "generator_prompt.md")
        config["target_brief_file"] = str(project / "briefs" / "target_brief.json")
        ai_mol_loop.write_json(project / "config.json", config)
        ai_mol_loop.write_text(
            project / "prompts" / "generator_prompt.md",
            "Base generator prompt for influenza neuraminidase.",
        )
        ai_mol_loop.write_json(
            project / "briefs" / "target_brief.json",
            {
                "target_catalog_id": "influenza_a_h1n1_na",
                "target_name": "Influenza A neuraminidase",
                "compliance_boundary": {
                    "scope": "computational screening only",
                    "no_claims": ["Do not claim clinical efficacy."],
                },
            },
        )
        ai_mol_loop.write_json(
            project / "evidence" / "stage2_closed_loop_assets.json",
            {
                "targets": [
                    {
                        "target_id": "influenza_a_h1n1_na",
                        "evidence_quality": {"readiness": "ready_for_mvp_closed_loop", "score": 0.91},
                        "closed_loop_assets": {
                            "primary_structure": {"pdb_id": "3TI6"},
                            "positive_controls": ["oseltamivir"],
                        },
                    }
                ]
            },
        )
        ai_mol_loop.write_csv(
            project / "filtered" / "round_1_filtered.csv",
            [
                {
                    "id": "r01_bad",
                    "smiles": "C(=O)Cl",
                    "passed_filter": "false",
                    "filter_reason": "structural_risk_flags",
                    "risk_flags": "acid_chloride",
                }
            ],
            ["id", "smiles", "passed_filter", "filter_reason", "risk_flags"],
        )
        ai_mol_loop.write_csv(
            project / "ranked" / "round_1_ranked.csv",
            [
                {
                    "rank": "1",
                    "id": "r01_good",
                    "smiles": "CCOC(=O)c1ccccc1",
                    "total_proxy": "0.81",
                    "docking_proxy": "0.74",
                    "pose_proxy": "1.0",
                    "score_source": "proxy+stage4_rdkit+external_docking",
                    "decision": "advance",
                },
                {
                    "rank": "2",
                    "id": "r01_hold",
                    "smiles": "CCS(=O)(=O)Nc1ccccc1",
                    "total_proxy": "0.52",
                    "docking_proxy": "0.38",
                    "pose_proxy": "0.0",
                    "score_source": "proxy+stage4_rdkit",
                    "decision": "hold",
                },
            ],
            ["rank", "id", "smiles", "total_proxy", "docking_proxy", "pose_proxy", "score_source", "decision"],
        )
        ai_mol_loop.write_json(
            project / "feedback" / "round_1_feedback.json",
            {
                "round": 1,
                "selected_smiles": ["CCOC(=O)c1ccccc1"],
                "recommended_next_steps": ["Bias toward morpholine-bearing analogs."],
            },
        )
        ai_mol_loop.write_csv(
            project / "stage4" / "round_1_docking_scores_template.csv",
            [
                {
                    "id": "r01_good",
                    "smiles": "CCOC(=O)c1ccccc1",
                    "docking_score": "-7.2",
                    "pose_pass": "true",
                    "backend": "vina",
                    "receptor": "3TI6_protein_only_obabel.pdbqt",
                    "notes": "posebusters=passed",
                }
            ],
            ["id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"],
        )
        ai_mol_loop.write_json(
            project / "stage4_5" / "round_1_control_validation.json",
            {
                "best_known_control": {"id": "oseltamivir", "docking_score": -7.6},
                "best_candidate": {"id": "r01_good", "docking_score": -7.2},
                "best_decoy": {"id": "decoy_001", "docking_score": -4.1},
                "candidate_vs_controls": {"best_candidate_delta_kcal_mol": 0.4},
            },
        )
        ai_mol_loop.write_json(
            project / "stage4_6" / "round_1_retrospective_benchmark.json",
            {"metrics": {"roc_auc": 0.88, "best_candidate": {"id": "r01_good", "rank": 2}}},
        )

    def test_openai_prompt_includes_prior_round_feedback_and_validation_evidence(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            self.build_prior_round_project(project)

            args = argparse.Namespace(
                round=2,
                n=12,
                prompt="Generate next round candidates.",
                context_url=[],
                timeout=10,
                max_url_bytes=800000,
                api_key="test-api-key-secret-should-not-appear",
                api_key_file=None,
                api_key_env="OPENAI_API_KEY",
            )

            prompt = ai_mol_loop.build_openai_candidate_prompt(project, args)

            self.assertIn("Prior Round Feedback Context", prompt)
            self.assertIn("round_1", prompt)
            self.assertIn("r01_good", prompt)
            self.assertIn("advance", prompt)
            self.assertIn("r01_hold", prompt)
            self.assertIn("hold", prompt)
            self.assertIn("structural_risk_flags", prompt)
            self.assertIn("acid_chloride", prompt)
            self.assertIn("Bias toward morpholine-bearing analogs", prompt)
            self.assertIn("-7.2", prompt)
            self.assertIn("pose_pass", prompt)
            self.assertIn("oseltamivir", prompt)
            self.assertIn("roc_auc", prompt)
            self.assertIn("source confidence", prompt.lower())
            self.assertIn("Do not claim clinical efficacy", prompt)
            self.assertNotIn("test-api-key-secret-should-not-appear", prompt)

    def test_stage3_screen_saves_openai_prompt_snapshot_without_api_key(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            self.build_prior_round_project(project)
            secret = "test-api-key-secret-should-not-be-written"
            original = ai_mol_loop.openai_generate_candidates

            def fake_openai_generate_candidates(prompt, api_key, model, timeout):
                self.assertEqual(api_key, secret)
                self.assertIn("Prior Round Feedback Context", prompt)
                return [
                    {
                        "id": "llm_001",
                        "smiles": "CCN(CC)C(=O)c1ccccc1",
                        "rationale": "Computational analog for screening.",
                        "expected_interaction": "High-level pocket fit hypothesis.",
                        "design_family": "amide aromatic",
                        "risk_note": "Computational prioritization only.",
                    }
                ]

            ai_mol_loop.openai_generate_candidates = fake_openai_generate_candidates
            try:
                ai_mol_loop.stage3_screen(
                    argparse.Namespace(
                        project=str(project),
                        round=2,
                        n=4,
                        top=1,
                        source_csv=None,
                        source_json=None,
                        source_url=[],
                        context_url=[],
                        use_openai=True,
                        openai_model="gpt-test",
                        api_key=secret,
                        api_key_file=None,
                        api_key_env="OPENAI_API_KEY",
                        prompt="Generate one safer analog.",
                        timeout=10,
                        max_url_bytes=800000,
                        max_heavy_atoms=55,
                        max_molecular_weight=550.0,
                        max_lipinski_violations=1,
                        allow_risk=False,
                        external_scores=None,
                        no_score=True,
                    )
                )
            finally:
                ai_mol_loop.openai_generate_candidates = original

            prompt_path = project / "prompts" / "round_2_openai_candidate_prompt.md"
            assets = json.loads((project / "stage3" / "round_2_stage3_assets.json").read_text(encoding="utf-8"))
            self.assertTrue(prompt_path.exists())
            prompt_text = prompt_path.read_text(encoding="utf-8")
            self.assertIn("Prior Round Feedback Context", prompt_text)
            self.assertIn("round_1", prompt_text)
            self.assertNotIn(secret, prompt_text)
            self.assertEqual(assets["openai"]["prompt_file"], str(prompt_path))
            self.assertTrue(assets["openai"]["context"]["uses_feedback_context"])
            self.assertEqual(assets["openai"]["context"]["prior_rounds"], [1])
            for path in project.rglob("*"):
                if path.is_file():
                    try:
                        content = path.read_text(encoding="utf-8")
                    except UnicodeDecodeError:
                        continue
                    self.assertNotIn(secret, content, str(path))


if __name__ == "__main__":
    unittest.main()
