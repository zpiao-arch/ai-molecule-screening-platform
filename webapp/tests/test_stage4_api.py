import tempfile
import unittest
from pathlib import Path

from fastapi.testclient import TestClient

import server


class Stage4ApiTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(dir=Path.cwd())
        self.projects_root = Path(self.tmp.name) / "projects"
        self.projects_root.mkdir(parents=True, exist_ok=True)
        self.original_projects_root = server.PROJECTS_ROOT
        server.PROJECTS_ROOT = self.projects_root
        self.client = TestClient(server.app)

    def tearDown(self):
        server.PROJECTS_ROOT = self.original_projects_root
        self.tmp.cleanup()

    def create_project(self, name="stage4_demo"):
        response = self.client.post("/api/projects", json={"name": name})
        self.assertEqual(response.status_code, 200, response.text)
        project = self.projects_root / name
        server.write_json(
            project / "briefs" / "target_brief.json",
            {
                "target_catalog_id": "influenza_a_h1n1_na",
                "target_name": "Influenza A(H1N1) neuraminidase",
                "disease_context": "甲流",
                "binding_site": {"reference_ligand": "oseltamivir", "key_residues": ["Arg118", "Glu276"]},
            },
        )
        server.write_csv(
            project / "evidence" / "stage2_target_sources.csv",
            [
                {
                    "target_id": "influenza_a_h1n1_na",
                    "target": "Influenza A(H1N1) neuraminidase",
                    "evidence_score": "0.9967",
                    "readiness": "ready_for_mvp_closed_loop",
                    "primary_pdb": "3TI6",
                    "positive_controls": "oseltamivir; zanamivir; peramivir",
                    "pubmed_articles": "8",
                }
            ],
            ["target_id", "target", "evidence_score", "readiness", "primary_pdb", "positive_controls", "pubmed_articles"],
        )
        server.write_json(
            project / "evidence" / "stage2_closed_loop_assets.json",
            {"targets": [{"target_id": "influenza_a_h1n1_na", "target": "Influenza A(H1N1) neuraminidase"}]},
        )
        server.write_csv(
            project / "candidates" / "round_1_candidates.csv",
            [
                {"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"},
                {"round": 1, "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "parent": "", "source": "unit"},
                {"round": 1, "id": "bad_001", "smiles": "not_a_smiles", "parent": "", "source": "unit"},
            ],
            ["round", "id", "smiles", "parent", "source"],
        )
        return name, project

    def test_stage4_real_generates_full_real_library_assets(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage4/real",
            json={
                "round": 1,
                "target": "influenza_a_h1n1_na",
                "top": 2,
                "decoys": 4,
                "max_conformers": 0,
                "no_sdf": True,
                "render_2d": True,
                "rescore": True,
                "rank_top": 2,
                "feedback_top": 2,
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 4)
        self.assertEqual(payload["command"], "stage4-real")
        self.assertTrue(payload["has_assets"])
        self.assertTrue(payload["has_report"])
        self.assertEqual(payload["assets"]["stage"], 4)
        self.assertEqual(payload["assets"]["target_id"], "influenza_a_h1n1_na")
        self.assertEqual(payload["assets"]["valid_count"], 2)
        self.assertEqual(payload["assets"]["invalid_count"], 1)
        self.assertTrue(payload["descriptors"])
        self.assertTrue(payload["similarity"])
        self.assertTrue(payload["diverse_selection"])
        self.assertTrue(payload["benchmark_panel"])
        self.assertTrue(payload["decoys"])
        self.assertIn("panel_counts", payload["validation_metrics"])
        self.assertIn("status", payload["docking_plan"])
        self.assertIn("preparation_status", payload["receptor_package"])
        self.assertIn("Stage 4 Real Library Validation", payload["report"])
        self.assertTrue((project / "stage4" / "round_1_real_descriptors.csv").exists())
        self.assertTrue((project / "stage4" / "round_1_stage4_assets.json").exists())
        self.assertTrue((project / "reports" / "stage4_round_1_report.md").exists())
        self.assertTrue((project / "scores" / "round_1_scores.csv").exists())
        self.assertTrue((project / "ranked" / "round_1_ranked.csv").exists())
        self.assertTrue((project / "feedback" / "round_1_feedback.json").exists())

    def test_stage4_status_reads_existing_assets_tables_and_image_urls(self):
        name, _project = self.create_project()
        self.client.post(
            f"/api/projects/{name}/stage4/real",
            json={"round": 1, "target": "influenza_a_h1n1_na", "top": 2, "decoys": 2, "max_conformers": 0},
        )

        response = self.client.get(f"/api/projects/{name}/stage4?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 4)
        self.assertTrue(payload["has_assets"])
        self.assertTrue(payload["has_report"])
        self.assertEqual(payload["assets"]["round"], 1)
        self.assertTrue(payload["descriptors"])
        self.assertTrue(payload["benchmark_panel"])
        self.assertTrue(payload["molecule_images"])
        self.assertTrue(payload["molecule_images"][0]["url"].startswith(f"/api/projects/{name}/stage4/images/1/"))

    def test_stage4_preflight_reports_missing_candidate_input_when_round_csv_is_absent(self):
        name, project = self.create_project()
        candidate_path = project / "candidates" / "round_1_candidates.csv"
        candidate_path.unlink()

        response = self.client.get(f"/api/projects/{name}/stage4?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        step = next(item for item in payload["stage4_preflight"] if item["step_id"] == "candidate_input")
        self.assertEqual(step["status"], "missing")
        self.assertIn("candidates/round_1_candidates.csv", step["evidence"])
        self.assertFalse((project / "candidates" / "round_1_candidates.csv").exists())

    def test_stage4_real_accepts_docking_box_and_persists_it_to_project_config(self):
        name, project = self.create_project()

        response = self.client.post(
            f"/api/projects/{name}/stage4/real",
            json={
                "round": 1,
                "target": "influenza_a_h1n1_na",
                "top": 1,
                "decoys": 0,
                "max_conformers": 0,
                "no_sdf": True,
                "render_2d": False,
                "pocket_center": "1.0, 2.0, 3.0",
                "pocket_size": "18, 19, 20",
                "pocket_source": "manual_unit_test",
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        config = server.read_json(project / "config.json")
        pocket = config["target"]["pocket"]
        self.assertEqual(pocket["center"], [1.0, 2.0, 3.0])
        self.assertEqual(pocket["size"], [18.0, 19.0, 20.0])
        self.assertEqual(pocket["source"], "manual_unit_test")
        plan = response.json()["docking_plan"]
        self.assertEqual(plan["docking_box"]["status"], "ready")
        self.assertEqual(plan["docking_box"]["center"], [1.0, 2.0, 3.0])
        self.assertEqual(plan["docking_box"]["size"], [18.0, 19.0, 20.0])

    def test_stage4_real_auto_uses_project_receptor_files_for_known_pdb_id(self):
        name, project = self.create_project()
        receptor_dir = project / "stage4" / "receptors"
        receptor_dir.mkdir(parents=True, exist_ok=True)
        receptor_pdb = receptor_dir / "3TI6.pdb"
        receptor_pdbqt = receptor_dir / "3TI6_protein_only_obabel.pdbqt"
        receptor_pdb.write_text(
            "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\nEND\n",
            encoding="utf-8",
        )
        receptor_pdbqt.write_text("REMARK receptor pdbqt\n", encoding="utf-8")

        response = self.client.post(
            f"/api/projects/{name}/stage4/real",
            json={
                "round": 1,
                "target": "influenza_a_h1n1_na",
                "pdb_id": "3TI6",
                "top": 1,
                "decoys": 0,
                "max_conformers": 0,
                "no_sdf": True,
                "render_2d": False,
                "run_docking": False,
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        rec = payload["receptor_package"]
        self.assertEqual(rec["local_receptor_pdb"], str(receptor_pdb.resolve()))
        self.assertEqual(rec["local_receptor_pdbqt"], str(receptor_pdbqt.resolve()))
        receptor_gate = next(item for item in payload["stage4_preflight"] if item["step_id"] == "receptor")
        self.assertEqual(receptor_gate["status"], "ready")
        viz = payload["visualization_assets"]
        self.assertEqual(viz["receptor"]["filename"], "3TI6.pdb")
        self.assertNotEqual(viz["status"], "missing_receptor")

    def test_stage4_capabilities_reports_installable_real_docking_stack(self):
        response = self.client.get("/api/stage4/capabilities")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 4)
        self.assertIn("modules", payload)
        self.assertIn("executables", payload)
        self.assertIn("docking_backend", payload)
        for module in ["rdkit", "meeko", "posebusters", "vina", "openbabel"]:
            self.assertIn(module, payload["modules"])
        for executable in ["vina", "gnina", "obabel", "bust"]:
            self.assertIn(executable, payload["executables"])
        self.assertIn("install_options", payload)
        self.assertIn("pip_user", payload["install_options"])
        self.assertIn("binary_tools", payload["install_options"])
        self.assertIn("docker_tools", payload["install_options"])
        self.assertIn("boundary", payload)

    def test_stage4_status_recovers_real_docking_evidence_from_stage45_assets(self):
        name, project = self.create_project()
        receptor_dir = project / "stage4" / "receptors"
        receptor_dir.mkdir(parents=True, exist_ok=True)
        receptor_pdb = receptor_dir / "3TI6.pdb"
        receptor_pdbqt = receptor_dir / "3TI6_protein_only_obabel.pdbqt"
        receptor_pdb.write_text(
            "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\nEND\n",
            encoding="utf-8",
        )
        receptor_pdbqt.write_text("REMARK receptor pdbqt\n", encoding="utf-8")
        run_dir = project / "stage4_5" / "control_docking_runs"
        run_dir.mkdir(parents=True, exist_ok=True)
        pose = run_dir / "mol_001_pose.sdf"
        log = run_dir / "mol_001.log"
        posebusters = run_dir / "mol_001_posebusters.csv"
        pose.write_text("mol_001\n  unit\n\nM  END\n$$$$\n", encoding="utf-8")
        log.write_text("REMARK vina log\n", encoding="utf-8")
        posebusters.write_text("check,status\ngeometry,passed\n", encoding="utf-8")
        server.write_json(
            project / "stage4" / "round_1_docking_plan.json",
            {
                "status": "skipped",
                "selected_backend": "vina",
                "receptor_status": "pdb_id_known_receptor_file_not_prepared",
                "docking_box": {"status": "missing", "center": [], "size": [], "source": "co_crystal"},
                "expected_scores_csv": str(project / "stage4" / "round_1_docking_scores_template.csv"),
            },
        )
        server.write_json(
            project / "stage4" / "round_1_receptor_package.json",
            {
                "target_id": "influenza_a_h1n1_na",
                "pdb_id": "3TI6",
                "local_receptor_pdb": "",
                "local_receptor_pdbqt": "",
                "preparation_status": "pdb_id_known_receptor_file_not_prepared",
                "binding_site": {"center": [], "size": []},
            },
        )
        server.write_json(
            project / "stage4" / "round_1_validation_metrics.json",
            {
                "panel_counts": {"candidate": 3, "positive_control": 0, "decoy": 0},
                "readiness": {"rdkit_descriptors": "ready", "control_panel": "missing", "decoy_panel": "missing", "docking": "skipped"},
            },
        )
        server.write_csv(
            project / "stage4_5" / "round_1_control_docking_inputs.csv",
            [
                {
                    "panel_type": "candidate",
                    "id": "mol_001",
                    "canonical_smiles": "CCOC(=O)c1ccccc1",
                    "ligand_pdbqt": str(project / "stage4_5" / "round_1_prepared_ligands" / "mol_001.pdbqt"),
                    "receptor_pdb": str(receptor_pdb),
                    "receptor_pdbqt": str(receptor_pdbqt),
                    "status": "ready_for_docking",
                    "note": "Prepared ligand PDBQT with OpenBabel.",
                }
            ],
            ["panel_type", "id", "canonical_smiles", "ligand_pdbqt", "receptor_pdb", "receptor_pdbqt", "status", "note"],
        )
        server.write_csv(
            project / "stage4_5" / "round_1_control_docking_scores.csv",
            [
                {
                    "panel_type": "candidate",
                    "id": "mol_001",
                    "smiles": "CCOC(=O)c1ccccc1",
                    "docking_score": "-6.54",
                    "pose_pass": "true",
                    "backend": "vina",
                    "receptor": str(receptor_pdbqt),
                    "notes": f"returncode=0;log={log};posebusters=passed;posebusters_report={posebusters}",
                }
            ],
            ["panel_type", "id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"],
        )
        server.write_json(
            project / "stage4_5" / "round_1_control_validation.json",
            {
                "docking": {
                    "status": "completed",
                    "backend": "vina",
                    "scored_count": 1,
                    "pose_pass_count": 1,
                    "docking_box": {"status": "ready", "center": [-28.914, 14.334, 20.794], "size": [23.585, 20.45, 24.18], "source": "binding_site"},
                }
            },
        )

        response = self.client.get(f"/api/projects/{name}/stage4?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["docking_plan"]["status"], "completed")
        self.assertEqual(payload["docking_plan"]["docking_box"]["status"], "ready")
        self.assertEqual(payload["receptor_package"]["local_receptor_pdb"], str(receptor_pdb))
        self.assertEqual(payload["receptor_package"]["local_receptor_pdbqt"], str(receptor_pdbqt))
        self.assertEqual(payload["stage4_readiness_guide"]["docking"], "completed")
        preflight = {row["step_id"]: row["status"] for row in payload["stage4_preflight"]}
        self.assertEqual(preflight["receptor"], "ready")
        self.assertEqual(preflight["docking_box"], "ready")
        self.assertEqual(preflight["ligand_preparation"], "ready")
        self.assertEqual(preflight["docking_scores"], "ready")
        self.assertEqual(payload["visualization_assets"]["status"], "ready")
        self.assertEqual(payload["visualization_assets"]["receptor"]["filename"], "3TI6.pdb")
        self.assertEqual(payload["visualization_assets"]["top_pose"]["filename"], "mol_001_pose.sdf")
        self.assertEqual(payload["docking_results"][0]["evidence_source"], "stage4_5")

    def test_stage4_repair_persists_recovered_docking_metadata(self):
        name, project = self.create_project("stage4_repair_demo")
        receptor_dir = project / "stage4" / "receptors"
        receptor_dir.mkdir(parents=True, exist_ok=True)
        receptor_pdb = receptor_dir / "3TI6.pdb"
        receptor_pdbqt = receptor_dir / "3TI6_protein_only_obabel.pdbqt"
        receptor_pdb.write_text("ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\nEND\n", encoding="utf-8")
        receptor_pdbqt.write_text("REMARK receptor pdbqt\n", encoding="utf-8")
        run_dir = project / "stage4_5" / "control_docking_runs"
        run_dir.mkdir(parents=True, exist_ok=True)
        pose = run_dir / "mol_001_pose.sdf"
        log = run_dir / "mol_001.log"
        posebusters = run_dir / "mol_001_posebusters.csv"
        pose.write_text("mol_001\n  unit\n\nM  END\n$$$$\n", encoding="utf-8")
        log.write_text("REMARK vina log\n", encoding="utf-8")
        posebusters.write_text("check,status\ngeometry,passed\n", encoding="utf-8")
        server.write_json(project / "stage4" / "round_1_docking_plan.json", {"status": "skipped", "selected_backend": "vina", "docking_box": {"status": "missing"}})
        server.write_json(project / "stage4" / "round_1_receptor_package.json", {"target_id": "influenza_a_h1n1_na", "pdb_id": "3TI6"})
        server.write_json(project / "stage4" / "round_1_validation_metrics.json", {"readiness": {"docking": "skipped"}})
        server.write_json(project / "stage4" / "round_1_stage4_assets.json", {"stage": 4, "round": 1, "docking_plan": {"status": "skipped"}})
        server.write_csv(project / "stage4" / "round_1_real_descriptors.csv", [{"id": "mol_001"}], ["id"])
        server.write_csv(
            project / "stage4_5" / "round_1_control_docking_inputs.csv",
            [{"panel_type": "candidate", "id": "mol_001", "receptor_pdb": str(receptor_pdb), "receptor_pdbqt": str(receptor_pdbqt), "status": "ready_for_docking"}],
            ["panel_type", "id", "receptor_pdb", "receptor_pdbqt", "status"],
        )
        server.write_csv(
            project / "stage4_5" / "round_1_control_docking_scores.csv",
            [{"panel_type": "candidate", "id": "mol_001", "smiles": "CCO", "docking_score": "-6.1", "pose_pass": "true", "backend": "vina", "receptor": str(receptor_pdbqt), "notes": f"returncode=0;log={log};posebusters_report={posebusters}"}],
            ["panel_type", "id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"],
        )
        server.write_json(project / "stage4_5" / "round_1_control_validation.json", {"docking": {"status": "completed", "backend": "vina", "docking_box": {"status": "ready", "center": [1, 2, 3], "size": [20, 20, 20], "source": "unit"}}})

        response = self.client.post(f"/api/projects/{name}/stage4/repair", json={"round": 1})

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["command"], "stage4-repair")
        self.assertTrue(payload["repaired"])
        self.assertEqual(payload["docking_plan"]["status"], "completed")
        self.assertEqual(server.read_json(project / "stage4" / "round_1_docking_plan.json")["status"], "completed")
        self.assertEqual(server.read_json(project / "stage4" / "round_1_validation_metrics.json")["readiness"]["docking"], "completed")
        self.assertEqual(server.read_json(project / "stage4" / "round_1_stage4_assets.json")["docking_results_count"], 1)

        doctor_response = self.client.get(f"/api/projects/{name}/doctor?round=1")
        self.assertEqual(doctor_response.status_code, 200, doctor_response.text)
        doctor = doctor_response.json()
        self.assertEqual(doctor["status"], "ok")
        self.assertEqual(doctor["recovered_docking_evidence"]["status"], "consistent")
        self.assertNotIn("recovered_docking", [issue.get("step_id") for issue in doctor["issues"]])

    def test_stage4_presets_returns_target_and_pocket_defaults(self):
        response = self.client.get("/api/stage4/presets")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 4)
        self.assertIn("presets", payload)
        presets = {item["target_id"]: item for item in payload["presets"]}
        self.assertIn("influenza_a_h1n1_na", presets)
        na = presets["influenza_a_h1n1_na"]
        self.assertEqual(na["recommended_pdb"], "3TI6")
        self.assertIn("oseltamivir", na["positive_controls"])
        self.assertIn("pocket", na)
        self.assertEqual(na["pocket"]["status"], "needs_curated_coordinates")
        smoke = {item["target_id"]: item for item in payload["presets"]}["vina_official_1iep_smoke"]
        self.assertEqual(smoke["pocket"]["center"], [15.19, 53.903, 16.917])
        self.assertEqual(smoke["pocket"]["size"], [20.0, 20.0, 20.0])
        self.assertTrue(smoke["local_receptor_pdbqt"].endswith("1iep_receptor.pdbqt"))

    def test_stage4_payload_includes_docking_results_and_readiness_guide(self):
        name, project = self.create_project()
        stage4_dir = project / "stage4"
        docking_dir = stage4_dir / "docking_runs"
        receptor_dir = stage4_dir / "receptors"
        docking_dir.mkdir(parents=True)
        receptor_dir.mkdir(parents=True)
        receptor_pdb = receptor_dir / "3TI6_protein_only.pdb"
        receptor_pdb.write_text("ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\nEND\n", encoding="utf-8")
        server.write_json(
            stage4_dir / "round_1_docking_plan.json",
            {
                "status": "completed",
                "selected_backend": "vina",
                "docking_box": {"status": "ready", "center": [1, 2, 3], "size": [20, 20, 20]},
                "expected_scores_csv": str(stage4_dir / "round_1_docking_scores_template.csv"),
            },
        )
        server.write_json(
            stage4_dir / "round_1_receptor_package.json",
            {
                "target_id": "influenza_a_h1n1_na",
                "pdb_id": "3TI6",
                "local_receptor_pdb": str(receptor_pdb),
                "local_receptor_pdbqt": "/tmp/receptor.pdbqt",
            },
        )
        server.write_json(
            stage4_dir / "round_1_stage4_assets.json",
            {"stage": 4, "round": 1, "target_id": "influenza_a_h1n1_na", "valid_count": 1},
        )
        server.write_csv(
            stage4_dir / "round_1_docking_scores_template.csv",
            [
                {
                    "id": "mol_001",
                    "smiles": "CCO",
                    "docking_score": "-7.8",
                    "pose_pass": "true",
                    "backend": "vina",
                    "receptor": "/tmp/receptor.pdbqt",
                    "notes": f"returncode=0;log={docking_dir / 'mol_001.log'};posebusters=passed;posebusters_report={docking_dir / 'mol_001_posebusters.csv'}",
                }
            ],
            ["id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"],
        )
        (docking_dir / "mol_001.log").write_text("mode | affinity\n1 -7.8\n", encoding="utf-8")
        (docking_dir / "mol_001_pose.pdbqt").write_text("REMARK pose\n", encoding="utf-8")
        (docking_dir / "mol_001_pose.sdf").write_text("pose\n$$$$\n", encoding="utf-8")
        (docking_dir / "mol_001_posebusters.csv").write_text("mol_pred,all_atoms_connected\nmol_001,true\n", encoding="utf-8")

        response = self.client.get(f"/api/projects/{name}/stage4?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertIn("docking_results", payload)
        self.assertEqual(payload["docking_results"][0]["id"], "mol_001")
        self.assertEqual(payload["docking_results"][0]["interpretation"]["score_band"], "strong")
        self.assertTrue(payload["docking_results"][0]["pose_file"].endswith("mol_001_pose.pdbqt"))
        self.assertTrue(payload["docking_results"][0]["log_file"].endswith("mol_001.log"))
        self.assertIn("stage4_readiness_guide", payload)
        self.assertEqual(payload["stage4_readiness_guide"]["docking"], "completed")
        self.assertTrue(payload["stage4_readiness_guide"]["next_actions"])
        self.assertIn("stage4_preflight", payload)
        preflight = {item["step_id"]: item for item in payload["stage4_preflight"]}
        self.assertEqual(preflight["docking_scores"]["status"], "ready")
        self.assertEqual(preflight["docking_backend"]["status"], "ready")
        self.assertIn("visualization_assets", payload)
        viz = payload["visualization_assets"]
        self.assertEqual(viz["viewer"], "3Dmol.js")
        self.assertEqual(viz["status"], "ready")
        self.assertEqual(viz["receptor"]["format"], "pdb")
        self.assertIn(f"/api/projects/{name}/artifact?path=", viz["receptor"]["url"])
        self.assertEqual(viz["top_pose"]["id"], "mol_001")
        self.assertEqual(viz["top_pose"]["format"], "sdf")
        self.assertIn(f"/api/projects/{name}/artifact?path=", viz["top_pose"]["url"])
        self.assertEqual(viz["docking_box"]["center"], [1, 2, 3])

    def test_stage4_real_imports_external_scores_into_stage4_results(self):
        name, project = self.create_project()
        external_scores = project / "external_docking_scores.csv"
        server.write_csv(
            external_scores,
            [
                {
                    "id": "mol_001",
                    "smiles": "CCOC(=O)c1ccccc1",
                    "docking_score": "-8.2",
                    "pose_pass": "true",
                    "backend": "external_vina",
                    "receptor": "/tmp/curated_receptor.pdbqt",
                    "notes": "curated external docking run",
                }
            ],
            ["id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"],
        )

        response = self.client.post(
            f"/api/projects/{name}/stage4/real",
            json={
                "round": 1,
                "target": "influenza_a_h1n1_na",
                "top": 2,
                "decoys": 2,
                "max_conformers": 0,
                "no_sdf": True,
                "render_2d": False,
                "rescore": True,
                "external_scores": str(external_scores),
            },
        )

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["docking_plan"]["status"], "imported_external_scores")
        self.assertEqual(payload["docking_results"][0]["id"], "mol_001")
        self.assertEqual(payload["docking_results"][0]["score_band"], "strong")
        self.assertEqual(payload["docking_results"][0]["backend"], "external_vina")
        imported = server.csv_to_dicts(project / "stage4" / "round_1_docking_scores_template.csv")
        self.assertEqual(imported[0]["docking_score"], "-8.2")
        ranked = server.csv_to_dicts(project / "ranked" / "round_1_ranked.csv")
        self.assertIn("external_docking", ranked[0]["score_source"])

    def test_stage4_smoke_test_route_returns_structured_toolchain_status(self):
        response = self.client.post("/api/stage4/smoke-test", json={"timeout": 1})

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 4)
        self.assertIn("status", payload)
        self.assertIn("message", payload)
        self.assertIn("example", payload)

    def test_stage45_status_reads_control_calibration_outputs(self):
        name, project = self.create_project()
        stage45 = project / "stage4_5"
        stage45.mkdir(parents=True, exist_ok=True)
        server.write_csv(
            stage45 / "round_1_control_docking_scores.csv",
            [
                {
                    "panel_type": "positive_control",
                    "id": "oseltamivir",
                    "smiles": "CCOC(=O)C1=C(C)OC(CC)C(C1NC(C)=O)N",
                    "docking_score": "-8.4",
                    "pose_pass": "true",
                    "backend": "vina",
                    "receptor": "/tmp/receptor.pdbqt",
                    "relative_to_best_control": "0.0",
                    "score_band": "strong",
                    "notes": "unit",
                },
                {
                    "panel_type": "candidate",
                    "id": "mol_001",
                    "smiles": "CCOC(=O)c1ccccc1",
                    "docking_score": "-6.7",
                    "pose_pass": "true",
                    "backend": "vina",
                    "receptor": "/tmp/receptor.pdbqt",
                    "relative_to_best_control": "1.7",
                    "score_band": "moderate",
                    "notes": "unit",
                },
            ],
            [
                "panel_type",
                "id",
                "smiles",
                "docking_score",
                "pose_pass",
                "backend",
                "receptor",
                "relative_to_best_control",
                "score_band",
                "notes",
            ],
        )
        server.write_json(
            stage45 / "round_1_control_validation.json",
            {
                "schema_version": "0.1",
                "stage": 4.5,
                "round": 1,
                "target_id": "influenza_a_h1n1_na",
                "counts": {"candidate": 1, "known_control": 1, "decoy": 0},
                "docking": {"status": "completed", "backend": "vina"},
                "best_known_control": {"id": "oseltamivir", "docking_score": -8.4},
                "best_candidate": {"id": "mol_001", "docking_score": -6.7},
                "candidate_vs_controls": {"best_candidate_delta_kcal_mol": 1.7},
                "redocking": {"status": "reference_pose_exported", "rmsd_status": "unavailable"},
                "boundary": ["Controls calibrate computation; this does not prove efficacy."],
            },
        )
        report_path = project / "reports" / "stage4_5_round_1_control_validation.md"
        report_path.write_text("# Stage 4.5 Control Calibration\n", encoding="utf-8")

        response = self.client.get(f"/api/projects/{name}/stage45?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 4.5)
        self.assertTrue(payload["has_validation"])
        self.assertEqual(payload["validation"]["best_known_control"]["id"], "oseltamivir")
        self.assertEqual(payload["scores"][0]["score_band"], "strong")
        self.assertIn("Stage 4.5 Control Calibration", payload["report"])
        self.assertTrue(payload["files"]["validation"].endswith("round_1_control_validation.json"))

    def test_stage45_status_ignores_placeholder_vina_scores_when_current_backend_invalid(self):
        name, project = self.create_project("stage45_invalid_vina_demo")
        stage4 = project / "stage4"
        stage45 = project / "stage4_5"
        stage4.mkdir(parents=True, exist_ok=True)
        stage45.mkdir(parents=True, exist_ok=True)
        server.write_json(
            stage4 / "round_1_docking_plan.json",
            {
                "status": "not_available",
                "selected_backend": "",
                "tool_status": {
                    "vina": {
                        "status": "invalid",
                        "validation": "placeholder_no_pose_output",
                    }
                },
                "docking_box": {"status": "ready", "center": [0, 0, 0], "size": [20, 20, 20]},
            },
        )
        server.write_csv(
            stage45 / "round_1_control_docking_scores.csv",
            [
                {
                    "panel_type": "positive_control",
                    "id": "oseltamivir",
                    "smiles": "CCO",
                    "docking_score": "-7.4",
                    "pose_pass": "",
                    "backend": "vina",
                    "receptor": "/tmp/receptor.pdbqt",
                    "relative_to_best_control": "0.0",
                    "score_band": "moderate",
                    "notes": "returncode=0;log=/tmp/oseltamivir.log;posebusters=missing_pose",
                }
            ],
            [
                "panel_type",
                "id",
                "smiles",
                "docking_score",
                "pose_pass",
                "backend",
                "receptor",
                "relative_to_best_control",
                "score_band",
                "notes",
            ],
        )
        server.write_json(
            stage45 / "round_1_control_validation.json",
            {
                "docking": {"status": "completed", "backend": "vina"},
                "best_known_control": {"id": "oseltamivir", "docking_score": -7.4},
            },
        )

        response = self.client.get(f"/api/projects/{name}/stage45?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertTrue(payload["has_scores"])
        self.assertEqual(payload["scores"], [])
        self.assertEqual(payload["score_filter"]["status"], "ignored_invalid_backend")
        self.assertEqual(payload["score_filter"]["raw_scores"], 1)
        self.assertIn("placeholder_no_pose_output", payload["score_filter"]["reason"])
        self.assertEqual(payload["validation"]["docking"]["status"], "not_available")

    def test_stage46_status_reads_retrospective_benchmark_outputs(self):
        name, project = self.create_project()
        stage46 = project / "stage4_6"
        stage46.mkdir(parents=True, exist_ok=True)
        server.write_csv(
            stage46 / "round_1_retrospective_ranking.csv",
            [
                {"rank": "1", "id": "oseltamivir", "panel_type": "positive_control", "label": "positive", "docking_score": "-8.4"},
                {"rank": "2", "id": "decoy_001", "panel_type": "decoy", "label": "negative", "docking_score": "-4.2"},
            ],
            ["rank", "id", "panel_type", "label", "docking_score"],
        )
        server.write_json(
            stage46 / "round_1_retrospective_benchmark.json",
            {
                "stage": 4.6,
                "round": 1,
                "counts": {"positives": 1, "negatives": 1, "candidates": 0},
                "metrics": {
                    "roc_auc": 1.0,
                    "top_k": {"top_1": {"control_hit_count": 1, "control_hit_rate": 1.0, "enrichment_factor": 2.0}},
                },
                "boundary": ["Retrospective benchmark is not wet-lab validation."],
            },
        )
        report_path = project / "reports" / "stage4_6_round_1_retrospective_benchmark.md"
        report_path.write_text("# Stage 4.6 Retrospective Benchmark\n", encoding="utf-8")

        response = self.client.get(f"/api/projects/{name}/stage46?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertEqual(payload["stage"], 4.6)
        self.assertTrue(payload["has_benchmark"])
        self.assertEqual(payload["benchmark"]["metrics"]["roc_auc"], 1.0)
        self.assertEqual(payload["ranking"][0]["label"], "positive")
        self.assertIn("Stage 4.6 Retrospective Benchmark", payload["report"])
        self.assertTrue(payload["files"]["benchmark"].endswith("round_1_retrospective_benchmark.json"))

    def test_stage46_status_ignores_benchmark_when_current_backend_invalid(self):
        name, project = self.create_project("stage46_invalid_vina_demo")
        stage4 = project / "stage4"
        stage46 = project / "stage4_6"
        stage4.mkdir(parents=True, exist_ok=True)
        stage46.mkdir(parents=True, exist_ok=True)
        server.write_json(
            stage4 / "round_1_docking_plan.json",
            {
                "status": "not_available",
                "tool_status": {
                    "vina": {
                        "status": "invalid",
                        "validation": "placeholder_no_pose_output",
                    }
                },
            },
        )
        server.write_csv(
            stage46 / "round_1_retrospective_ranking.csv",
            [
                {"rank": "1", "id": "decoy_001", "panel_type": "decoy", "label": "negative", "docking_score": "-7.4", "backend": "vina"},
                {"rank": "2", "id": "oseltamivir", "panel_type": "positive_control", "label": "positive", "docking_score": "-7.4", "backend": "vina"},
            ],
            ["rank", "id", "panel_type", "label", "docking_score", "backend"],
        )
        server.write_json(
            stage46 / "round_1_retrospective_benchmark.json",
            {
                "stage": 4.6,
                "round": 1,
                "metrics": {"roc_auc": 0.5},
                "counts": {"positives": 1, "negatives": 1},
            },
        )

        response = self.client.get(f"/api/projects/{name}/stage46?round=1")

        self.assertEqual(response.status_code, 200, response.text)
        payload = response.json()
        self.assertTrue(payload["has_benchmark"])
        self.assertEqual(payload["ranking"], [])
        self.assertEqual(payload["benchmark"], {})
        self.assertEqual(payload["benchmark_filter"]["status"], "ignored_invalid_backend")
        self.assertEqual(payload["benchmark_filter"]["raw_ranking"], 2)
        self.assertIn("placeholder_no_pose_output", payload["benchmark_filter"]["reason"])


if __name__ == "__main__":
    unittest.main()
