import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "webapp" / "server.py"
SPEC = importlib.util.spec_from_file_location("webapp_server_module", MODULE_PATH)
assert SPEC and SPEC.loader
import sys
sys.modules.pop("ai_mol_loop", None)
server = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(server)


def write_csv(path: Path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


class Stage4WebRecoveryTests(unittest.TestCase):
    def test_existing_receptor_metadata_is_consistent_when_stage45_scores_are_ignored(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            stage4 = project / "stage4"
            stage45 = project / "stage4_5"
            receptor_dir = stage4 / "receptors"
            receptor_dir.mkdir(parents=True, exist_ok=True)
            receptor_pdb = receptor_dir / "3TI6.pdb"
            receptor_pdbqt = receptor_dir / "3TI6_protein_only_obabel.pdbqt"
            receptor_pdb.write_text("ATOM      1  N   SER A  82     -47.086  46.266  27.606  1.00 28.97           N\nEND\n", encoding="utf-8")
            receptor_pdbqt.write_text("REMARK receptor pdbqt\nATOM      1  N   SER A  82     -47.086  46.266  27.606  0.00  0.00    +0.000 NA\n", encoding="utf-8")
            server.ensure_project_dirs(project)
            server.write_json(
                stage4 / "round_1_receptor_package.json",
                {
                    "local_receptor_pdb": str(receptor_pdb),
                    "local_receptor_pdbqt": str(receptor_pdbqt),
                    "prepared_receptor_pdb": str(receptor_pdb),
                    "preparation_status": "receptor_pdbqt_available",
                    "binding_site": {"center": [1, 2, 3], "size": [20, 20, 20]},
                },
            )
            server.write_json(
                stage4 / "round_1_docking_plan.json",
                {
                    "status": "not_available",
                    "selected_backend": "",
                    "requested_backend": "auto",
                    "tool_status": {
                        "vina": {
                            "status": "invalid",
                            "validation": "placeholder_no_pose_output",
                        }
                    },
                    "docking_box": {"status": "ready", "center": [1, 2, 3], "size": [20, 20, 20]},
                    "expected_scores_csv": str(stage4 / "round_1_docking_scores_template.csv"),
                },
            )
            server.write_json(stage4 / "round_1_validation_metrics.json", {"readiness": {"docking": "not_available"}})
            write_csv(stage4 / "round_1_real_descriptors.csv", [{"id": "mol_001"}], ["id"])
            write_csv(
                stage45 / "round_1_control_docking_scores.csv",
                [
                    {
                        "panel_type": "candidate",
                        "id": "mol_001",
                        "smiles": "CCO",
                        "docking_score": "-7.4",
                        "pose_pass": "",
                        "backend": "vina",
                        "receptor": str(receptor_pdbqt),
                        "notes": "returncode=0;posebusters=missing_pose",
                    }
                ],
                ["panel_type", "id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"],
            )

            payload = server.stage4_payload("project", project, 1)

            recovered = payload["recovered_docking_evidence"]
            self.assertEqual(recovered["status"], "none")
            self.assertFalse(recovered["metadata_changed"])
            self.assertEqual(payload["visualization_assets"]["status"], "missing_pose")
            self.assertEqual(payload["visualization_assets"]["receptor"]["filename"], "3TI6.pdb")

    def test_missing_pose_visualization_exposes_reference_ligand_when_available(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            stage4 = project / "stage4"
            receptor_dir = stage4 / "receptors"
            receptor_dir.mkdir(parents=True, exist_ok=True)
            receptor_pdb = receptor_dir / "3TI6.pdb"
            receptor_pdb.write_text(
                "\n".join(
                    [
                        "ATOM      1  N   SER A  82     -47.086  46.266  27.606  1.00 28.97           N",
                        "HETATM    2  C1  G39 A 801     -28.899  12.358  24.407  1.00  9.56           C",
                        "HETATM    3  O1A G39 A 801     -29.755  12.004  25.199  1.00 10.07           O",
                        "HETATM    4  O1B G39 A 801     -27.565  12.195  24.683  1.00 10.25           O",
                        "HETATM    5  C2  G39 A 801     -29.231  13.908  24.492  1.00  9.76           C",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            server.ensure_project_dirs(project)
            server.write_json(
                stage4 / "round_1_receptor_package.json",
                {
                    "local_receptor_pdb": str(receptor_pdb),
                    "binding_site": {
                        "reference_ligand": "oseltamivir",
                        "detected_ligand": {"resname": "G39", "chain": "A", "resseq": 801, "icode": ""},
                    },
                },
            )
            server.write_json(
                stage4 / "round_1_docking_plan.json",
                {"status": "not_available", "docking_box": {"status": "ready", "center": [-28.9, 12.3, 24.8], "size": [20, 20, 20]}},
            )

            payload = server.stage4_payload("project", project, 1)

            viz = payload["visualization_assets"]
            self.assertEqual(viz["status"], "missing_pose")
            self.assertEqual(viz["reference_ligand"]["filename"], "reference_ligand_G39_A_801.pdb")
            self.assertEqual(viz["reference_ligand"]["role"], "co_crystal_reference_ligand")
            self.assertTrue(Path(viz["reference_ligand"]["path"]).exists())
            self.assertFalse(viz["poses"])

    def test_missing_pose_visualization_detects_reference_ligand_from_receptor_when_metadata_is_sparse(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            stage4 = project / "stage4"
            receptor_dir = stage4 / "receptors"
            receptor_dir.mkdir(parents=True, exist_ok=True)
            receptor_pdb = receptor_dir / "3TI6.pdb"
            receptor_pdb.write_text(
                "\n".join(
                    [
                        "ATOM      1  N   SER A  82     -47.086  46.266  27.606  1.00 28.97           N",
                        "HETATM    2  C1  G39 A 801     -28.899  12.358  24.407  1.00  9.56           C",
                        "HETATM    3  O1A G39 A 801     -29.755  12.004  25.199  1.00 10.07           O",
                        "HETATM    4  O1B G39 A 801     -27.565  12.195  24.683  1.00 10.25           O",
                        "HETATM    5  C2  G39 A 801     -29.231  13.908  24.492  1.00  9.76           C",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            server.ensure_project_dirs(project)
            server.write_json(
                stage4 / "round_1_receptor_package.json",
                {
                    "local_receptor_pdb": str(receptor_pdb),
                    "binding_site": {"reference_ligand": "oseltamivir"},
                },
            )
            server.write_json(
                stage4 / "round_1_docking_plan.json",
                {"status": "not_available", "docking_box": {"status": "ready", "center": [-28.9, 12.3, 24.8], "size": [20, 20, 20]}},
            )

            payload = server.stage4_payload("project", project, 1)

            ref = payload["visualization_assets"]["reference_ligand"]
            self.assertEqual(ref["resname"], "G39")
            self.assertEqual(ref["filename"], "reference_ligand_G39_A_801.pdb")

    def test_invalid_current_vina_blocks_stage45_recovered_completed_status(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            stage4 = project / "stage4"
            stage45 = project / "stage4_5"
            receptor = stage4 / "receptors" / "receptor.pdbqt"
            receptor.parent.mkdir(parents=True, exist_ok=True)
            receptor.write_text("REMARK RECEPTOR\nEND\n", encoding="utf-8")
            server.ensure_project_dirs(project)
            server.write_json(
                stage4 / "round_1_receptor_package.json",
                {"local_receptor_pdbqt": str(receptor), "preparation_status": "receptor_pdb_available"},
            )
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
                    "expected_scores_csv": str(stage4 / "round_1_docking_scores_template.csv"),
                },
            )
            server.write_json(stage4 / "round_1_validation_metrics.json", {"readiness": {"docking": "not_available"}})
            write_csv(
                stage4 / "round_1_real_descriptors.csv",
                [{"id": "mol_001", "canonical_smiles": "CCO"}],
                ["id", "canonical_smiles"],
            )
            pose = stage45 / "control_docking_runs" / "control_pose.pdbqt"
            pose.parent.mkdir(parents=True, exist_ok=True)
            pose.write_text("REMARK OLD POSE\nEND\n", encoding="utf-8")
            log = stage45 / "control_docking_runs" / "control.log"
            log.write_text("$ placeholder vina\n", encoding="utf-8")
            write_csv(
                stage45 / "round_1_control_docking_scores.csv",
                [
                    {
                        "panel_type": "positive_control",
                        "id": "control",
                        "smiles": "CCO",
                        "docking_score": "-7.4",
                        "pose_pass": "true",
                        "backend": "vina",
                        "receptor": str(receptor),
                        "notes": f"returncode=0;log={log};posebusters=passed",
                    }
                ],
                ["panel_type", "id", "smiles", "docking_score", "pose_pass", "backend", "receptor", "notes"],
            )
            server.write_json(
                stage45 / "round_1_control_validation.json",
                {"docking": {"status": "completed", "backend": "vina", "docking_box": {"status": "ready", "center": [0, 0, 0], "size": [20, 20, 20]}}},
            )

            payload = server.stage4_payload("project", project, 1)

            self.assertEqual(payload["docking_plan"]["status"], "not_available")
            self.assertEqual(payload["stage4_readiness_guide"]["docking"], "not_available")
            self.assertEqual(payload["docking_results"], [])
            self.assertEqual(payload["recovered_docking_evidence"]["ignored_stage45_scores"], 1)


if __name__ == "__main__":
    unittest.main()
