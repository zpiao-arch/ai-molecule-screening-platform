import argparse
import csv
import os
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


def pdb_line(record: str, serial: int, atom: str, resname: str, chain: str, resseq: int, x: float, y: float, z: float, element: str) -> str:
    return (
        f"{record:<6}{serial:5d} {atom:<4s} {resname:>3s} {chain:1s}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{20.00:6.2f}          {element:>2s}\n"
    )


class Stage4RealLibraryTests(unittest.TestCase):
    def test_stage4_real_writes_rdkit_assets_for_candidate_round(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            ai_mol_loop.ensure_project_dirs(project)
            ai_mol_loop.write_json(project / "config.json", ai_mol_loop.default_config())
            ai_mol_loop.write_csv(
                project / "candidates" / "round_1_candidates.csv",
                [
                    {"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"},
                    {"round": 1, "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "parent": "", "source": "unit"},
                    {"round": 1, "id": "bad_001", "smiles": "not_a_smiles", "parent": "", "source": "unit"},
                ],
                ["round", "id", "smiles", "parent", "source"],
            )

            ai_mol_loop.stage4_real(
                argparse.Namespace(
                    project=str(project),
                    round=1,
                    target="influenza_a_h1n1_na",
                    input_csv=None,
                    controls_csv=None,
                    top=2,
                    max_conformers=1,
                    seed=7,
                    no_sdf=False,
                )
            )

            descriptors_path = project / "stage4" / "round_1_real_descriptors.csv"
            sdf_path = project / "stage4" / "round_1_ligands.sdf"
            similarity_path = project / "stage4" / "round_1_similarity_to_controls.csv"
            assets_path = project / "stage4" / "round_1_stage4_assets.json"
            report_path = project / "reports" / "stage4_round_1_report.md"

            self.assertTrue(descriptors_path.exists())
            self.assertTrue(sdf_path.exists())
            self.assertTrue(similarity_path.exists())
            self.assertTrue(assets_path.exists())
            self.assertTrue(report_path.exists())

            descriptors = rows_from_csv(descriptors_path)
            self.assertEqual(len(descriptors), 3)
            valid_rows = [row for row in descriptors if row["valid"] == "1"]
            self.assertEqual(len(valid_rows), 2)
            self.assertTrue(all(row["descriptor_source"] == "rdkit" for row in valid_rows))
            self.assertTrue(all(row["qed"] for row in valid_rows))
            self.assertEqual([row["stage4_status"] for row in valid_rows], ["ready_for_similarity_and_sdf", "ready_for_similarity_and_sdf"])

            sdf_text = sdf_path.read_text(encoding="utf-8")
            self.assertIn("$$$$", sdf_text)
            self.assertIn("mol_001", sdf_text)

            similarity = rows_from_csv(similarity_path)
            self.assertTrue(similarity)
            self.assertTrue(any(row["control_drug"] == "oseltamivir" for row in similarity))
            self.assertTrue(all(row["fingerprint"] == "rdkit_morgan_2048_r2" for row in similarity))

            assets = json.loads(assets_path.read_text(encoding="utf-8"))
            self.assertEqual(assets["stage"], 4)
            self.assertEqual(assets["real_libraries"]["rdkit"]["status"], "available")
            self.assertIn("docking_backend", assets)
            self.assertEqual(assets["files"]["real_descriptors"], str(descriptors_path))

            report = report_path.read_text(encoding="utf-8")
            self.assertIn("Stage 4 Real Library Validation", report)
            self.assertIn("RDKit", report)

    def test_score_uses_stage4_rdkit_descriptors_when_available(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            ai_mol_loop.ensure_project_dirs(project)
            ai_mol_loop.write_json(project / "config.json", ai_mol_loop.default_config())
            ai_mol_loop.write_csv(
                project / "candidates" / "round_1_candidates.csv",
                [
                    {"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"},
                    {"round": 1, "id": "bad_001", "smiles": "not_a_smiles", "parent": "", "source": "unit"},
                ],
                ["round", "id", "smiles", "parent", "source"],
            )

            ai_mol_loop.stage4_real(
                argparse.Namespace(
                    project=str(project),
                    round=1,
                    target="influenza_a_h1n1_na",
                    input_csv=None,
                    controls_csv=None,
                    top=2,
                    max_conformers=0,
                    seed=7,
                    no_sdf=True,
                )
            )
            ai_mol_loop.score_candidates(
                argparse.Namespace(project=str(project), round=1, external_scores=None, real_descriptors=None)
            )

            scores = rows_from_csv(project / "scores" / "round_1_scores.csv")
            by_id = {row["id"]: row for row in scores}
            self.assertIn("stage4_rdkit", by_id["mol_001"]["score_source"])
            self.assertEqual(by_id["mol_001"]["validity_proxy"], "1.0")
            self.assertNotEqual(by_id["mol_001"]["qed_proxy"], "")
            self.assertIn("stage4_rdkit", by_id["bad_001"]["score_source"])
            self.assertEqual(by_id["bad_001"]["validity_proxy"], "0.0")

    def test_stage4_real_can_rescore_rank_and_feedback_in_one_command(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            ai_mol_loop.ensure_project_dirs(project)
            ai_mol_loop.write_json(project / "config.json", ai_mol_loop.default_config())
            ai_mol_loop.write_csv(
                project / "candidates" / "round_1_candidates.csv",
                [
                    {"round": 1, "id": "mol_001", "smiles": "CC(=O)Oc1ccccc1C(=O)ON1CCOCC1N", "parent": "", "source": "unit"},
                    {"round": 1, "id": "mol_002", "smiles": "CC(=O)Oc1ccccc1C(=O)ON1CCOCC1O", "parent": "", "source": "unit"},
                    {"round": 1, "id": "mol_003", "smiles": "CCS(=O)(=O)Nc1ccccc1", "parent": "", "source": "unit"},
                ],
                ["round", "id", "smiles", "parent", "source"],
            )

            ai_mol_loop.stage4_real(
                argparse.Namespace(
                    project=str(project),
                    round=1,
                    target="influenza_a_h1n1_na",
                    input_csv=None,
                    controls_csv=None,
                    top=2,
                    max_conformers=0,
                    seed=7,
                    no_sdf=True,
                    rescore=True,
                    rank_top=2,
                    feedback_top=2,
                    external_scores=None,
                )
            )

            self.assertTrue((project / "scores" / "round_1_scores.csv").exists())
            self.assertTrue((project / "ranked" / "round_1_ranked.csv").exists())
            self.assertTrue((project / "feedback" / "round_1_feedback.json").exists())
            scores = rows_from_csv(project / "scores" / "round_1_scores.csv")
            self.assertTrue(all("stage4_rdkit" in row["score_source"] for row in scores))
            feedback = json.loads((project / "feedback" / "round_1_feedback.json").read_text(encoding="utf-8"))
            self.assertGreaterEqual(feedback["selected_count"], 1)

    def test_stage4_real_writes_receptor_decoy_visual_and_docking_plan_assets(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            ai_mol_loop.ensure_project_dirs(project)
            config = ai_mol_loop.default_config()
            config["target"] = {"target_catalog_id": "influenza_a_h1n1_na", "pdb_id": "3TI6"}
            ai_mol_loop.write_json(project / "config.json", config)
            ai_mol_loop.write_csv(
                project / "candidates" / "round_1_candidates.csv",
                [
                    {"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"},
                    {"round": 1, "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "parent": "", "source": "unit"},
                ],
                ["round", "id", "smiles", "parent", "source"],
            )

            ai_mol_loop.stage4_real(
                argparse.Namespace(
                    project=str(project),
                    round=1,
                    target="influenza_a_h1n1_na",
                    input_csv=None,
                    controls_csv=None,
                    top=2,
                    max_conformers=0,
                    seed=7,
                    no_sdf=True,
                    rescore=False,
                    rank_top=0,
                    feedback_top=0,
                    external_scores=None,
                    receptor_pdb=None,
                    pdb_id="3TI6",
                    fetch_receptor=False,
                    docking_backend="auto",
                    run_docking=False,
                    decoys=4,
                    render_2d=True,
                )
            )

            receptor_path = project / "stage4" / "round_1_receptor_package.json"
            benchmark_path = project / "stage4" / "round_1_benchmark_panel.csv"
            benchmark_sdf_path = project / "stage4" / "round_1_benchmark_panel.sdf"
            decoy_path = project / "stage4" / "round_1_decoys.csv"
            docking_plan_path = project / "stage4" / "round_1_docking_plan.json"
            docking_inputs_path = project / "stage4" / "round_1_docking_inputs.csv"
            validation_path = project / "stage4" / "round_1_validation_metrics.json"
            images_dir = project / "stage4" / "round_1_2d"
            assets_path = project / "stage4" / "round_1_stage4_assets.json"

            self.assertTrue(receptor_path.exists())
            self.assertTrue(benchmark_path.exists())
            self.assertTrue(benchmark_sdf_path.exists())
            self.assertTrue(decoy_path.exists())
            self.assertTrue(docking_plan_path.exists())
            self.assertTrue(docking_inputs_path.exists())
            self.assertTrue(validation_path.exists())
            self.assertTrue(images_dir.exists())

            receptor = json.loads(receptor_path.read_text(encoding="utf-8"))
            self.assertEqual(receptor["pdb_id"], "3TI6")
            self.assertEqual(receptor["target_id"], "influenza_a_h1n1_na")
            self.assertIn("binding_site", receptor)
            self.assertIn("preparation_status", receptor)

            benchmark = rows_from_csv(benchmark_path)
            self.assertTrue(any(row["panel_type"] == "candidate" for row in benchmark))
            self.assertTrue(any(row["panel_type"] == "positive_control" for row in benchmark))
            self.assertTrue(any(row["panel_type"] == "decoy" for row in benchmark))
            benchmark_sdf = benchmark_sdf_path.read_text(encoding="utf-8")
            self.assertIn("mol_001", benchmark_sdf)
            self.assertIn("oseltamivir", benchmark_sdf)
            self.assertIn("decoy_001", benchmark_sdf)

            decoys = rows_from_csv(decoy_path)
            self.assertEqual(len(decoys), 4)
            self.assertTrue(all(row["id"].startswith("decoy_") for row in decoys))

            docking_plan = json.loads(docking_plan_path.read_text(encoding="utf-8"))
            self.assertIn(docking_plan["status"], {"not_available", "planned", "completed", "skipped"})
            self.assertIn("required_tools", docking_plan)
            self.assertIn("commands", docking_plan)

            validation = json.loads(validation_path.read_text(encoding="utf-8"))
            self.assertIn("panel_counts", validation)
            self.assertIn("control_similarity", validation)
            self.assertIn("docking_status", validation)

            pngs = list(images_dir.glob("*.png"))
            self.assertGreaterEqual(len(pngs), 1)

            assets = json.loads(assets_path.read_text(encoding="utf-8"))
            self.assertIn("receptor_package", assets["files"])
            self.assertIn("benchmark_panel", assets["files"])
            self.assertIn("benchmark_sdf", assets["files"])
            self.assertIn("docking_plan", assets["files"])
            self.assertIn("validation_metrics", assets["files"])

    def test_stage4_receptor_package_falls_back_to_project_receptor_when_fetch_fails(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            ai_mol_loop.ensure_project_dirs(project)
            receptor_dir = project / "stage4" / "receptors"
            receptor_dir.mkdir(parents=True, exist_ok=True)
            receptor = receptor_dir / "3TI6.pdb"
            receptor.write_text(
                "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 10.00           N\nEND\n",
                encoding="utf-8",
            )
            args = argparse.Namespace(receptor_pdb=None, fetch_receptor=True, pdb_id="3TI6")
            original_fetch = ai_mol_loop.fetch_pdb_file
            ai_mol_loop.fetch_pdb_file = lambda pdb_id, output: (False, "network unavailable")
            try:
                package = ai_mol_loop.write_stage4_receptor_package(
                    project / "stage4" / "round_1_receptor_package.json",
                    project,
                    1,
                    "influenza_a_h1n1_na",
                    {"target": {"target_catalog_id": "influenza_a_h1n1_na"}},
                    args,
                )
            finally:
                ai_mol_loop.fetch_pdb_file = original_fetch

            self.assertEqual(package["local_receptor_pdb"], str(receptor.resolve()))
            self.assertEqual(package["preparation_status"], "receptor_pdb_available")
            self.assertIn("local_project_receptor_found", package["fetch_status"])

    def test_stage4_real_extracts_docking_box_from_cocrystal_ligand_pdb(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            receptor = Path(raw_tmp) / "3TI6_minimal.pdb"
            receptor.write_text(
                "".join(
                    [
                        pdb_line("ATOM", 1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"),
                        pdb_line("HETATM", 2, "C1", "G39", "A", 501, 10.0, 20.0, 30.0, "C"),
                        pdb_line("HETATM", 3, "C2", "G39", "A", 501, 12.0, 20.0, 30.0, "C"),
                        pdb_line("HETATM", 4, "O1", "G39", "A", 501, 10.0, 24.0, 30.0, "O"),
                        pdb_line("HETATM", 5, "O", "HOH", "A", 601, 99.0, 99.0, 99.0, "O"),
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            ai_mol_loop.ensure_project_dirs(project)
            config = ai_mol_loop.default_config()
            config["target"] = {"target_catalog_id": "influenza_a_h1n1_na", "pdb_id": "3TI6"}
            ai_mol_loop.write_json(project / "config.json", config)
            ai_mol_loop.write_csv(
                project / "candidates" / "round_1_candidates.csv",
                [{"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"}],
                ["round", "id", "smiles", "parent", "source"],
            )

            ai_mol_loop.stage4_real(
                argparse.Namespace(
                    project=str(project),
                    round=1,
                    target="influenza_a_h1n1_na",
                    input_csv=None,
                    controls_csv=None,
                    top=1,
                    max_conformers=0,
                    seed=7,
                    no_sdf=True,
                    rescore=False,
                    rank_top=0,
                    feedback_top=0,
                    external_scores=None,
                    receptor_pdb=str(receptor),
                    pdb_id="3TI6",
                    fetch_receptor=False,
                    docking_backend="vina",
                    run_docking=False,
                    docking_timeout=5,
                    decoys=0,
                    render_2d=False,
                )
            )

            receptor_package = json.loads((project / "stage4" / "round_1_receptor_package.json").read_text(encoding="utf-8"))
            site = receptor_package["binding_site"]
            self.assertEqual(site["source"], "co_crystal_ligand")
            self.assertEqual(site["reference_ligand"], "oseltamivir")
            self.assertEqual(site["detected_ligand"]["resname"], "G39")
            self.assertEqual(site["detected_ligand"]["chain"], "A")
            self.assertEqual(site["detected_ligand"]["resseq"], 501)
            self.assertEqual(site["center"], [10.667, 21.333, 30.0])
            self.assertEqual(site["size"], [20.0, 20.0, 20.0])
            self.assertEqual(site["status"], "extracted_from_receptor")

            docking_plan = json.loads((project / "stage4" / "round_1_docking_plan.json").read_text(encoding="utf-8"))
            self.assertEqual(docking_plan["docking_box"]["status"], "ready")
            self.assertEqual(docking_plan["docking_box"]["center"], [10.667, 21.333, 30.0])
            self.assertEqual(docking_plan["docking_box"]["size"], [20.0, 20.0, 20.0])

    def test_stage4_real_runs_vina_when_backend_and_prepared_inputs_exist(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            project = root / "project"
            fake_bin = root / "bin"
            receptor = root / "receptor.pdb"
            receptor_pdbqt = root / "receptor.pdbqt"
            fake_bin.mkdir(parents=True)
            receptor.write_text("HEADER TEST RECEPTOR\nEND\n", encoding="utf-8")
            receptor_pdbqt.write_text("REMARK TEST RECEPTOR PDBQT\nEND\n", encoding="utf-8")
            vina = fake_bin / "vina"
            vina.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '--log' ]; then echo 'unexpected --log' >&2; exit 64; fi\n"
                "  if [ \"$1\" = '--out' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK TEST POSE\\nEND\\n' > \"$out\"\n"
                "echo '-----+------------+----------+----------'\n"
                "echo '   1        -7.4      0.000      0.000'\n"
                "exit 0\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            ligand_prep = fake_bin / "mk_prepare_ligand.py"
            ligand_prep.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-o' ] || [ \"$1\" = '--out' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK TEST LIGAND PDBQT\\nEND\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            ligand_prep.chmod(0o755)
            obabel = fake_bin / "obabel"
            obabel.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'POSE SDF\\n$$$$\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            obabel.chmod(0o755)
            bust = fake_bin / "bust"
            bust.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '--output' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'mol_pred,all_atoms_connected\\npose,false\\n' > \"$out\"\n"
                "exit 1\n",
                encoding="utf-8",
            )
            bust.chmod(0o755)
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = f"{fake_bin}{os.pathsep}{old_path}"
            try:
                ai_mol_loop.ensure_project_dirs(project)
                config = ai_mol_loop.default_config()
                config["target"] = {
                    "target_catalog_id": "influenza_a_h1n1_na",
                    "pocket": {"center": [0.0, 0.0, 0.0], "size": [20.0, 20.0, 20.0]},
                }
                ai_mol_loop.write_json(project / "config.json", config)
                ai_mol_loop.write_csv(
                    project / "candidates" / "round_1_candidates.csv",
                    [{"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"}],
                    ["round", "id", "smiles", "parent", "source"],
                )

                ai_mol_loop.stage4_real(
                    argparse.Namespace(
                        project=str(project),
                        round=1,
                        target="influenza_a_h1n1_na",
                        input_csv=None,
                        controls_csv=None,
                        top=1,
                        max_conformers=0,
                        seed=7,
                        no_sdf=False,
                        rescore=False,
                        rank_top=0,
                        feedback_top=0,
                        external_scores=None,
                        receptor_pdb=str(receptor),
                        pdb_id="",
                        fetch_receptor=False,
                        docking_backend="vina",
                        run_docking=True,
                        docking_timeout=5,
                        decoys=0,
                        render_2d=False,
                    )
                )
            finally:
                os.environ["PATH"] = old_path

            plan = json.loads((project / "stage4" / "round_1_docking_plan.json").read_text(encoding="utf-8"))
            self.assertEqual(plan["status"], "completed")
            self.assertEqual(plan["selected_backend"], "vina")
            scores = rows_from_csv(project / "stage4" / "round_1_docking_scores_template.csv")
            self.assertTrue(scores)
            self.assertEqual(scores[0]["docking_score"], "-7.4")
            self.assertEqual(scores[0]["pose_pass"], "")
            self.assertIn("posebusters=", scores[0]["notes"])
            self.assertIn("log=", scores[0]["notes"])
            log_path = project / "stage4" / "docking_runs" / "mol_001.log"
            self.assertTrue(log_path.exists())
            self.assertIn("-7.4", log_path.read_text(encoding="utf-8"))

    def test_stage4_real_rescore_uses_generated_docking_scores_by_default(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            project = root / "project"
            fake_bin = root / "bin"
            receptor = root / "receptor.pdb"
            receptor_pdbqt = root / "receptor.pdbqt"
            fake_bin.mkdir(parents=True)
            receptor.write_text("HEADER TEST RECEPTOR\nEND\n", encoding="utf-8")
            receptor_pdbqt.write_text("REMARK TEST RECEPTOR PDBQT\nEND\n", encoding="utf-8")
            vina = fake_bin / "vina"
            vina.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '--out' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK TEST POSE\\nEND\\n' > \"$out\"\n"
                "echo '-----+------------+----------+----------'\n"
                "echo '   1        -8.8      0.000      0.000'\n"
                "exit 0\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            ligand_prep = fake_bin / "mk_prepare_ligand.py"
            ligand_prep.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-o' ] || [ \"$1\" = '--out' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK TEST LIGAND PDBQT\\nEND\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            ligand_prep.chmod(0o755)
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = f"{fake_bin}{os.pathsep}{old_path}"
            try:
                ai_mol_loop.ensure_project_dirs(project)
                config = ai_mol_loop.default_config()
                config["target"] = {
                    "target_catalog_id": "influenza_a_h1n1_na",
                    "pocket": {"center": [0.0, 0.0, 0.0], "size": [20.0, 20.0, 20.0]},
                }
                ai_mol_loop.write_json(project / "config.json", config)
                ai_mol_loop.write_csv(
                    project / "candidates" / "round_1_candidates.csv",
                    [{"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"}],
                    ["round", "id", "smiles", "parent", "source"],
                )

                ai_mol_loop.stage4_real(
                    argparse.Namespace(
                        project=str(project),
                        round=1,
                        target="influenza_a_h1n1_na",
                        input_csv=None,
                        controls_csv=None,
                        top=1,
                        max_conformers=0,
                        seed=7,
                        no_sdf=False,
                        rescore=True,
                        rank_top=1,
                        feedback_top=1,
                        external_scores=None,
                        receptor_pdb=str(receptor),
                        pdb_id="",
                        fetch_receptor=False,
                        docking_backend="vina",
                        run_docking=True,
                        docking_timeout=5,
                        decoys=0,
                        render_2d=False,
                    )
                )
            finally:
                os.environ["PATH"] = old_path

            score_rows = rows_from_csv(project / "scores" / "round_1_scores.csv")
            self.assertTrue(score_rows)
            self.assertIn("external_docking", score_rows[0]["score_source"])
            self.assertEqual(score_rows[0]["raw_docking_kcal_mol"], "-8.8")
            self.assertTrue((project / "ranked" / "round_1_ranked.csv").exists())
            self.assertTrue((project / "feedback" / "round_1_feedback.json").exists())

    def test_stage4_real_runs_posebusters_after_vina_pose_output(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            project = root / "project"
            fake_bin = root / "bin"
            receptor = root / "receptor.pdb"
            receptor_pdbqt = root / "receptor.pdbqt"
            fake_bin.mkdir(parents=True)
            receptor.write_text("HEADER TEST RECEPTOR\nEND\n", encoding="utf-8")
            receptor_pdbqt.write_text("REMARK TEST RECEPTOR PDBQT\nEND\n", encoding="utf-8")
            vina = fake_bin / "vina"
            vina.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '--out' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK TEST POSE PDBQT\\nEND\\n' > \"$out\"\n"
                "echo '-----+------------+----------+----------'\n"
                "echo '   1        -7.7      0.000      0.000'\n"
                "exit 0\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            ligand_prep = fake_bin / "mk_prepare_ligand.py"
            ligand_prep.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-o' ] || [ \"$1\" = '--out' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK TEST LIGAND PDBQT\\nEND\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            ligand_prep.chmod(0o755)
            obabel = fake_bin / "obabel"
            obabel.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'pose sdf\\n$$$$\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            obabel.chmod(0o755)
            bust = fake_bin / "bust"
            bust.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '--output' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'file,all_checks\\npose.sdf,True\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            bust.chmod(0o755)
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = f"{fake_bin}{os.pathsep}{old_path}"
            try:
                ai_mol_loop.ensure_project_dirs(project)
                config = ai_mol_loop.default_config()
                config["target"] = {
                    "target_catalog_id": "influenza_a_h1n1_na",
                    "pocket": {"center": [0.0, 0.0, 0.0], "size": [20.0, 20.0, 20.0]},
                }
                ai_mol_loop.write_json(project / "config.json", config)
                ai_mol_loop.write_csv(
                    project / "candidates" / "round_1_candidates.csv",
                    [{"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"}],
                    ["round", "id", "smiles", "parent", "source"],
                )

                ai_mol_loop.stage4_real(
                    argparse.Namespace(
                        project=str(project),
                        round=1,
                        target="influenza_a_h1n1_na",
                        input_csv=None,
                        controls_csv=None,
                        top=1,
                        max_conformers=0,
                        seed=7,
                        no_sdf=False,
                        rescore=False,
                        rank_top=0,
                        feedback_top=0,
                        external_scores=None,
                        receptor_pdb=str(receptor),
                        pdb_id="",
                        fetch_receptor=False,
                        docking_backend="vina",
                        run_docking=True,
                        docking_timeout=5,
                        decoys=0,
                        render_2d=False,
                    )
                )
            finally:
                os.environ["PATH"] = old_path

            scores = rows_from_csv(project / "stage4" / "round_1_docking_scores_template.csv")
            self.assertEqual(scores[0]["docking_score"], "-7.7")
            self.assertEqual(scores[0]["pose_pass"], "true")
            self.assertIn("posebusters=passed", scores[0]["notes"])
            self.assertTrue((project / "stage4" / "docking_runs" / "mol_001_posebusters.csv").exists())

    def test_stage4_rejects_vina_score_when_pose_output_is_missing(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            fake_bin = root / "bin"
            fake_bin.mkdir(parents=True)
            vina = fake_bin / "vina"
            vina.write_text(
                "#!/bin/sh\n"
                "echo '-----+------------+----------+----------'\n"
                "echo '   1        -7.4      0.000      0.000'\n"
                "exit 0\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            plan = {
                "status": "planned",
                "selected_backend": "vina",
                "tool_status": {"vina": {"status": "found", "path": str(vina)}},
                "docking_box": {"status": "ready", "center": [0.0, 0.0, 0.0], "size": [20.0, 20.0, 20.0]},
            }
            receptor = root / "receptor.pdbqt"
            ligand = root / "ligand.pdbqt"
            receptor.write_text("REMARK RECEPTOR\nEND\n", encoding="utf-8")
            ligand.write_text("REMARK LIGAND\nEND\n", encoding="utf-8")
            output_csv = root / "scores.csv"

            plan, results = ai_mol_loop.run_stage4_docking_if_possible(
                plan,
                [
                    {
                        "panel_type": "candidate",
                        "id": "mol_001",
                        "canonical_smiles": "CCO",
                        "ligand_pdbqt": str(ligand),
                        "receptor_pdbqt": str(receptor),
                    }
                ],
                [{"id": "mol_001", "canonical_smiles": "CCO"}],
                output_csv,
                timeout=5,
            )

            self.assertEqual(plan["status"], "attempted_no_scores")
            self.assertEqual(results, [])
            self.assertFalse(output_csv.exists())
            self.assertIn("missing_pose_output", plan["message"])

    def test_stage4_capabilities_marks_placeholder_vina_invalid(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            fake_bin = root / "bin"
            fake_bin.mkdir(parents=True)
            vina = fake_bin / "vina"
            vina.write_text(
                "#!/bin/sh\n"
                "echo '-----+------------+----------+----------'\n"
                "echo '   1        -7.4      0.000      0.000'\n"
                "exit 0\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = f"{fake_bin}{os.pathsep}{old_path}"
            try:
                capabilities = ai_mol_loop.stage4_capabilities()
            finally:
                os.environ["PATH"] = old_path

            vina_status = capabilities["executables"]["vina"]
            self.assertEqual(vina_status["status"], "invalid")
            self.assertEqual(vina_status["validation"], "placeholder_no_pose_output")
            self.assertNotIn("vina", capabilities["docking_backend"]["available_backends"])

    def test_stage4_ligand_preparation_falls_back_to_openbabel_when_meeko_fails(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            fake_bin = root / "bin"
            output_dir = root / "prepared"
            fake_bin.mkdir(parents=True)
            meeko = fake_bin / "mk_prepare_ligand.py"
            meeko.write_text("#!/bin/sh\necho meeko failed >&2\nexit 2\n", encoding="utf-8")
            meeko.chmod(0o755)
            obabel = fake_bin / "obabel"
            obabel.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK OPENBABEL PDBQT\\nEND\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            obabel.chmod(0o755)
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = f"{fake_bin}{os.pathsep}{old_path}"
            try:
                libs = ai_mol_loop.require_rdkit()
                capabilities = ai_mol_loop.stage4_capabilities()
                records, attempts = ai_mol_loop.stage4_prepare_ligand_records(
                    [
                        {
                            "panel_type": "candidate",
                            "id": "mol_001",
                            "canonical_smiles": "CCOC(=O)c1ccccc1",
                        }
                    ],
                    libs,
                    capabilities,
                    output_dir,
                    seed=7,
                    timeout=5,
                    enabled=True,
                )
            finally:
                os.environ["PATH"] = old_path

            self.assertEqual(records["mol_001"]["status"], "ready_for_docking")
            self.assertEqual(records["mol_001"]["tool"], "openbabel")
            self.assertTrue(Path(records["mol_001"]["ligand_pdbqt"]).exists())
            self.assertTrue(any(item["status"] == "ligand_prepare_failed" and item["tool"] == "meeko" for item in attempts))
            self.assertTrue(any(item["status"] == "ready_for_docking" and item["tool"] == "openbabel" for item in attempts))

    def test_stage4_receptor_preparation_falls_back_to_clean_protein_pdb(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            fake_bin = root / "bin"
            receptor = root / "3TI6_raw.pdb"
            fake_bin.mkdir(parents=True)
            receptor.write_text(
                "".join(
                    [
                        pdb_line("ATOM", 1, "N", "SER", "A", 82, -47.086, 46.266, 27.606, "N"),
                        pdb_line("ATOM", 2, "CA", "SER", "A", 82, -45.951, 45.359, 27.470, "C"),
                        pdb_line("HETATM", 3, "C1", "G39", "A", 801, -28.899, 12.358, 24.407, "C"),
                        pdb_line("HETATM", 4, "O", "HOH", "A", 901, 99.0, 99.0, 99.0, "O"),
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            receptor_tool = fake_bin / "mk_prepare_receptor.py"
            receptor_tool.write_text(
                "#!/bin/sh\n"
                "input=''\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '--read_pdb' ]; then shift; input=\"$1\"; fi\n"
                "  if [ \"$1\" = '-p' ] || [ \"$1\" = '--write_pdbqt' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "if grep -q '^HETATM' \"$input\"; then echo 'unknown hetero residue' >&2; exit 2; fi\n"
                "printf 'REMARK CLEAN RECEPTOR PDBQT\\nEND\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            receptor_tool.chmod(0o755)
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = f"{fake_bin}{os.pathsep}{old_path}"
            try:
                capabilities = ai_mol_loop.stage4_capabilities()
                receptor_package = {"local_receptor_pdb": str(receptor), "local_receptor_pdbqt": "", "preparation_status": "receptor_pdb_available"}
                result = ai_mol_loop.stage4_prepare_receptor_pdbqt(
                    receptor_package,
                    capabilities,
                    {"status": "ready", "center": [1.0, 2.0, 3.0], "size": [20.0, 20.0, 20.0]},
                    timeout=5,
                )
            finally:
                os.environ["PATH"] = old_path

            self.assertEqual(result["status"], "ready")
            self.assertEqual(result["fallback"], "protein_only_clean_pdb")
            self.assertTrue(Path(result["receptor_pdbqt"]).exists())
            cleaned = Path(result["cleaned_receptor_pdb"])
            self.assertTrue(cleaned.exists())
            cleaned_text = cleaned.read_text(encoding="utf-8")
            self.assertIn("ATOM", cleaned_text)
            self.assertNotIn("HETATM", cleaned_text)
            self.assertEqual(receptor_package["preparation_status"], "receptor_pdbqt_available")

    def test_stage4_receptor_preparation_falls_back_to_openbabel_when_meeko_receptor_fails(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            fake_bin = root / "bin"
            receptor = root / "3TI6_raw.pdb"
            fake_bin.mkdir(parents=True)
            receptor.write_text(
                "".join(
                    [
                        pdb_line("ATOM", 1, "N", "SER", "A", 82, -47.086, 46.266, 27.606, "N"),
                        pdb_line("ATOM", 2, "CA", "SER", "A", 82, -45.951, 45.359, 27.470, "C"),
                        pdb_line("HETATM", 3, "C1", "G39", "A", 801, -28.899, 12.358, 24.407, "C"),
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            receptor_tool = fake_bin / "mk_prepare_receptor.py"
            receptor_tool.write_text("#!/bin/sh\necho meeko receptor failed >&2\nexit 2\n", encoding="utf-8")
            receptor_tool.chmod(0o755)
            obabel = fake_bin / "obabel"
            obabel.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK OPENBABEL RECEPTOR PDBQT\\nEND\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            obabel.chmod(0o755)
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = f"{fake_bin}{os.pathsep}{old_path}"
            try:
                capabilities = ai_mol_loop.stage4_capabilities()
                receptor_package = {"local_receptor_pdb": str(receptor), "local_receptor_pdbqt": "", "preparation_status": "receptor_pdb_available"}
                result = ai_mol_loop.stage4_prepare_receptor_pdbqt(
                    receptor_package,
                    capabilities,
                    {"status": "ready", "center": [1.0, 2.0, 3.0], "size": [20.0, 20.0, 20.0]},
                    timeout=5,
                )
            finally:
                os.environ["PATH"] = old_path

            self.assertEqual(result["status"], "ready")
            self.assertEqual(result["tool"], "openbabel")
            self.assertEqual(result["fallback"], "openbabel_protein_only_pdbqt")
            self.assertTrue(Path(result["receptor_pdbqt"]).exists())
            self.assertTrue(Path(result["cleaned_receptor_pdb"]).exists())
            self.assertEqual(receptor_package["local_receptor_pdbqt"], result["receptor_pdbqt"])

    def test_stage45_validate_controls_builds_candidate_control_decoy_calibration_assets(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            root = Path(raw_tmp)
            project = root / "project"
            fake_bin = root / "bin"
            receptor = root / "3TI6_minimal.pdb"
            fake_bin.mkdir(parents=True)
            receptor.write_text(
                "".join(
                    [
                        pdb_line("ATOM", 1, "N", "SER", "A", 82, -47.086, 46.266, 27.606, "N"),
                        pdb_line("ATOM", 2, "CA", "SER", "A", 82, -45.951, 45.359, 27.470, "C"),
                        pdb_line("HETATM", 3, "C1", "G39", "A", 801, -28.899, 12.358, 24.407, "C"),
                        pdb_line("HETATM", 4, "C2", "G39", "A", 801, -27.899, 12.358, 24.407, "C"),
                        pdb_line("HETATM", 5, "O1", "G39", "A", 801, -28.899, 13.358, 24.407, "O"),
                        pdb_line("HETATM", 6, "N1", "G39", "A", 801, -28.899, 12.358, 25.407, "N"),
                        "END\n",
                    ]
                ),
                encoding="utf-8",
            )
            (root / "3TI6_minimal.pdbqt").write_text("REMARK RECEPTOR PDBQT\nEND\n", encoding="utf-8")
            vina = fake_bin / "vina"
            vina.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "lig=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '--out' ]; then shift; out=\"$1\"; fi\n"
                "  if [ \"$1\" = '--ligand' ]; then shift; lig=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK TEST POSE\\nEND\\n' > \"$out\"\n"
                "case \"$lig\" in\n"
                "  *oseltamivir*) score='-8.4' ;;\n"
                "  *decoy*) score='-4.2' ;;\n"
                "  *) score='-6.7' ;;\n"
                "esac\n"
                "echo '-----+------------+----------+----------'\n"
                "echo \"   1        $score      0.000      0.000\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            vina.chmod(0o755)
            ligand_prep = fake_bin / "mk_prepare_ligand.py"
            ligand_prep.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-o' ] || [ \"$1\" = '--out' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'REMARK TEST LIGAND PDBQT\\nEND\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            ligand_prep.chmod(0o755)
            obabel = fake_bin / "obabel"
            obabel.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '-O' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'pose sdf\\n$$$$\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            obabel.chmod(0o755)
            bust = fake_bin / "bust"
            bust.write_text(
                "#!/bin/sh\n"
                "out=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  if [ \"$1\" = '--output' ]; then shift; out=\"$1\"; fi\n"
                "  shift\n"
                "done\n"
                "printf 'file,all_checks\\npose.sdf,True\\n' > \"$out\"\n"
                "exit 0\n",
                encoding="utf-8",
            )
            bust.chmod(0o755)
            old_path = os.environ.get("PATH", "")
            os.environ["PATH"] = f"{fake_bin}{os.pathsep}{old_path}"
            try:
                ai_mol_loop.ensure_project_dirs(project)
                config = ai_mol_loop.default_config()
                config["target"] = {"target_catalog_id": "influenza_a_h1n1_na"}
                ai_mol_loop.write_json(project / "config.json", config)
                ai_mol_loop.write_csv(
                    project / "candidates" / "round_1_candidates.csv",
                    [
                        {"round": 1, "id": "mol_001", "smiles": "CCOC(=O)c1ccccc1", "parent": "", "source": "unit"},
                        {"round": 1, "id": "mol_002", "smiles": "CCS(=O)(=O)Nc1ccccc1", "parent": "", "source": "unit"},
                    ],
                    ["round", "id", "smiles", "parent", "source"],
                )
                ai_mol_loop.stage4_real(
                    argparse.Namespace(
                        project=str(project),
                        round=1,
                        target="influenza_a_h1n1_na",
                        input_csv=None,
                        controls_csv=None,
                        top=2,
                        max_conformers=0,
                        seed=7,
                        no_sdf=False,
                        rescore=True,
                        rank_top=2,
                        feedback_top=2,
                        external_scores=None,
                        receptor_pdb=str(receptor),
                        pdb_id="3TI6",
                        fetch_receptor=False,
                        docking_backend="vina",
                        run_docking=True,
                        docking_timeout=5,
                        decoys=2,
                        render_2d=False,
                    )
                )

                ai_mol_loop.stage45_validate_controls(
                    argparse.Namespace(
                        project=str(project),
                        round=1,
                        target="influenza_a_h1n1_na",
                        top_candidates=2,
                        decoys=2,
                        controls_csv=None,
                        docking_backend="vina",
                        docking_timeout=5,
                        seed=11,
                        no_docking=False,
                    )
                )
            finally:
                os.environ["PATH"] = old_path

            scores_path = project / "stage4_5" / "round_1_control_docking_scores.csv"
            validation_path = project / "stage4_5" / "round_1_control_validation.json"
            report_path = project / "reports" / "stage4_5_round_1_control_validation.md"
            self.assertTrue(scores_path.exists())
            self.assertTrue(validation_path.exists())
            self.assertTrue(report_path.exists())

            rows = rows_from_csv(scores_path)
            panel_types = {row["panel_type"] for row in rows}
            self.assertIn("candidate", panel_types)
            self.assertIn("positive_control", panel_types)
            self.assertIn("decoy", panel_types)
            self.assertTrue(any(row["id"] == "oseltamivir" and row["docking_score"] == "-8.4" for row in rows))
            self.assertTrue(any(row["id"] == "mol_001" and row["docking_score"] == "-6.7" for row in rows))
            self.assertTrue(all("relative_to_best_control" in row for row in rows))

            validation = json.loads(validation_path.read_text(encoding="utf-8"))
            self.assertEqual(validation["stage"], 4.5)
            self.assertEqual(validation["target_id"], "influenza_a_h1n1_na")
            self.assertEqual(validation["counts"]["candidate"], 2)
            self.assertGreaterEqual(validation["counts"]["known_control"], 4)
            self.assertEqual(validation["counts"]["decoy"], 2)
            self.assertEqual(validation["docking"]["status"], "completed")
            self.assertEqual(validation["best_known_control"]["id"], "oseltamivir")
            self.assertEqual(validation["best_candidate"]["id"], "mol_001")
            self.assertIn("candidate_vs_controls", validation)
            self.assertIn("redocking", validation)
            self.assertIn("boundary", validation)
            self.assertTrue((project / "stage4_5" / "reference_ligand_G39_A_801.pdb").exists())

            report = report_path.read_text(encoding="utf-8")
            self.assertIn("Stage 4.5 Control Calibration", report)
            self.assertIn("known controls calibrate the docking workflow", report)
            self.assertIn("does not prove antiviral efficacy", report)

    def test_stage46_retrospective_benchmark_scores_controls_against_decoys(self):
        with tempfile.TemporaryDirectory(dir=Path.cwd()) as raw_tmp:
            project = Path(raw_tmp) / "project"
            ai_mol_loop.ensure_project_dirs(project)
            stage45 = project / "stage4_5"
            stage45.mkdir(parents=True, exist_ok=True)
            ai_mol_loop.write_csv(
                stage45 / "round_1_control_docking_scores.csv",
                [
                    {
                        "panel_type": "positive_control",
                        "id": "oseltamivir",
                        "smiles": "CCOC(=O)C1=C(C)OC(CC)C(C1NC(C)=O)N",
                        "docking_score": "-8.4",
                        "pose_pass": "true",
                        "backend": "vina",
                        "receptor": "3TI6",
                        "relative_to_best_control": "0.0",
                        "score_band": "strong",
                        "notes": "known active",
                    },
                    {
                        "panel_type": "positive_control",
                        "id": "zanamivir",
                        "smiles": "CC(=O)NC1C(N)C(O)OC(C(O)C(O)CO)C1C(=O)O",
                        "docking_score": "-7.9",
                        "pose_pass": "true",
                        "backend": "vina",
                        "receptor": "3TI6",
                        "relative_to_best_control": "0.5",
                        "score_band": "strong",
                        "notes": "known active",
                    },
                    {
                        "panel_type": "reference_control",
                        "id": "laninamivir",
                        "smiles": "CC(=O)NC1C(N)C(OC)OC(C(O)C(O)CO)C1C(=O)O",
                        "docking_score": "-7.6",
                        "pose_pass": "true",
                        "backend": "vina",
                        "receptor": "3TI6",
                        "relative_to_best_control": "0.8",
                        "score_band": "strong",
                        "notes": "known active",
                    },
                    {
                        "panel_type": "candidate",
                        "id": "mol_001",
                        "smiles": "CCOC(=O)c1ccccc1",
                        "docking_score": "-6.7",
                        "pose_pass": "true",
                        "backend": "vina",
                        "receptor": "3TI6",
                        "relative_to_best_control": "1.7",
                        "score_band": "moderate",
                        "notes": "candidate",
                    },
                    {
                        "panel_type": "decoy",
                        "id": "decoy_001",
                        "smiles": "CCO",
                        "docking_score": "-4.0",
                        "pose_pass": "true",
                        "backend": "vina",
                        "receptor": "3TI6",
                        "relative_to_best_control": "4.4",
                        "score_band": "weak",
                        "notes": "decoy",
                    },
                    {
                        "panel_type": "decoy",
                        "id": "decoy_002",
                        "smiles": "c1ccccc1",
                        "docking_score": "-5.0",
                        "pose_pass": "true",
                        "backend": "vina",
                        "receptor": "3TI6",
                        "relative_to_best_control": "3.4",
                        "score_band": "weak",
                        "notes": "decoy",
                    },
                    {
                        "panel_type": "decoy",
                        "id": "decoy_003",
                        "smiles": "CCN(CC)CC",
                        "docking_score": "-3.5",
                        "pose_pass": "true",
                        "backend": "vina",
                        "receptor": "3TI6",
                        "relative_to_best_control": "4.9",
                        "score_band": "weak",
                        "notes": "decoy",
                    },
                ],
                ai_mol_loop.STAGE45_DOCKING_SCORE_FIELDS,
            )

            ai_mol_loop.stage46_retrospective_benchmark(
                argparse.Namespace(project=str(project), round=1, positive_types="positive_control,reference_control,control", negative_types="decoy", top_k="1,3,5")
            )

            benchmark_path = project / "stage4_6" / "round_1_retrospective_benchmark.json"
            ranking_path = project / "stage4_6" / "round_1_retrospective_ranking.csv"
            report_path = project / "reports" / "stage4_6_round_1_retrospective_benchmark.md"
            self.assertTrue(benchmark_path.exists())
            self.assertTrue(ranking_path.exists())
            self.assertTrue(report_path.exists())

            benchmark = json.loads(benchmark_path.read_text(encoding="utf-8"))
            self.assertEqual(benchmark["stage"], 4.6)
            self.assertEqual(benchmark["counts"]["positives"], 3)
            self.assertEqual(benchmark["counts"]["negatives"], 3)
            self.assertEqual(benchmark["metrics"]["roc_auc"], 1.0)
            self.assertEqual(benchmark["metrics"]["top_k"]["top_1"]["control_hit_count"], 1)
            self.assertEqual(benchmark["metrics"]["top_k"]["top_3"]["control_hit_count"], 3)
            self.assertEqual(benchmark["metrics"]["best_candidate"]["id"], "mol_001")
            self.assertEqual(benchmark["metrics"]["best_candidate"]["rank"], 4)

            ranking = rows_from_csv(ranking_path)
            self.assertEqual(ranking[0]["id"], "oseltamivir")
            self.assertEqual(ranking[0]["label"], "positive")
            self.assertEqual(ranking[-1]["label"], "negative")

            report = report_path.read_text(encoding="utf-8")
            self.assertIn("Stage 4.6 Retrospective Benchmark", report)
            self.assertIn("ROC-AUC", report)
            self.assertIn("not wet-lab validation", report)


if __name__ == "__main__":
    unittest.main()
