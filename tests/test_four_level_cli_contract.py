import importlib.util
import importlib
import builtins
import inspect
import json
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
PACKAGE_DIR = ROOT / "scientific_validation" / "four_level_cli_1kx10k"
SCORING_DIR = ROOT / "scoring"


def _load_runtime_doctor():
    module_path = PACKAGE_DIR / "runtime_doctor.py"
    assert module_path.is_file(), f"runtime doctor missing: {module_path}"
    spec = importlib.util.spec_from_file_location("four_level_runtime_doctor", module_path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_asset_integrity():
    module_path = SCORING_DIR / "asset_integrity.py"
    assert module_path.is_file(), f"asset integrity module missing: {module_path}"
    spec = importlib.util.spec_from_file_location("four_level_asset_integrity", module_path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _load_module_from_file(name: str, module_path: Path):
    spec = importlib.util.spec_from_file_location(name, module_path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_asset_manifest(root: Path, relative: str, digest: str) -> Path:
    manifest = root / "ASSET_MANIFEST.json"
    manifest.write_text(
        json.dumps({"schema": "test", "files": {relative: digest}}),
        encoding="utf-8",
    )
    return manifest


def test_asset_verification_accepts_manifested_file(tmp_path, monkeypatch):
    asset_integrity = _load_asset_integrity()
    model = tmp_path / "scoring" / "models" / "bindingdb_l2" / "model.joblib"
    model.parent.mkdir(parents=True)
    model.write_bytes(b"trusted model")
    digest = asset_integrity.sha256_file(model)
    manifest = _write_asset_manifest(
        tmp_path,
        "scoring/models/bindingdb_l2/model.joblib",
        digest,
    )
    monkeypatch.setenv("FOUR_LEVEL_ASSET_MANIFEST", str(manifest))

    result = asset_integrity.verify_asset(model)

    assert result["verified"] is True
    assert result["sha256"] == digest
    assert result["relative_path"] == "scoring/models/bindingdb_l2/model.joblib"


def test_asset_verification_rejects_mismatch(tmp_path, monkeypatch):
    asset_integrity = _load_asset_integrity()
    model = tmp_path / "scoring" / "models" / "bindingdb_l2" / "model.joblib"
    model.parent.mkdir(parents=True)
    model.write_bytes(b"tampered model")
    manifest = _write_asset_manifest(
        tmp_path,
        "scoring/models/bindingdb_l2/model.joblib",
        "0" * 64,
    )
    monkeypatch.setenv("FOUR_LEVEL_ASSET_MANIFEST", str(manifest))

    with pytest.raises(asset_integrity.AssetIntegrityError, match="sha256 mismatch"):
        asset_integrity.verify_asset(model)


def test_asset_verification_rejects_unmanifested_file(tmp_path, monkeypatch):
    asset_integrity = _load_asset_integrity()
    model = tmp_path / "custom.joblib"
    model.write_bytes(b"custom model")
    manifest = _write_asset_manifest(tmp_path, "different.joblib", "0" * 64)
    monkeypatch.setenv("FOUR_LEVEL_ASSET_MANIFEST", str(manifest))

    with pytest.raises(asset_integrity.AssetIntegrityError, match="not listed"):
        asset_integrity.verify_asset(model)


def test_asset_verification_requires_explicit_unsafe_override(tmp_path, monkeypatch):
    asset_integrity = _load_asset_integrity()
    model = tmp_path / "custom.joblib"
    model.write_bytes(b"custom model")
    monkeypatch.delenv("FOUR_LEVEL_ASSET_MANIFEST", raising=False)
    monkeypatch.delenv("FOUR_LEVEL_ALLOW_UNVERIFIED_ASSETS", raising=False)

    with pytest.raises(asset_integrity.AssetIntegrityError, match="manifest"):
        asset_integrity.verify_asset(model)

    monkeypatch.setenv("FOUR_LEVEL_ALLOW_UNVERIFIED_ASSETS", "1")
    result = asset_integrity.verify_asset(model)
    assert result["verified"] is False
    assert result["unsafe_override"] is True


def test_bindingdb_load_checks_asset_before_joblib_deserialization(tmp_path, monkeypatch):
    sys.path.insert(0, str(SCORING_DIR))
    import l2_bindingdb

    model = tmp_path / "model.joblib"
    model.write_bytes(b"not a joblib bundle")
    asset_integrity = _load_asset_integrity()
    monkeypatch.setattr(
        l2_bindingdb,
        "verify_asset",
        lambda path: (_ for _ in ()).throw(asset_integrity.AssetIntegrityError("blocked before load")),
    )
    scorer = l2_bindingdb.Layer2BindingDB(model_path=str(model))

    with pytest.raises(asset_integrity.AssetIntegrityError, match="blocked before load"):
        scorer._ensure_model()


def test_runtime_doctor_honors_individual_binary_paths(tmp_path, monkeypatch):
    runtime_doctor = _load_runtime_doctor()
    scoring_dir = tmp_path / "repo" / "scoring"
    scoring_dir.mkdir(parents=True)
    binaries = {}
    for name in ("smina", "obabel"):
        binary = tmp_path / name
        binary.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        binary.chmod(0o755)
        binaries[name] = binary
        monkeypatch.setenv(f"{name.upper()}_BIN", str(binary))
    monkeypatch.delenv("DOCK_BIN_DIR", raising=False)

    result = runtime_doctor.inspect_runtime(scoring_dir=scoring_dir)

    assert result["binaries"]["smina"]["path"] == str(binaries["smina"].resolve())
    assert result["binaries"]["obabel"]["path"] == str(binaries["obabel"].resolve())


@pytest.mark.integration_assets
@pytest.mark.integration_docking
def test_runtime_doctor_requires_all_four_backends():
    runtime_doctor = _load_runtime_doctor()
    result = runtime_doctor.inspect_runtime(scoring_dir=SCORING_DIR)

    assert result["python"][:2] == [3, 11]
    assert result["models"]["l2"]["exists"] is True
    assert result["models"]["l2_compat_sklearn_1_7_2"]["exists"] is True
    assert result["models"]["admet"]["count"] == 4
    assert result["models"]["unimol_weights"]["exists"] is True
    assert result["imports"]["unimol_tools"] == "ok"
    assert result["binaries"]["smina"]["exists"] is True
    assert result["binaries"]["obabel"]["exists"] is True


@pytest.mark.integration_assets
@pytest.mark.integration_docking
def test_runtime_doctor_strict_probe_scores_all_four_layers(monkeypatch):
    runtime_doctor = _load_runtime_doctor()
    probe = {
        "smiles": "CCO",
        "statuses": {"l1": "ok", "l2": "ok", "l3": "ok", "l4": "ok"},
        "overall": "ok",
    }
    monkeypatch.setattr(runtime_doctor, "run_backend_probe", lambda _path: probe)

    result = runtime_doctor.inspect_runtime(SCORING_DIR, probe_backends=True)

    assert result["backend_probe"] == probe
    assert result["overall"] == "ok"


def test_runtime_doctor_default_bin_dir_is_package_relative(tmp_path, monkeypatch):
    runtime_doctor = _load_runtime_doctor()
    scoring_dir = tmp_path / "repo" / "scoring"
    scoring_dir.mkdir(parents=True)
    monkeypatch.delenv("DOCK_BIN_DIR", raising=False)

    result = runtime_doctor.inspect_runtime(scoring_dir=scoring_dir)

    assert result["binaries"]["smina"]["path"] == str(scoring_dir.parent / "bin" / "smina")


def test_runtime_doctor_marks_non_executable_binaries_unavailable(tmp_path, monkeypatch):
    runtime_doctor = _load_runtime_doctor()
    scoring_dir = tmp_path / "repo" / "scoring"
    bin_dir = scoring_dir.parent / "bin"
    scoring_dir.mkdir(parents=True)
    bin_dir.mkdir()
    for name in ("smina", "obabel"):
        binary = bin_dir / name
        binary.write_text("not executable", encoding="utf-8")
        binary.chmod(0o644)
    monkeypatch.delenv("DOCK_BIN_DIR", raising=False)

    result = runtime_doctor.inspect_runtime(scoring_dir=scoring_dir)

    assert result["binaries"]["smina"]["executable"] is False
    assert result["binaries"]["obabel"]["executable"] is False
    assert result["overall"] == "failed"


def _import_scoring_modules():
    sys.path.insert(0, str(SCORING_DIR))
    import pipeline_router
    import scoring

    return pipeline_router, scoring


def test_asset_root_resolves_external_model_directories(tmp_path):
    asset_paths = _load_module_from_file("four_level_asset_paths_injected", SCORING_DIR / "asset_paths.py")

    model_root = tmp_path / "scoring" / "models"
    l2_dir = model_root / "bindingdb_l2"
    l2_dir.mkdir(parents=True)
    compat_model = l2_dir / "l2_model_sklearn_1_7_2.joblib"
    compat_model.write_bytes(b"compat")
    (l2_dir / "l2_model.joblib").write_bytes(b"fallback")
    (tmp_path / "ASSET_MANIFEST.json").write_text("{}\n", encoding="utf-8")

    paths = asset_paths.resolve_asset_paths(tmp_path)

    assert paths is not None
    assert paths.root == tmp_path.resolve()
    assert paths.manifest == (tmp_path / "ASSET_MANIFEST.json").resolve()
    assert paths.model_root == model_root.resolve()
    assert paths.l2_model == compat_model.resolve()
    assert paths.admet_model_dir == (model_root / "admet").resolve()
    assert paths.unimol_model_dir == model_root.resolve()
    assert paths.receptor_root == (tmp_path / "scoring" / "receptors").resolve()


def test_router_resolves_relative_receptor_from_external_asset_root(monkeypatch, tmp_path):
    pipeline_router, _ = _import_scoring_modules()
    registry = tmp_path / "source" / "scoring" / "receptor_registry.json"
    registry.parent.mkdir(parents=True)
    registry.write_text(
        json.dumps(
            {
                "entries": {
                    "CHEMBL2051": {
                        "target_name": "Neuraminidase",
                        "pdbqt": "receptors/test.pdbqt",
                        "box_center": [1, 2, 3],
                        "box_size": [20, 20, 20],
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    receptor = tmp_path / "assets" / "scoring" / "receptors" / "test.pdbqt"
    receptor.parent.mkdir(parents=True)
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    monkeypatch.setattr(pipeline_router, "REGISTRY_PATH", registry)
    monkeypatch.setenv("FOUR_LEVEL_ASSET_ROOT", str(tmp_path / "assets"))

    decision = pipeline_router.route("CHEMBL2051", chembl_id="CHEMBL2051")

    assert decision["branch"] == "cascade"
    assert decision["receptor"] == str(receptor.resolve())


def test_molecule_scorer_threads_external_asset_paths(monkeypatch, tmp_path):
    scoring = _load_module_from_file("four_level_scoring_asset_injected", SCORING_DIR / "scoring.py")
    import l2_bindingdb

    captured = {}

    class FakeLayer3:
        def __init__(self, **kwargs):
            captured["l3"] = kwargs

    class FakeLayer2:
        model_kind = "fake"

        def __init__(self, **kwargs):
            captured["l2"] = kwargs

    monkeypatch.setattr(scoring, "Layer3Scorer", FakeLayer3)
    monkeypatch.setattr(l2_bindingdb, "Layer2BindingDB", FakeLayer2)

    scorer = scoring.MoleculeScorer(asset_root=tmp_path, use_unimol=False)

    assert captured["l3"]["model_dir"] == tmp_path / "scoring" / "models" / "admet"
    assert captured["l2"]["model_path"] == str(
        tmp_path / "scoring" / "models" / "bindingdb_l2" / "l2_model.joblib"
    )


def test_unimol_uses_injected_model_directory(tmp_path, monkeypatch):
    import pickle

    unimol_scorer = _load_module_from_file(
        "four_level_unimol_scorer_injected",
        SCORING_DIR / "scripts" / "unimol_scorer.py",
    )

    model_dir = tmp_path / "scoring" / "models"
    model_dir.mkdir(parents=True)
    embeddings = model_dir / "ref_embeddings.npz"
    smiles = model_dir / "ref_smiles.pkl"
    np.savez(embeddings, embeddings=np.ones((2, 4), dtype=np.float32))
    smiles.write_bytes(pickle.dumps(["CCO", "CCN"]))
    integrity = _load_asset_integrity()
    manifest = {
        "schema": "test",
        "files": {
            "scoring/models/ref_embeddings.npz": integrity.sha256_file(embeddings),
            "scoring/models/ref_smiles.pkl": integrity.sha256_file(smiles),
        },
    }
    (tmp_path / "ASSET_MANIFEST.json").write_text(json.dumps(manifest), encoding="utf-8")
    monkeypatch.setenv("FOUR_LEVEL_ASSET_MANIFEST", str(tmp_path / "ASSET_MANIFEST.json"))

    scorer = unimol_scorer.UniMolScorer(device="cpu", model_dir=model_dir)

    assert scorer.model_dir == model_dir.resolve()
    assert scorer.ref_embeddings.shape == (2, 4)
    assert scorer.ref_names == ["FDA_1", "FDA_2"]


def test_scoring_help_exposes_asset_root():
    result = subprocess.run(
        [sys.executable, str(SCORING_DIR / "scoring.py"), "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert result.returncode == 0
    assert "--asset-root" in result.stdout


def _write_base_scores(path: Path, rows: list[tuple[str, str]]) -> Path:
    frame = pd.DataFrame([
        {
            "id": mol_id,
            "smiles": smiles,
            "layer1_status": "ok",
            "layer2_status": "ok",
            "layer3_status": "ok",
            "layer4_status": "ok",
            "docking_normalized": 0.2 + index * 0.1,
            "final_score": 0.4 + index * 0.1,
        }
        for index, (mol_id, smiles) in enumerate(rows)
    ])
    frame.to_csv(path, index=False)
    return path


def test_read_base_scores_requires_exact_input_identity(tmp_path):
    _, scoring = _import_scoring_modules()
    expected = [("mol-1", "CCO"), ("mol-2", "CCN")]
    base = _write_base_scores(tmp_path / "base.csv", expected[:1])

    with pytest.raises(ValueError, match="identity mismatch"):
        scoring.read_base_scores_csv(base, expected)


def test_read_base_scores_rejects_smiles_mutation(tmp_path):
    _, scoring = _import_scoring_modules()
    base = _write_base_scores(tmp_path / "base.csv", [("mol-1", "CCC")])

    with pytest.raises(ValueError, match="identity mismatch"):
        scoring.read_base_scores_csv(base, [("mol-1", "CCO")])


def test_validate_base_columns_rejects_docking_mutation():
    _, scoring = _import_scoring_modules()
    before = [{"id": "mol-1", "smiles": "CCO", "final_score": 0.4}]
    after = [{
        "id": "mol-1",
        "smiles": "CCO",
        "final_score": 0.5,
        "final_score_dock": 0.6,
    }]

    with pytest.raises(ValueError, match="base score mutated"):
        scoring.validate_base_columns_unchanged(before, after)


def test_load_or_score_results_skips_scorer_for_base_scores(tmp_path, monkeypatch):
    from types import SimpleNamespace

    _, scoring = _import_scoring_modules()
    base = _write_base_scores(tmp_path / "base.csv", [("mol-1", "CCO")])

    class ForbiddenScorer:
        def __init__(self, **kwargs):
            del kwargs
            raise AssertionError("scorer initialized")

    monkeypatch.setattr(scoring, "MoleculeScorer", ForbiddenScorer)
    args = SimpleNamespace(
        base_scores=str(base),
        target="CHEMBL2051",
        target_seq=None,
        l2="bindingdb",
        strict_backends=True,
        asset_root=None,
    )

    rows, snapshot = scoring.load_or_score_results(
        args,
        [("mol-1", "CCO")],
        {"l2_model_path": None},
        "Neuraminidase",
    )

    assert rows == snapshot
    assert rows[0]["final_score"] == 0.4


def test_scoring_help_exposes_base_scores():
    result = subprocess.run(
        [sys.executable, str(SCORING_DIR / "scoring.py"), "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert result.returncode == 0
    assert "--base-scores" in result.stdout


def test_open_molecule_bridge_counts_complete_structure_docking_success(tmp_path):
    bridge = _load_module_from_file(
        "open_molecule_python_bridge",
        ROOT / "apps" / "open-molecule-lab" / "server" / "python-bridge.py",
    )
    result_path = tmp_path / "scores.csv"
    pd.DataFrame(
        {
            "id": ["mol-1", "mol-2", "mol-3"],
            "smiles": ["CCO", "CCN", "CCC"],
            "layer1_status": ["ok", "ok", "ok"],
            "layer2_status": ["ok", "ok", "ok"],
            "layer3_status": ["ok", "ok", "ok"],
            "layer4_status": ["ok", "ok", "ok"],
            "final_score": [0.7, 0.6, 0.5],
            "structure_docking_status": ["ok", "prep_failed", "ok"],
            "final_score_dock": [0.9, 0.6, 0.8],
        }
    ).to_csv(result_path, index=False)

    summary = bridge.summarize_results(result_path)

    assert summary["structureDockingOk"] == 2
    assert summary["rankingScoreField"] == "final_score_dock"


def test_molecule_scorer_defaults_to_packaged_bindingdb_backend():
    _, scoring = _import_scoring_modules()

    default = inspect.signature(scoring.MoleculeScorer).parameters["l2_method"].default

    assert default == "bindingdb"


def test_router_imports_under_system_python():
    code = f"import sys; sys.path.insert(0, {str(SCORING_DIR)!r}); import pipeline_router"
    result = subprocess.run(
        ["/usr/bin/python3", "-c", code],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr


def test_router_rejects_registered_receptor_with_missing_asset(monkeypatch, tmp_path):
    pipeline_router, _ = _import_scoring_modules()
    registry = tmp_path / "receptors.json"
    registry.write_text(
        json.dumps(
            {
                "entries": {
                    "CHEMBL2051": {
                        "target_name": "Neuraminidase",
                        "pdbqt": str(tmp_path / "missing.pdbqt"),
                        "box_center": [1, 2, 3],
                        "box_size": [20, 20, 20],
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(pipeline_router, "REGISTRY_PATH", registry)

    decision = pipeline_router.route("CHEMBL2051", chembl_id="CHEMBL2051")

    assert decision["branch"] == "library"
    assert "受体资产无效" in decision["rationale"]


def test_forced_cascade_rejects_missing_receptor(monkeypatch, tmp_path):
    pipeline_router, _ = _import_scoring_modules()
    registry = tmp_path / "receptors.json"
    registry.write_text(
        json.dumps(
            {
                "entries": {
                    "CHEMBL2051": {
                        "target_name": "Neuraminidase",
                        "pdbqt": str(tmp_path / "missing.pdbqt"),
                        "box_center": [1, 2, 3],
                        "box_size": [20, 20, 20],
                    }
                }
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(pipeline_router, "REGISTRY_PATH", registry)

    with pytest.raises(RuntimeError, match="cascade|级联|受体"):
        pipeline_router.route("CHEMBL2051", chembl_id="CHEMBL2051", mode="cascade")


def test_router_rationale_depends_on_assets_not_benchmark_auc():
    pipeline_router, _ = _import_scoring_modules()

    for decision in (
        pipeline_router.route("CHEMBL2051", chembl_id="CHEMBL2051"),
        pipeline_router.route("CHEMBL999999", chembl_id="CHEMBL999999"),
    ):
        assert "AUC" not in decision["rationale"]
        assert "0.935" not in decision["rationale"]


def test_strict_l4_rejects_missing_backend(monkeypatch):
    _, scoring = _import_scoring_modules()
    real_import = builtins.__import__

    class FakeLayer3:
        def __init__(self, **kwargs):
            del kwargs

    monkeypatch.setattr(scoring, "Layer3Scorer", FakeLayer3)

    def fail_unimol_import(name, *args, **kwargs):
        if name == "scripts.unimol_scorer":
            raise ModuleNotFoundError("unimol backend unavailable")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fail_unimol_import)
    with pytest.raises(RuntimeError, match="L4"):
        scoring.MoleculeScorer(use_unimol=True, strict_backends=True)


def test_strict_l2_rejects_missing_backend(monkeypatch):
    _, scoring = _import_scoring_modules()
    import l2_bindingdb

    class FakeLayer3:
        def __init__(self, **kwargs):
            del kwargs

    monkeypatch.setattr(scoring, "Layer3Scorer", FakeLayer3)

    monkeypatch.setattr(
        l2_bindingdb.Layer2BindingDB,
        "_ensure_model",
        lambda self: (_ for _ in ()).throw(FileNotFoundError("missing L2")),
    )

    with pytest.raises(RuntimeError, match="L2"):
        scoring.MoleculeScorer(use_unimol=False, l2_method="bindingdb", strict_backends=True)


def test_non_strict_missing_l2_reports_explicit_unavailable_backend(monkeypatch):
    _, scoring = _import_scoring_modules()
    import l2_bindingdb

    monkeypatch.setattr(
        l2_bindingdb.Layer2BindingDB,
        "_ensure_model",
        lambda self: (_ for _ in ()).throw(FileNotFoundError("missing L2")),
    )

    class FakeDeepPurpose:
        def __init__(self, **kwargs):
            del kwargs
            raise AssertionError("formal bindingdb fallback must not construct DeepPurpose")

    monkeypatch.setattr(scoring, "Layer2DeepPurpose", FakeDeepPurpose)
    scorer = scoring.MoleculeScorer(
        use_unimol=False,
        l2_method="bindingdb",
        default_target_text="target",
        strict_backends=False,
    )

    result = scorer.score_one("CCO")

    assert scorer.l2_method == "bindingdb"
    assert result["layer2_status"].startswith("failed:backend_unavailable")
    assert result["layer2_backend"] == "BindingDB-L2-unavailable"
    assert result["docking_normalized"] == 0.0


@pytest.mark.parametrize("value", ["-0.1", "1.1", "nan", "inf"])
def test_dock_fusion_parser_rejects_values_outside_unit_interval(value):
    _, scoring = _import_scoring_modules()

    with pytest.raises(Exception, match="0.*1|finite|有限|范围"):
        scoring._parse_unit_interval(value)

    assert scoring._parse_unit_interval("0") == 0.0
    assert scoring._parse_unit_interval("1") == 1.0


def test_ranking_displays_fused_score_when_results_are_fusion_sorted(capsys):
    _, scoring = _import_scoring_modules()
    results = [{
        "id": "mol-x",
        "final_score": 0.1111,
        "final_score_dock": 0.9876,
        "layer1_score": 0.2,
        "docking_normalized": 0.3,
        "admet_score": 0.4,
        "unimol_score": 0.5,
        "gate_status": "PASS",
        "docking_method": "test",
    }]

    scoring.MoleculeScorer.print_ranking(results)

    output = capsys.readouterr().out
    assert "级联分" in output
    assert "0.9876" in output
    assert "0.1111" not in output


def test_bindingdb_dispatch_never_substitutes_molecule_id_for_target_text():
    _, scoring = _import_scoring_modules()
    call = {}

    class FakeL1:
        def score(self, smiles, all_smiles=None):
            del smiles, all_smiles
            return {"layer1_score": 0.8, "valid": 1}

    class FakeL2:
        model_path = Path("l2.joblib")

        def score(self, smiles, target_text="", mol_id="mol"):
            call.update(smiles=smiles, target_text=target_text, mol_id=mol_id)
            return {
                "docking_normalized": 0.6,
                "docking_score_kcal_mol": "",
                "docking_status": "success",
                "docking_method": "BindingDB-L2",
            }

    class FakeL3:
        def score(self, smiles):
            del smiles
            return {
                "admet_score": 0.7,
                "toxicity_count": 0,
                "toxicity_severe": 0,
                "toxicity_flags": "low",
                "bbb_prob": 0.4,
            }

    scorer = object.__new__(scoring.MoleculeScorer)
    scorer.l1, scorer.l2, scorer.l3 = FakeL1(), FakeL2(), FakeL3()
    scorer.l4, scorer.unimol = scoring.Layer4Aggregator(), None
    scorer.l2_method = "bindingdb"
    scorer.default_target_text = None
    scorer.strict_backends = False

    scorer.score_one("CCO", mol_id="candidate-17")

    assert call == {"smiles": "CCO", "target_text": "", "mol_id": "candidate-17"}


def test_strict_score_one_rejects_non_ok_l2_status():
    _, scoring = _import_scoring_modules()

    class FakeL1:
        def score(self, smiles, all_smiles=None):
            del smiles, all_smiles
            return {"layer1_score": 0.8, "valid": 1}

    class FakeL2:
        model_path = Path("missing.joblib")

        def score(self, smiles, target_text="", mol_id="mol"):
            del smiles, target_text, mol_id
            return {
                "docking_normalized": 0.0,
                "docking_score_kcal_mol": "",
                "docking_status": "failed:missing",
                "docking_method": "BindingDB-L2",
            }

    class FakeL3:
        def score(self, smiles):
            del smiles
            return {"admet_score": 0.7, "toxicity_count": 0, "toxicity_severe": 0, "bbb_prob": 0.4}

    class FakeL4:
        def score(self, smiles):
            del smiles
            return {"unimol_score": 0.5, "pos_similarity": 0.5, "neg_similarity": 0.0}

    scorer = object.__new__(scoring.MoleculeScorer)
    scorer.l1, scorer.l2, scorer.l3 = FakeL1(), FakeL2(), FakeL3()
    scorer.l4, scorer.unimol = scoring.Layer4Aggregator(), FakeL4()
    scorer.l2_method = "bindingdb"
    scorer.default_target_text = "target"
    scorer.strict_backends = True

    with pytest.raises(RuntimeError, match="L2"):
        scorer.score_one("CCO")


def test_strict_l3_requires_all_model_assets(tmp_path):
    _, scoring = _import_scoring_modules()

    with pytest.raises(RuntimeError, match="L3"):
        scoring.Layer3Scorer(strict_backends=True, model_dir=tmp_path)


def test_shared_four_level_combiner_matches_cli_formula():
    _, scoring = _import_scoring_modules()
    l1 = np.array([0.2, 0.8])
    l2 = np.array([0.4, 0.6])
    l3 = np.array([0.5, 0.7])
    l4 = np.array([0.9, 0.1])

    result = scoring.Layer4Aggregator.combine(l1, l2, l3, l4)

    expected = np.round(0.20 * l1 + 0.50 * l2 + 0.20 * l3 + 0.10 * l4, 4)
    assert np.array_equal(result, expected)


def test_cli_benchmark_is_not_evaluated_without_explicit_labels():
    _, scoring = _import_scoring_modules()
    rows = [
        {"id": "a", "docking_score_kcal_mol": "", "final_score": 0.9},
        {"id": "b", "docking_score_kcal_mol": "", "final_score": 0.1},
    ]

    result = scoring.Layer4Aggregator.vs_benchmark(rows)

    assert result["status"] == "not_evaluated"
    assert result["positive"] == 0
    assert result["negative"] == 0


def test_cli_benchmark_uses_only_explicit_positive_ids():
    _, scoring = _import_scoring_modules()
    rows = [
        {"id": "active", "docking_score_kcal_mol": -8.0, "final_score": 0.9},
        {"id": "decoy", "docking_score_kcal_mol": -5.0, "final_score": 0.1},
    ]

    result = scoring.Layer4Aggregator.vs_benchmark(rows, positive_ids={"active"})

    assert result["status"] == "evaluated"
    assert result["positive"] == 1
    assert result["negative"] == 1


def test_cli_benchmark_uses_higher_binding_affinity_and_normalized_bedroc():
    _, scoring = _import_scoring_modules()
    rows = [
        {"id": "active", "docking_score_kcal_mol": 8.0, "final_score": 0.9},
        {"id": "decoy", "docking_score_kcal_mol": 5.0, "final_score": 0.1},
    ]

    result = scoring.Layer4Aggregator.vs_benchmark(rows, positive_ids={"active"})

    assert result["roc_auc"] == 1.0
    assert result["ef_10pct"] == 2.0
    assert result["bedroc_alpha20"] == 1.0


def test_enrichment_factor_counts_tied_positive_and_negative_scores_once_each():
    _, scoring = _import_scoring_modules()

    result = scoring.Layer4Aggregator.enrichment_factor(
        [8.0], [8.0, 5.0], frac=2 / 3, higher_is_better=True
    )

    assert result == 1.5


def test_score_one_emits_auditable_status_for_every_layer():
    _, scoring = _import_scoring_modules()

    class FakeL1:
        def score(self, smiles, all_smiles=None):
            return {"layer1_score": 0.8}

    class FakeL2:
        model_path = Path("models/bindingdb_l2/l2_model.joblib")

        def score(self, smiles, target_text="", mol_id="mol"):
            return {
                "docking_normalized": 0.6,
                "docking_score_kcal_mol": "",
                "docking_status": "success",
                "docking_method": "BindingDB-L2-mlp",
            }

    class FakeL3:
        def score(self, smiles):
            return {"admet_score": 0.7, "toxicity_count": 0, "toxicity_severe": 0, "bbb_prob": 0.4}

    class FakeL4:
        def score(self, smiles):
            return {"unimol_score": 0.5, "pos_similarity": 0.5, "neg_similarity": 0.0}

    scorer = object.__new__(scoring.MoleculeScorer)
    scorer.l1 = FakeL1()
    scorer.l2 = FakeL2()
    scorer.l3 = FakeL3()
    scorer.l4 = scoring.Layer4Aggregator()
    scorer.unimol = FakeL4()
    scorer.l2_method = "bindingdb"
    scorer.default_target_text = "CHEMBL1 target"
    scorer.strict_backends = True

    result = scorer.score_one("CCO", mol_id="ethanol")

    for layer in range(1, 5):
        assert result[f"layer{layer}_status"] == "ok"
        assert result[f"layer{layer}_backend"]
        assert result[f"layer{layer}_model_asset_id"]
    assert result["bbb_prob"] == 0.4
    assert "bbb_logBB" not in result


def test_save_csv_creates_missing_parent(tmp_path):
    _, scoring = _import_scoring_modules()
    destination = tmp_path / "nested" / "scores.csv"

    scoring.MoleculeScorer.save_csv([{"id": "x", "bbb_prob": 0.4}], str(destination))

    assert destination.is_file()


def test_read_molecules_csv_enforces_required_values_and_unique_ids(tmp_path):
    _, scoring = _import_scoring_modules()
    cases = {
        "missing_header.csv": "name,smiles\na,CCO\n",
        "blank_value.csv": "id,smiles\na,\n",
        "duplicate_id.csv": "id,smiles\na,CCO\na,CCN\n",
    }

    for filename, content in cases.items():
        path = tmp_path / filename
        path.write_text(content, encoding="utf-8")
        with pytest.raises(ValueError):
            scoring.read_molecules_csv(path)


def test_read_molecules_csv_accepts_excel_utf8_bom(tmp_path):
    _, scoring = _import_scoring_modules()
    path = tmp_path / "excel.csv"
    path.write_bytes("\ufeffid,smiles\nmol-1,CCO\n".encode("utf-8"))

    assert scoring.read_molecules_csv(path) == [("mol-1", "CCO")]


def test_batch_cli_rejects_non_positive_numeric_arguments():
    from scientific_validation.four_level_cli_1kx10k import batch_cli

    invalid = (
        ["prepare", "--targets", "0"],
        ["prepare", "--pool-size", "0"],
        ["cache-layers", "--batch-size", "0"],
        ["cache-layers", "--shard-size", "-1"],
        ["dock", "--top-n", "0"],
        ["dock", "--workers", "0"],
        ["dock", "--timeout", "0"],
        ["dock", "--cpu", "0"],
    )
    for argv in invalid:
        with pytest.raises(SystemExit) as exc_info:
            batch_cli._build_parser().parse_args(argv)
        assert exc_info.value.code == 2


def test_find_binary_accepts_direct_executable_override(tmp_path, monkeypatch):
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    binary = tmp_path / "custom-smina"
    binary.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    binary.chmod(0o755)
    monkeypatch.setenv("SMINA_BIN", str(binary))
    monkeypatch.delenv("DOCK_BIN_DIR", raising=False)

    assert dock_rerank.find_binary("smina") == str(binary)


def test_prep_ligand_rejects_stale_output_after_obabel_failure(tmp_path, monkeypatch):
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    output = tmp_path / "ligand.pdbqt"
    output.write_text("stale ligand", encoding="utf-8")

    class FailedProcess:
        returncode = 1

    monkeypatch.setattr(dock_rerank.subprocess, "run", lambda *args, **kwargs: FailedProcess())

    assert dock_rerank.prep_ligand("CCO", str(output), "obabel") is False
    assert not output.exists()


def test_protonated_receptor_rejects_failure_without_using_adjacent_stale_file(tmp_path, monkeypatch):
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    stale = Path(str(receptor) + ".protonated.pdb")
    stale.write_text("STALE\n", encoding="utf-8")

    class FailedProcess:
        returncode = 1
        stderr = "conversion failed"
        stdout = ""

    monkeypatch.setattr(dock_rerank.subprocess, "run", lambda *args, **kwargs: FailedProcess())

    with pytest.raises(RuntimeError, match="conversion failed"):
        with dock_rerank.protonated_receptor(receptor, "obabel"):
            raise AssertionError("failed protonation must not yield a receptor")
    assert stale.read_text(encoding="utf-8") == "STALE\n"


def test_protonated_receptor_cleans_temporary_output(tmp_path, monkeypatch):
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    receptor = tmp_path / "receptor.pdb"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    created = {}

    class SuccessfulProcess:
        returncode = 0
        stderr = ""
        stdout = ""

    def successful_run(command, **kwargs):
        output = Path(command[command.index("-O") + 1])
        output.write_text("PROTONATED\n", encoding="utf-8")
        created["path"] = output
        return SuccessfulProcess()

    monkeypatch.setattr(dock_rerank.subprocess, "run", successful_run)

    with dock_rerank.protonated_receptor(receptor, "obabel") as protonated:
        assert Path(protonated).is_file()
    assert not created["path"].exists()


def test_dock_one_keeps_untrusted_molecule_id_inside_workdir(tmp_path, monkeypatch):
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    workdir = tmp_path / "work"
    workdir.mkdir()
    captured = {}
    monkeypatch.setattr(
        dock_rerank,
        "prep_ligand",
        lambda smiles, path, *args, **kwargs: captured.update(path=path) or False,
    )
    reranker = object.__new__(dock_rerank.DockingReranker)
    reranker.workdir = str(workdir)
    reranker.obabel_bin = "obabel"
    reranker.smina_bin = "smina"
    reranker.receptor = "receptor"
    reranker.center = (1, 2, 3)
    reranker.size = (20, 20, 20)
    reranker.pH = 7.4
    reranker.exhaustiveness = 4
    reranker.cpu = 1
    reranker.num_modes = 1
    reranker.seed = 42
    reranker.flexres = ""

    reranker.dock_one("../../outside", "CCO")

    assert Path(captured["path"]).resolve().is_relative_to(workdir.resolve())


def test_docking_request_requires_receptor_and_valid_manual_box(tmp_path):
    _, scoring = _import_scoring_modules()
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")

    with pytest.raises(ValueError, match="receptor|受体"):
        scoring.validate_docking_request(
            receptor=None,
            box_center=None,
            box_size=None,
            dock_rerank=True,
            manual_receptor=False,
        )
    with pytest.raises(ValueError, match="box|对接盒"):
        scoring.validate_docking_request(
            receptor=str(receptor),
            box_center=None,
            box_size=None,
            dock_rerank=False,
            manual_receptor=True,
        )
    with pytest.raises(ValueError, match="三个|3"):
        scoring.validate_docking_request(
            receptor=str(receptor),
            box_center="1,2",
            box_size="20,20,20",
            dock_rerank=True,
            manual_receptor=True,
        )

    center, size = scoring.validate_docking_request(
        receptor=str(receptor),
        box_center="1,2,3",
        box_size="20,21,22",
        dock_rerank=True,
        manual_receptor=True,
    )
    assert center == (1.0, 2.0, 3.0)
    assert size == (20.0, 21.0, 22.0)


def test_forced_cascade_propagates_docking_backend_failure(tmp_path, monkeypatch):
    pipeline_router, scoring = _import_scoring_modules()
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("id,smiles\nmol-1,CCO\n", encoding="utf-8")

    monkeypatch.setattr(
        pipeline_router,
        "route",
        lambda *args, **kwargs: {
            "branch": "cascade",
            "target_text": "target",
            "method": "test",
            "receptor": str(receptor),
            "box_center": [1, 2, 3],
            "box_size": [20, 20, 20],
            "l2_model_path": None,
            "rationale": "test cascade",
        },
    )

    class FakeScorer:
        def __init__(self, **kwargs):
            del kwargs

        def score_batch(self, rows):
            mid, smiles = rows[0]
            return [{
                "id": mid,
                "smiles": smiles,
                "layer1_score": 0.8,
                "docking_normalized": 0.6,
                "admet_score": 0.7,
                "unimol_score": 0.5,
                "final_score": 0.67,
                "gate_status": "PASS",
                "docking_method": "fake",
            }]

        save_csv = staticmethod(lambda results, path: None)
        print_ranking = staticmethod(lambda results: None)

    class FailingReranker:
        def __init__(self, **kwargs):
            del kwargs
            raise RuntimeError("docking binary unavailable")

    monkeypatch.setattr(scoring, "MoleculeScorer", FakeScorer)
    monkeypatch.setattr(dock_rerank, "DockingReranker", FailingReranker)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scoring.py",
            "--input", str(candidates),
            "--output", str(tmp_path / "scores.csv"),
            "--target", "CHEMBL2051",
            "--mode", "cascade",
        ],
    )

    with pytest.raises(RuntimeError, match="docking binary unavailable"):
        scoring.main()


def test_explicit_receptor_and_box_take_precedence_over_registry_route(tmp_path, monkeypatch):
    pipeline_router, scoring = _import_scoring_modules()
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    explicit_receptor = tmp_path / "explicit.pdbqt"
    registry_receptor = tmp_path / "registry.pdbqt"
    explicit_receptor.write_text("EXPLICIT\n", encoding="utf-8")
    registry_receptor.write_text("REGISTRY\n", encoding="utf-8")
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("id,smiles\nmol-1,CCO\n", encoding="utf-8")
    captured = {}

    monkeypatch.setattr(
        pipeline_router,
        "route",
        lambda *args, **kwargs: {
            "branch": "cascade",
            "target_text": "target",
            "method": "test",
            "receptor": str(registry_receptor),
            "box_center": [9, 9, 9],
            "box_size": [30, 30, 30],
            "l2_model_path": None,
            "rationale": "registered cascade",
        },
    )

    class FakeScorer:
        def __init__(self, **kwargs):
            del kwargs

        def score_batch(self, rows):
            mid, smiles = rows[0]
            return [{
                "id": mid, "smiles": smiles, "layer1_score": 0.8,
                "docking_normalized": 0.6, "admet_score": 0.7,
                "unimol_score": 0.5, "final_score": 0.67,
                "gate_status": "PASS", "docking_method": "fake",
            }]

        save_csv = staticmethod(lambda results, path: None)
        print_ranking = staticmethod(lambda results: None)

    class CapturingReranker:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def dock_all(self, rows, mode="le"):
            del mode
            mid, smiles = rows[0]
            return [{
                "id": mid, "smiles": smiles, "affinity": -7.0,
                "heavy_atoms": 3, "ligand_efficiency": -2.3333,
                "dock_rerank_rank": 1, "status": "ok",
            }]

    monkeypatch.setattr(scoring, "MoleculeScorer", FakeScorer)
    monkeypatch.setattr(dock_rerank, "DockingReranker", CapturingReranker)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scoring.py", "--input", str(candidates),
            "--output", str(tmp_path / "scores.csv"),
            "--target", "CHEMBL2051", "--mode", "cascade",
            "--receptor", str(explicit_receptor),
            "--box-center", "1,2,3", "--box-size", "20,21,22",
        ],
    )

    scoring.main()

    assert captured["receptor"] == str(explicit_receptor)
    assert captured["center"] == (1.0, 2.0, 3.0)
    assert captured["size"] == (20.0, 21.0, 22.0)


def test_forced_cascade_rejects_zero_successful_docking_results(tmp_path, monkeypatch):
    pipeline_router, scoring = _import_scoring_modules()
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("id,smiles\nmol-1,CCO\n", encoding="utf-8")
    monkeypatch.setattr(
        pipeline_router,
        "route",
        lambda *args, **kwargs: {
            "branch": "cascade", "target_text": "target", "method": "test",
            "receptor": str(receptor), "box_center": [1, 2, 3],
            "box_size": [20, 20, 20], "l2_model_path": None,
            "rationale": "test cascade",
        },
    )

    class FakeScorer:
        def __init__(self, **kwargs):
            del kwargs

        def score_batch(self, rows):
            mid, smiles = rows[0]
            return [{
                "id": mid, "smiles": smiles, "layer1_score": 0.8,
                "docking_normalized": 0.6, "admet_score": 0.7,
                "unimol_score": 0.5, "final_score": 0.67,
                "gate_status": "PASS", "docking_method": "fake",
            }]

        save_csv = staticmethod(lambda results, path: None)
        print_ranking = staticmethod(lambda results: None)

    class ZeroSuccessReranker:
        def __init__(self, **kwargs):
            del kwargs

        def dock_all(self, rows, mode="le"):
            del mode
            mid, smiles = rows[0]
            return [{
                "id": mid, "smiles": smiles, "affinity": None,
                "heavy_atoms": 3, "ligand_efficiency": None,
                "dock_rerank_rank": None, "status": "dock_no_score",
            }]

    monkeypatch.setattr(scoring, "MoleculeScorer", FakeScorer)
    monkeypatch.setattr(dock_rerank, "DockingReranker", ZeroSuccessReranker)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scoring.py", "--input", str(candidates),
            "--output", str(tmp_path / "scores.csv"),
            "--target", "CHEMBL2051", "--mode", "cascade",
        ],
    )

    with pytest.raises(RuntimeError, match="zero|0|零|没有.*成功|无.*成功"):
        scoring.main()


def test_explicit_zero_dock_fusion_disables_cascade_score_fusion(tmp_path, monkeypatch):
    pipeline_router, scoring = _import_scoring_modules()
    sys.path.insert(0, str(SCORING_DIR))
    import dock_rerank

    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("id,smiles\nmol-1,CCO\n", encoding="utf-8")
    saved = {}

    monkeypatch.setattr(
        pipeline_router,
        "route",
        lambda *args, **kwargs: {
            "branch": "cascade", "target_text": "target", "method": "test",
            "receptor": str(receptor), "box_center": [1, 2, 3],
            "box_size": [20, 20, 20], "l2_model_path": None,
            "rationale": "test cascade",
        },
    )

    class FakeScorer:
        def __init__(self, **kwargs):
            del kwargs

        def score_batch(self, rows):
            mid, smiles = rows[0]
            return [{
                "id": mid, "smiles": smiles, "layer1_score": 0.8,
                "docking_normalized": 0.6, "admet_score": 0.7,
                "unimol_score": 0.5, "final_score": 0.67,
                "gate_status": "PASS", "docking_method": "fake",
            }]

        @staticmethod
        def save_csv(results, path):
            saved[path] = [dict(row) for row in results]

        print_ranking = staticmethod(lambda results: None)

    class SuccessfulReranker:
        def __init__(self, **kwargs):
            del kwargs

        def dock_all(self, rows, mode="le"):
            del mode
            mid, smiles = rows[0]
            return [{
                "id": mid, "smiles": smiles, "affinity": -7.0,
                "heavy_atoms": 3, "ligand_efficiency": -2.3333,
                "dock_rerank_rank": 1, "status": "ok",
            }]

    monkeypatch.setattr(scoring, "MoleculeScorer", FakeScorer)
    monkeypatch.setattr(dock_rerank, "DockingReranker", SuccessfulReranker)
    monkeypatch.setattr(
        dock_rerank,
        "cascade_corrected_fusion",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("explicit zero must not invoke score fusion")
        ),
    )
    output = tmp_path / "scores.csv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scoring.py", "--input", str(candidates), "--output", str(output),
            "--target", "CHEMBL2051", "--mode", "cascade",
            "--dock-fusion", "0",
        ],
    )

    scoring.main()

    assert "final_score_dock" not in saved[str(output)][0]


def test_forced_library_rejects_manual_docking_arguments(tmp_path, monkeypatch, capsys):
    _, scoring = _import_scoring_modules()
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("RECEPTOR\n", encoding="utf-8")
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("id,smiles\nmol-1,CCO\n", encoding="utf-8")

    class UnexpectedScorer:
        def __init__(self, **kwargs):
            del kwargs
            raise AssertionError("library/docking conflict must fail before scoring")

    monkeypatch.setattr(scoring, "MoleculeScorer", UnexpectedScorer)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scoring.py", "--input", str(candidates),
            "--target", "CHEMBL2051", "--mode", "library",
            "--receptor", str(receptor),
            "--box-center", "1,2,3", "--box-size", "20,20,20",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        scoring.main()

    assert exc_info.value.code == 2
    error = capsys.readouterr().err
    assert "library" in error
    assert "receptor" in error or "受体" in error


@pytest.mark.parametrize(
    "extra",
    [["--dock-mode", "raw"], ["--cascade-top-n", "10"]],
)
def test_forced_library_rejects_explicit_docking_controls(tmp_path, monkeypatch, capsys, extra):
    _, scoring = _import_scoring_modules()
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("id,smiles\nmol-1,CCO\n", encoding="utf-8")

    class UnexpectedScorer:
        def __init__(self, **kwargs):
            del kwargs
            raise AssertionError("library/docking control must fail before scoring")

    monkeypatch.setattr(scoring, "MoleculeScorer", UnexpectedScorer)
    monkeypatch.setattr(
        sys,
        "argv",
        ["scoring.py", "--input", str(candidates), "--mode", "library", *extra],
    )

    with pytest.raises(SystemExit) as exc_info:
        scoring.main()

    assert exc_info.value.code == 2
    assert "library" in capsys.readouterr().err


def test_auto_library_rejects_explicit_fusion_before_scoring(tmp_path, monkeypatch, capsys):
    _, scoring = _import_scoring_modules()
    candidates = tmp_path / "candidates.csv"
    candidates.write_text("id,smiles\nmol-1,CCO\n", encoding="utf-8")

    class UnexpectedScorer:
        def __init__(self, **kwargs):
            del kwargs
            raise AssertionError("auto library fusion must fail before scoring")

    monkeypatch.setattr(scoring, "MoleculeScorer", UnexpectedScorer)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scoring.py", "--input", str(candidates), "--target", "CHEMBL999999",
            "--mode", "auto", "--dock-fusion", "0.5",
        ],
    )

    with pytest.raises(SystemExit) as exc_info:
        scoring.main()

    assert exc_info.value.code == 2
    assert "library" in capsys.readouterr().err


def test_importing_scoring_does_not_patch_global_os_directory_functions():
    code = (
        "import os, sys; before=(os.mkdir, os.makedirs); "
        f"sys.path.insert(0, {str(SCORING_DIR)!r}); import scoring; "
        "assert (os.mkdir, os.makedirs) == before"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr


def test_cascade_docks_only_precomputed_l2_top_n():
    _, scoring = _import_scoring_modules()
    rows = [("a", "CC"), ("b", "CCC"), ("c", "CCCC"), ("d", "CCO")]
    scores = {"a": 0.1, "b": 0.9, "c": 0.5, "d": 0.5}

    selected = scoring.select_cascade_candidates(rows, scores, top_n=3)

    assert [mid for mid, _ in selected] == ["b", "c", "d"]


def test_cli_exposes_strict_backends_and_cascade_top_n():
    result = subprocess.run(
        [sys.executable, str(SCORING_DIR / "scoring.py"), "--help"],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert "--strict-backends" in result.stdout
    assert "--cascade-top-n" in result.stdout
    assert "--l2 {bindingdb}" in result.stdout
    assert "{bindingdb,bindingdb_seq,deeppurpose}" not in result.stdout


def test_source_release_builder_preserves_utf8_names_and_excludes_junk(tmp_path):
    module_path = ROOT / "scripts" / "build_source_release.py"
    spec = importlib.util.spec_from_file_location("four_level_build_source_release", module_path)
    assert spec and spec.loader
    build_source_release = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(build_source_release)

    source = tmp_path / "four-level-molecule-cli"
    (source / "scoring").mkdir(parents=True)
    (source / "scoring" / "流感案例.md").write_text("case\n", encoding="utf-8")
    (source / "README.md").write_text("readme\n", encoding="utf-8")
    (source / "legacy").mkdir()
    (source / "legacy" / "research.py").write_text("print('source')\n", encoding="utf-8")
    (source / "legacy" / "report.md").write_text("aggregate report\n", encoding="utf-8")
    (source / "legacy" / "candidate_rows.csv").write_text("id,smiles\n1,CCO\n", encoding="utf-8")
    (source / "legacy" / "results.json").write_text("{}\n", encoding="utf-8")
    (source / "legacy" / "scores.parquet").write_bytes(b"row-level data")
    (source / "__pycache__").mkdir()
    (source / "__pycache__" / "cached.pyc").write_bytes(b"bytecode")
    (source / "logs").mkdir()
    (source / "logs" / "runtime.log").write_text("host path\n", encoding="utf-8")
    (source / "build" / "lib").mkdir(parents=True)
    (source / "build" / "lib" / "generated.py").write_text("generated\n", encoding="utf-8")
    (source / "package.egg-info").mkdir()
    (source / "package.egg-info" / "PKG-INFO").write_text("generated metadata\n", encoding="utf-8")
    first = tmp_path / "first.zip"
    second = tmp_path / "second.zip"

    build_source_release.build_source_zip(source, first)
    build_source_release.build_source_zip(source, second)

    assert first.read_bytes() == second.read_bytes()
    with zipfile.ZipFile(first) as archive:
        names = archive.namelist()
        chinese = "four-level-molecule-cli/scoring/流感案例.md"
        assert chinese in names
        assert archive.getinfo(chinese).flag_bits & 0x800
        assert not any("__pycache__" in name or name.endswith(".pyc") for name in names)
        assert not any("/logs/" in name or name.endswith(".log") for name in names)
        assert not any("/build/" in name or ".egg-info/" in name for name in names)
        assert "four-level-molecule-cli/legacy/research.py" in names
        assert "four-level-molecule-cli/legacy/report.md" in names
        assert "four-level-molecule-cli/legacy/candidate_rows.csv" not in names
        assert "four-level-molecule-cli/legacy/results.json" not in names
        assert "four-level-molecule-cli/legacy/scores.parquet" not in names


def test_scoring_is_importable_as_a_package():
    result = subprocess.run(
        [sys.executable, "-c", "import scoring.scoring as module; assert callable(module.main)"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == 0, result.stderr


def test_asset_release_builder_is_deterministic_and_excludes_appledouble(tmp_path):
    module_path = ROOT / "scripts" / "build_source_release.py"
    spec = importlib.util.spec_from_file_location("four_level_build_asset_release", module_path)
    assert spec and spec.loader
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)

    source = tmp_path / "four-level-molecule-cli-offline-assets"
    (source / "scoring" / "models").mkdir(parents=True)
    (source / "scoring" / "models" / "model.pt").write_bytes(b"weights")
    (source / "scoring" / "models" / "._model.pt").write_bytes(b"appledouble")
    (source / ".DS_Store").write_bytes(b"finder")
    (source / "logs").mkdir()
    (source / "logs" / "runtime.log").write_text("log\n", encoding="utf-8")
    first = tmp_path / "first.tar.gz"
    second = tmp_path / "second.tar.gz"

    builder.build_asset_tar(source, first)
    builder.build_asset_tar(source, second)

    assert first.read_bytes() == second.read_bytes()
    with tarfile.open(first, "r:gz") as archive:
        names = archive.getnames()
        assert "four-level-molecule-cli-offline-assets/scoring/models/model.pt" in names
        assert not any(Path(name).name.startswith("._") for name in names)
        assert not any(Path(name).name == ".DS_Store" for name in names)
        assert not any("/logs/" in name or name.endswith(".log") for name in names)


@pytest.mark.integration_assets
def test_unimol_batch_matches_scalar_and_reports_status(monkeypatch):
    weight_dir = (
        SCORING_DIR
        / "models"
        / "unimol"
        / "models--dptech--Uni-Mol-Models"
        / "snapshots"
        / "9f19c45c718192888a1c8a1c905f69f0755ea502"
    )
    monkeypatch.setenv("UNIMOL_WEIGHT_DIR", str(weight_dir))
    sys.path.insert(0, str(SCORING_DIR))
    for name in ("scripts.unimol_scorer", "scripts.unimol_embedding"):
        sys.modules.pop(name, None)
    module = importlib.import_module("scripts.unimol_scorer")

    smiles = ["CCO", "CC(=O)Oc1ccccc1C(=O)O"]
    scorer = module.UniMolScorer(device="cpu")
    scalar = [scorer.score(smi)["unimol_score"] for smi in smiles]
    scorer._cache.clear()
    batch = scorer.score_many(smiles, batch_size=2)

    assert np.allclose(scalar, [row["unimol_score"] for row in batch], atol=1e-4)
    assert all(row["status"] == "ok" for row in batch)
    assert np.std([row["unimol_score"] for row in batch]) > 0


@pytest.mark.skipif(not __import__("torch").backends.mps.is_available(), reason="MPS unavailable")
@pytest.mark.integration_assets
def test_unimol_mps_embeddings_match_cpu(monkeypatch):
    weight_dir = (
        SCORING_DIR
        / "models"
        / "unimol"
        / "models--dptech--Uni-Mol-Models"
        / "snapshots"
        / "9f19c45c718192888a1c8a1c905f69f0755ea502"
    )
    monkeypatch.setenv("UNIMOL_WEIGHT_DIR", str(weight_dir))
    sys.path.insert(0, str(SCORING_DIR))
    from scripts.unimol_embedding import UniMolEmbedding

    mps = UniMolEmbedding(
        str(weight_dir / "mol_pre_all_h_220816.pt"),
        str(weight_dir / "mol.dict.txt"),
        device="mps",
        batch_size=2,
    )
    assert mps.device == "mps"
    mps_values = mps.extract(["CCO", "c1ccccc1"], batch_size=2)

    cpu = UniMolEmbedding(
        str(weight_dir / "mol_pre_all_h_220816.pt"),
        str(weight_dir / "mol.dict.txt"),
        device="cpu",
        batch_size=2,
    )
    cpu_values = cpu.extract(["CCO", "c1ccccc1"], batch_size=2)

    assert np.allclose(mps_values, cpu_values, atol=1e-3)


@pytest.mark.skipif(not __import__("torch").backends.mps.is_available(), reason="MPS unavailable")
@pytest.mark.integration_assets
def test_unimol_mps_representation_is_batch_context_stable(monkeypatch):
    weight_dir = (
        SCORING_DIR
        / "models"
        / "unimol"
        / "models--dptech--Uni-Mol-Models"
        / "snapshots"
        / "9f19c45c718192888a1c8a1c905f69f0755ea502"
    )
    monkeypatch.setenv("UNIMOL_WEIGHT_DIR", str(weight_dir))
    sys.path.insert(0, str(SCORING_DIR))
    from scripts.unimol_embedding import UniMolEmbedding

    extractor = UniMolEmbedding(
        str(weight_dir / "mol_pre_all_h_220816.pt"),
        str(weight_dir / "mol.dict.txt"),
        device="mps",
        batch_size=2,
    )
    single = extractor.extract(["CCO"], batch_size=1)[0]
    mixed = extractor.extract(["CCO", "Nc1ccccc1C(=O)NCCCCCCCCCCCCCCCCCC"], batch_size=2)[0]

    assert np.allclose(single, mixed, atol=1e-4)


def test_unimol_mps_padding_uses_stable_length_buckets():
    sys.path.insert(0, str(SCORING_DIR))
    from scripts.unimol_embedding import _stable_padding_length

    assert _stable_padding_length(12, max_tokens=258) == 32
    assert _stable_padding_length(32, max_tokens=258) == 32
    assert _stable_padding_length(33, max_tokens=258) == 64
    assert _stable_padding_length(100, max_tokens=258) == 128
    assert _stable_padding_length(129, max_tokens=258) == 258
    assert _stable_padding_length(258, max_tokens=258) == 258


def test_unimol_oversized_atom_crop_is_deterministic():
    sys.path.insert(0, str(SCORING_DIR))
    from scripts.unimol_embedding import _deterministic_crop_indices

    first = _deterministic_crop_indices(400, 256)
    second = _deterministic_crop_indices(400, 256)
    assert np.array_equal(first, second)
    assert len(first) == 256
    assert first[0] == 0
    assert first[-1] == 399


@pytest.mark.integration_assets
def test_layer_cache_matches_production_scorers(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import layer_cache

    cache = layer_cache.build(
        ["CCO", "c1ccccc1"],
        tmp_path,
        batch_size=2,
        scoring_dir=SCORING_DIR,
    )

    assert set(cache.columns) >= {
        "smiles", "l1", "l3", "l4", "mol_features", "molfeat_status",
        "l1_status", "l3_status", "l4_status",
    }
    assert (cache[["l1_status", "l3_status", "l4_status"]] == "ok").all().all()
    assert cache["molfeat_status"].eq("ok").all()
    assert all(len(features) == 520 for features in cache["mol_features"])
    assert cache["l4"].nunique() > 1


def test_layer_cache_metadata_rejects_changed_asset_namespace(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import layer_cache

    layer_cache.write_cache_metadata(
        tmp_path,
        namespace="namespace-a",
        contract={"model_sha256": "aaa"},
        canonical_smiles=["CCO", "c1ccccc1"],
        shard_size=2,
        status="in_progress",
    )
    assert layer_cache.require_cache_metadata(
        tmp_path,
        namespace="namespace-a",
        canonical_smiles=["CCO", "c1ccccc1"],
        shard_size=2,
    )["ok"] is True

    with pytest.raises(RuntimeError, match="namespace mismatch"):
        layer_cache.require_cache_metadata(
            tmp_path,
            namespace="namespace-b",
            canonical_smiles=["CCO", "c1ccccc1"],
            shard_size=2,
        )


def test_layer3_batch_matches_scalar_fields():
    _, scoring = _import_scoring_modules()
    scorer = scoring.Layer3Scorer()
    assert hasattr(scorer, "score_many")
    smiles = [
        "CCO",
        "CC(=O)Oc1ccccc1C(=O)O",
        "O=C1NC(=O)C2=C1CCN(C2=O)c1ccccc1",
        "C[N+]1=CC=C(C=C1)c2cc[n+](C)cc2",
    ]

    scalar = [scorer.score(smi) for smi in smiles]
    batch = scorer.score_many(smiles)

    assert batch == scalar


@pytest.mark.integration_assets
def test_batch_l2_and_final_match_production_formula():
    from scientific_validation.four_level_cli_1kx10k import batch_cli

    _, scoring = _import_scoring_modules()
    from l2_bindingdb import BindingDBFeature, Layer2BindingDB

    smiles = [
        "C", "CC", "CCC", "CCCC", "CCO", "CCN", "COC", "CNC", "CCCl", "CCBr",
        "c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1", "CC(=O)O",
        "CC(=O)N", "C1CCCCC1", "c1ccncc1", "CC(C)O", "CC(C)N",
    ]
    target_text = "CHEMBL2051 Neuraminidase SINGLE PROTEIN Influenza A virus"
    featurizer = BindingDBFeature()
    cache = pd.DataFrame(
        {
            "smiles": smiles,
            "mol_features": [featurizer.mol_features(smi).tolist() for smi in smiles],
            "l1": np.linspace(0.35, 0.85, len(smiles)),
            "l3": np.linspace(0.45, 0.9, len(smiles)),
            "l4": np.linspace(0.2, 0.75, len(smiles)),
            "l1_status": "ok",
            "l3_status": "ok",
            "l4_status": "ok",
            "l3_details": [{} for _ in smiles],
        }
    )

    actual = batch_cli.score_target_frame(
        cache,
        target_id="CHEMBL2051",
        target_text=target_text,
        scoring_dir=SCORING_DIR,
    )
    production_l2 = Layer2BindingDB(prefer="mlp")
    expected_l2 = np.array(
        [production_l2.score(smi, target_text=target_text)["docking_normalized"] for smi in smiles]
    )
    expected_final = scoring.Layer4Aggregator.combine(
        cache["l1"], expected_l2, cache["l3"], cache["l4"]
    )

    assert np.allclose(actual["l2"], expected_l2, atol=1e-4)
    assert np.allclose(actual["final_score"], expected_final, atol=1e-4)
    assert actual["layer2_status"].eq("ok").all()
    assert set(actual["gate_status"]) <= {"PASS", "FAIL"}


@pytest.mark.integration_assets
def test_l2_compat_model_matches_frozen_original():
    import joblib
    import warnings

    sys.path.insert(0, str(SCORING_DIR))
    from l2_bindingdb import BindingDBFeature

    model_dir = SCORING_DIR / "models" / "bindingdb_l2"
    compat_path = model_dir / "l2_model_sklearn_1_7_2.joblib"
    assert compat_path.is_file()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        original = joblib.load(model_dir / "l2_model.joblib")["mlp"]
    compat = joblib.load(compat_path)["mlp"]
    featurizer = BindingDBFeature()
    features = np.stack(
        [
            featurizer.features("CCO", "CHEMBL2051 Neuraminidase"),
            featurizer.features("c1ccccc1", "CHEMBL2051 Neuraminidase"),
        ]
    )

    assert np.array_equal(original.predict_proba(features), compat.predict_proba(features))


@pytest.mark.integration_assets
def test_default_l2_uses_compatible_asset_without_version_warning():
    import warnings
    from sklearn.exceptions import InconsistentVersionWarning

    sys.path.insert(0, str(SCORING_DIR))
    from l2_bindingdb import Layer2BindingDB

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        model = Layer2BindingDB(prefer="mlp")
        model._ensure_model()

    assert model.model_path.name == "l2_model_sklearn_1_7_2.joblib"
    assert not any(isinstance(item.message, InconsistentVersionWarning) for item in caught)


def test_parity_comparison_covers_every_four_level_field():
    from scientific_validation.four_level_cli_1kx10k import batch_cli

    common = {
        "canonical_smiles": ["CCO", "c1ccccc1"],
        "l1": [0.8, 0.7],
        "l2": [0.6, 0.5],
        "l3": [0.9, 0.75],
        "l4": [0.4, 0.3],
        "final_score": [0.68, 0.565],
        "gate_status": ["PASS", "PASS"],
        "l1_status": ["ok", "ok"],
        "layer2_status": ["ok", "ok"],
        "l3_status": ["ok", "ok"],
        "l4_status": ["ok", "ok"],
        "l1_model_asset_id": ["rdkit-runtime"] * 2,
        "layer2_model_asset_id": ["l2_model_sklearn_1_7_2.joblib"] * 2,
        "l3_model_asset_id": ["tox21.pkl+bbbp.pkl+clintox.pkl+sider.pkl"] * 2,
        "l4_model_asset_id": ["mol_pre_all_h_220816.pt+ref_embeddings.npz"] * 2,
    }
    batch = pd.DataFrame(common)
    scalar = pd.DataFrame(common)

    evidence = batch_cli.compare_parity_outputs(batch, scalar, tolerance=1e-4)

    assert evidence["passed"] is True
    assert set(evidence["numeric_fields"]) == {"l1", "l2", "l3", "l4", "final_score"}
    assert set(evidence["categorical_fields"]) >= {
        "gate_status", "l1_status", "layer2_status", "l3_status", "l4_status",
        "l1_model_asset_id", "layer2_model_asset_id", "l3_model_asset_id", "l4_model_asset_id",
    }
    scalar.loc[1, "l4"] += 0.01
    assert batch_cli.compare_parity_outputs(batch, scalar, tolerance=1e-4)["passed"] is False


def test_batch_cli_exposes_checkpointed_run_stages():
    result = subprocess.run(
        [sys.executable, "-m", "scientific_validation.four_level_cli_1kx10k.batch_cli", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert result.returncode == 0, result.stderr
    for command in ("prepare", "cache-layers", "score", "dock", "report"):
        assert command in result.stdout


def test_report_module_exposes_cli_help():
    result = subprocess.run(
        [sys.executable, "-m", "scientific_validation.four_level_cli_1kx10k.report", "--help"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert result.returncode == 0, result.stderr
    assert "--run-dir" in result.stdout


def test_docking_selection_uses_only_stable_l2_top_n():
    from scientific_validation.four_level_cli_1kx10k import batch_cli

    scores = pd.DataFrame(
        {
            "candidate_index": [0, 1, 2, 3],
            "canonical_smiles": ["C", "CC", "CCC", "CCCC"],
            "l2": [0.2, 0.9, 0.5, 0.5],
        }
    )

    selected = batch_cli.select_docking_rows(scores, top_n=3)

    assert selected["candidate_index"].tolist() == [1, 2, 3]
    assert selected["selection_rank"].tolist() == [1, 2, 3]


def test_docking_fusion_changes_only_successful_rows_and_recomputes_metrics():
    from scientific_validation.four_level_cli_1kx10k import batch_cli

    scores = _frame = pd.DataFrame(
        {
            "candidate_index": [0, 1, 2, 3, 4, 5],
            "canonical_smiles": ["C", "CC", "CCC", "CCCC", "CCO", "CCN"],
            "label_role": [
                "heldout_positive", "heldout_positive", "heldout_negative",
                "heldout_negative", "unlabeled_background", "unlabeled_background",
            ],
            "l2": [0.8, 0.7, 0.6, 0.5, 0.4, 0.3],
            "final_score": [0.75, 0.65, 0.55, 0.45, 0.35, 0.25],
        }
    )
    docking = pd.DataFrame(
        {
            "candidate_index": [0, 1, 2],
            "status": ["ok", "ok", "prep_failed"],
            "ligand_efficiency": [-0.20, -0.50, None],
        }
    )

    fused, comparison = batch_cli.fuse_docking_results(scores, docking, weight=0.3)

    assert fused.loc[2, "fused_score"] == scores.loc[2, "l2"]
    assert fused.loc[3, "fused_score"] == scores.loc[3, "l2"]
    assert fused.loc[0, "fused_score"] != scores.loc[0, "l2"]
    assert fused.loc[1, "fused_score"] != scores.loc[1, "l2"]
    assert comparison["before"].get("auc_fused_labeled") is None
    assert comparison["after"]["auc_fused_labeled"] is not None
