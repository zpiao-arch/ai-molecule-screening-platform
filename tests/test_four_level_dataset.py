from pathlib import Path

import pandas as pd
import pytest

from scientific_validation.four_level_cli_1kx10k.contracts import BenchmarkConfig
from scientific_validation.four_level_cli_1kx10k import dataset


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.integration_assets
def test_reconstructed_split_matches_production_model_metadata():
    config = BenchmarkConfig(root=ROOT)

    split = dataset.reconstruct_combined_split(config)

    assert len(split.fit_rows) == 62999
    assert len(split.heldout_rows) == 15750
    assert split.exact_key_overlap_count == 0
    assert split.eligible_fit_pair_overlap_count == 0
    assert set(split.eligible_heldout_rows["label"].unique()) <= {0, 1}


def test_target_pool_is_exact_deterministic_and_excludes_fit_pairs():
    fit_rows = pd.DataFrame(
        [
            {"target_id": "CHEMBL1", "canonical_smiles": "C", "label": 1, "target_text": "t1"},
            {"target_id": "CHEMBL1", "canonical_smiles": "CO", "label": 0, "target_text": "t1"},
        ]
    )
    heldout_rows = pd.DataFrame(
        [
            {"target_id": "CHEMBL1", "canonical_smiles": "CC", "label": 1, "target_text": "t1"},
            {"target_id": "CHEMBL1", "canonical_smiles": "CCC", "label": 0, "target_text": "t1"},
        ]
    )
    library = pd.DataFrame(
        {"canonical_smiles": ["C", "CO", "CC", "CCC", "CCCC", "CCO", "CCN", "CN", "CNC", "COC", "CCCl", "CCBr"]}
    )

    first = dataset.build_target_pool(
        target_id="CHEMBL1",
        target_index=0,
        fit_rows=fit_rows,
        heldout_rows=heldout_rows,
        library=library,
        pool_size=8,
        seed=42,
    )
    second = dataset.build_target_pool(
        target_id="CHEMBL1",
        target_index=0,
        fit_rows=fit_rows,
        heldout_rows=heldout_rows,
        library=library,
        pool_size=8,
        seed=42,
    )

    assert len(first) == 8
    assert first.equals(second)
    assert first["canonical_smiles"].is_unique
    assert not {"C", "CO"} & set(first["canonical_smiles"])
    roles = dict(zip(first["canonical_smiles"], first["label_role"]))
    assert roles["CC"] == "heldout_positive"
    assert roles["CCC"] == "heldout_negative"
    assert set(first["label_role"]) <= {
        "heldout_positive",
        "heldout_negative",
        "unlabeled_background",
    }


def test_select_targets_prioritizes_strictly_evaluable_then_real_fallbacks():
    rows = []
    for target_id, n_pos, n_neg in (("CHEMBL1", 5, 5), ("CHEMBL2", 2, 5), ("CHEMBL3", 1, 0)):
        rows.extend(
            {"target_id": target_id, "canonical_smiles": f"{target_id}P{i}", "label": 1, "target_text": target_id}
            for i in range(n_pos)
        )
        rows.extend(
            {"target_id": target_id, "canonical_smiles": f"{target_id}N{i}", "label": 0, "target_text": target_id}
            for i in range(n_neg)
        )
    heldout = pd.DataFrame(rows)

    targets = dataset.select_targets(
        heldout,
        n_targets=5,
        fallback_target_ids=["CHEMBL4", "CHEMBL5", "CHEMBL6"],
    )

    assert targets[0]["target_id"] == "CHEMBL1"
    assert targets[0]["evaluation_tier"] == "strict_labeled"
    assert [row["target_id"] for row in targets] == ["CHEMBL1", "CHEMBL2", "CHEMBL3", "CHEMBL4", "CHEMBL5"]


def test_materialized_prepare_run_writes_exact_independent_pools(tmp_path):
    from scientific_validation.four_level_cli_1kx10k import batch_cli

    fit_rows = pd.DataFrame(
        [
            {"target_id": "CHEMBL1", "canonical_smiles": "C", "label": 1, "target_text": "t1"},
            {"target_id": "CHEMBL2", "canonical_smiles": "CO", "label": 0, "target_text": "t2"},
        ]
    )
    heldout_rows = pd.DataFrame(
        [
            {"target_id": "CHEMBL1", "canonical_smiles": "CC", "label": 1, "target_text": "t1"},
            {"target_id": "CHEMBL1", "canonical_smiles": "CCC", "label": 0, "target_text": "t1"},
            {"target_id": "CHEMBL2", "canonical_smiles": "CCO", "label": 1, "target_text": "t2"},
        ]
    )
    split = dataset.SplitResult(fit_rows, heldout_rows, heldout_rows, 0, 0)
    targets = [
        {"target_id": "CHEMBL1", "target_text": "t1", "evaluation_tier": "strict_labeled"},
        {"target_id": "CHEMBL2", "target_text": "t2", "evaluation_tier": "limited_positive"},
    ]
    library = pd.DataFrame(
        {"canonical_smiles": ["C", "CO", "CC", "CCC", "CCO", "CCN", "CN", "CNC", "COC", "CCCl", "CCBr"]}
    )

    summary = batch_cli.materialize_prepared_run(
        run_dir=tmp_path,
        split=split,
        targets=targets,
        library=library,
        pool_size=6,
        seed=42,
    )

    assert summary["n_targets"] == 2
    assert summary["n_pairs"] == 12
    manifest = pd.read_parquet(tmp_path / "target_manifest.parquet")
    assert len(manifest) == 2
    pools = sorted((tmp_path / "pool_manifest").glob("*.parquet"))
    assert len(pools) == 2
    assert all(len(pd.read_parquet(path)) == 6 for path in pools)
