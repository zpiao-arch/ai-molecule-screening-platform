from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

from . import batch_cli, layer_cache, provenance
from .contracts import BenchmarkConfig


def _atomic_parquet(frame: pd.DataFrame, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    frame.to_parquet(temporary, index=False, engine="pyarrow", compression="zstd")
    temporary.replace(destination)


def repair(run_dir: str | Path) -> dict[str, object]:
    root = Path(run_dir).resolve()
    config = BenchmarkConfig(root=Path(__file__).resolve().parents[2])
    library = pd.read_parquet(root / "library_manifest.parquet")
    oversized_indices: list[int] = []
    oversized_smiles: list[str] = []
    for index, smiles in enumerate(library["canonical_smiles"].astype(str)):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None and Chem.AddHs(molecule).GetNumAtoms() > 256:
            oversized_indices.append(index)
            oversized_smiles.append(smiles)

    scoring_dir = config.scoring_dir.resolve()
    if str(scoring_dir) not in sys.path:
        sys.path.insert(0, str(scoring_dir))
    from scripts.unimol_scorer import UniMolScorer

    scorer = UniMolScorer(device="mps", multi_process=True)
    l4_rows = scorer.score_many(oversized_smiles, batch_size=128)
    if any(str(row.get("status")) != "ok" for row in l4_rows):
        raise RuntimeError("deterministic oversized-molecule L4 repair contains failures")
    by_index = dict(zip(oversized_indices, l4_rows))

    namespace, _ = layer_cache.cache_namespace(scoring_dir)
    shard_dir = root / "layer_cache" / "layer_cache_shards"
    repaired = 0
    for shard_path in sorted(shard_dir.glob("part-*.parquet")):
        frame = pd.read_parquet(shard_path)
        for row_index, record in frame.iterrows():
            input_index = int(record["input_index"])
            if input_index not in by_index:
                continue
            l4 = by_index[input_index]
            frame.at[row_index, "l4"] = l4["unimol_score"]
            frame.at[row_index, "l4_status"] = l4["status"]
            frame.at[row_index, "l4_pos_similarity"] = l4["pos_similarity"]
            frame.at[row_index, "l4_top_refs"] = json.dumps(l4.get("top_refs", []), ensure_ascii=True)
            repaired += 1
        frame["cache_key"] = [layer_cache._cache_key(str(smiles), namespace) for smiles in frame["smiles"]]
        _atomic_parquet(frame, shard_path)

    if repaired != len(oversized_indices):
        raise RuntimeError(f"repaired rows {repaired} != oversized molecules {len(oversized_indices)}")

    layer_cache.seal_existing_cache(
        root / "layer_cache",
        scoring_dir=scoring_dir,
        canonical_smiles=library["canonical_smiles"].astype(str).tolist(),
        shard_size=4096,
    )
    cache = pd.concat(
        (pd.read_parquet(path) for path in sorted(shard_dir.glob("part-*.parquet"))),
        ignore_index=True,
    )
    summary = {
        "n_molecules": int(len(cache)),
        "l1_failures": int((cache["l1_status"] != "ok").sum()),
        "l3_failures": int((cache["l3_status"] != "ok").sum()),
        "l4_failures": int((cache["l4_status"] != "ok").sum()),
        "l4_std": float(pd.to_numeric(cache["l4"], errors="coerce").std()),
        "deterministic_oversized_repaired": int(repaired),
    }
    summary_path = root / "layer_cache_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")

    inputs = {
        "library_manifest_sha256": batch_cli._sha256(root / "library_manifest.parquet"),
        "n_molecules": int(len(library)),
        "batch_size": 128,
        "shard_size": 4096,
        "multi_process": True,
        "unimol_device": "mps",
        "cache_namespace": namespace,
    }
    provenance.write_stage_checkpoint(
        root,
        "layer_cache",
        inputs=inputs,
        status="in_progress",
        metadata={"deterministic_oversized_repaired": int(repaired)},
    )
    result = batch_cli.cache_run_layers(
        config,
        root,
        batch_size=128,
        shard_size=4096,
        resume=True,
        strict=True,
        multi_process=True,
        unimol_device="mps",
    )
    result["deterministic_oversized_repaired"] = int(repaired)
    return result


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args()
    print(json.dumps(repair(args.run_dir), indent=2, sort_keys=True))
