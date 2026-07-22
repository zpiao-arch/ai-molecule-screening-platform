from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
import sys
from pathlib import Path
from typing import Iterable, Sequence

import pandas as pd
from rdkit import Chem

from . import provenance


def _canonicalize(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)


def _cache_key(smiles: str, namespace: str) -> str:
    return hashlib.sha256(f"{smiles}\0{namespace}".encode("utf-8")).hexdigest()


def cache_contract(scoring_dir: str | Path) -> dict[str, object]:
    scoring_dir = Path(scoring_dir).resolve()
    snapshot = (
        scoring_dir
        / "models"
        / "unimol"
        / "models--dptech--Uni-Mol-Models"
        / "snapshots"
        / "9f19c45c718192888a1c8a1c905f69f0755ea502"
    )
    paths = [
        scoring_dir / "scoring.py",
        scoring_dir / "l2_bindingdb.py",
        scoring_dir / "scripts" / "unimol_embedding.py",
        scoring_dir / "scripts" / "unimol_scorer.py",
        *(scoring_dir / "models" / "admet" / f"{name}.pkl" for name in ("tox21", "bbbp", "clintox", "sider")),
        snapshot / "mol_pre_all_h_220816.pt",
        snapshot / "mol.dict.txt",
        scoring_dir / "models" / "ref_embeddings.npz",
        scoring_dir / "models" / "ref_smiles.pkl",
    ]
    assets = {}
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(f"cache asset missing: {path}")
        assets[str(path.relative_to(scoring_dir))] = provenance.sha256_file(path)
    return {
        "assets": assets,
        "packages": {
            "rdkit": importlib.metadata.version("rdkit"),
            "scikit-learn": importlib.metadata.version("scikit-learn"),
            "torch": importlib.metadata.version("torch"),
            "unimol-tools": importlib.metadata.version("unimol-tools"),
        },
    }


def cache_namespace(scoring_dir: str | Path) -> tuple[str, dict[str, object]]:
    contract = cache_contract(scoring_dir)
    return provenance.payload_sha256(contract), contract


def _cache_metadata_path(output: str | Path) -> Path:
    return Path(output) / "cache_metadata.json"


def write_cache_metadata(
    output: str | Path,
    *,
    namespace: str,
    contract: dict[str, object],
    canonical_smiles: Sequence[str],
    shard_size: int,
    status: str,
    n_shards: int | None = None,
) -> Path:
    if status not in {"in_progress", "complete"}:
        raise ValueError("cache metadata status must be in_progress or complete")
    destination = _cache_metadata_path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "status": status,
        "namespace": namespace,
        "contract": contract,
        "canonical_smiles_sha256": provenance.payload_sha256(list(canonical_smiles)),
        "n_molecules": len(canonical_smiles),
        "shard_size": int(shard_size),
        "n_shards": n_shards,
    }
    temporary = destination.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    os.replace(temporary, destination)
    return destination


def validate_cache_metadata(
    output: str | Path,
    *,
    namespace: str,
    canonical_smiles: Sequence[str],
    shard_size: int,
) -> dict[str, object]:
    path = _cache_metadata_path(output)
    if not path.is_file():
        return {"ok": False, "reason": "metadata missing", "path": str(path)}
    payload = json.loads(path.read_text(encoding="utf-8"))
    reasons = []
    if payload.get("namespace") != namespace:
        reasons.append("namespace mismatch")
    if payload.get("canonical_smiles_sha256") != provenance.payload_sha256(list(canonical_smiles)):
        reasons.append("canonical SMILES fingerprint mismatch")
    if int(payload.get("shard_size", -1)) != int(shard_size):
        reasons.append("shard size mismatch")
    return {"ok": not reasons, "reason": "; ".join(reasons) if reasons else "ok", "metadata": payload}


def require_cache_metadata(
    output: str | Path,
    *,
    namespace: str,
    canonical_smiles: Sequence[str],
    shard_size: int,
) -> dict[str, object]:
    result = validate_cache_metadata(
        output,
        namespace=namespace,
        canonical_smiles=canonical_smiles,
        shard_size=shard_size,
    )
    if not result["ok"]:
        raise RuntimeError(f"cache metadata {result['reason']}")
    return result


def _atomic_parquet(frame: pd.DataFrame, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    frame.to_parquet(temporary, index=False, engine="pyarrow", compression="zstd")
    os.replace(temporary, destination)


def _import_production(scoring_dir: Path):
    path = str(scoring_dir.resolve())
    if path not in sys.path:
        sys.path.insert(0, path)
    import scoring
    from l2_bindingdb import BindingDBFeature
    from scripts.unimol_scorer import UniMolScorer

    return scoring, BindingDBFeature, UniMolScorer


def _score_shard(
    smiles: Sequence[str],
    *,
    start_index: int,
    l1_scorer,
    l3_scorer,
    mol_featurizer,
    unimol_scorer,
    batch_size: int,
    namespace: str,
) -> pd.DataFrame:
    l3_results = l3_scorer.score_many(list(smiles))
    base_rows = []
    for offset, smi in enumerate(smiles):
        l1 = l1_scorer.score(smi)
        l3 = l3_results[offset]
        mol_features = mol_featurizer.mol_features(smi)
        base_rows.append(
            {
                "input_index": start_index + offset,
                "smiles": smi,
                "cache_key": _cache_key(smi, namespace),
                "l1": l1.get("layer1_score"),
                "l1_status": "ok" if l1.get("valid") else "error:invalid_smiles",
                "l1_backend": "RDKit",
                "l1_model_asset_id": "rdkit-runtime",
                "l1_details": l1,
                "l3": l3.get("admet_score") if l3.get("toxicity_flags") != "invalid" else None,
                "l3_status": "ok" if l3.get("toxicity_flags") != "invalid" else "error:invalid_smiles",
                "l3_backend": "ADMET-RF+RDKit",
                "l3_model_asset_id": "tox21.pkl+bbbp.pkl+clintox.pkl+sider.pkl",
                "l3_details": l3,
                "mol_features": mol_features.tolist() if mol_features is not None else None,
                "molfeat_status": "ok" if mol_features is not None else "error:invalid_smiles",
                "molfeat_backend": "BindingDBFeature-RDKit",
            }
        )

    l4_rows = unimol_scorer.score_many(list(smiles), batch_size=batch_size)
    if len(l4_rows) != len(base_rows):
        raise RuntimeError("Uni-Mol batch result count does not match input count")
    for row, l4 in zip(base_rows, l4_rows):
        status = str(l4.get("status", "error:missing_status"))
        row.update(
            {
                "l4": l4.get("unimol_score") if status == "ok" else None,
                "l4_status": status,
                "l4_backend": "UniMolRepr-mol_pre_all_h_220816",
                "l4_model_asset_id": "mol_pre_all_h_220816.pt+ref_embeddings.npz",
                "l4_pos_similarity": l4.get("pos_similarity"),
                "l4_top_refs": json.dumps(l4.get("top_refs", []), ensure_ascii=True),
            }
        )
    return pd.DataFrame(base_rows)


def build(
    smiles: Iterable[str],
    out_dir: str | Path,
    *,
    batch_size: int = 32,
    shard_size: int = 512,
    resume: bool = True,
    strict: bool = True,
    multi_process: bool = False,
    unimol_device: str = "mps",
    scoring_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Build or resume the target-independent four-level cache."""
    if batch_size < 1 or shard_size < 1:
        raise ValueError("batch_size and shard_size must be positive")
    scoring_dir = Path(scoring_dir or os.environ.get(
        "FOUR_LEVEL_SCORING_DIR",
        Path(__file__).resolve().parents[2] / "scoring",
    ))
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    shard_dir = output / "layer_cache_shards"
    shard_dir.mkdir(parents=True, exist_ok=True)

    canonical = list(dict.fromkeys(_canonicalize(smi) for smi in smiles))
    namespace, contract = cache_namespace(scoring_dir)
    existing_shards = sorted(shard_dir.glob("part-*.parquet"))
    if existing_shards:
        if not resume:
            raise RuntimeError("layer cache shards already exist; use --resume with matching metadata")
        require_cache_metadata(
            output,
            namespace=namespace,
            canonical_smiles=canonical,
            shard_size=shard_size,
        )
    else:
        write_cache_metadata(
            output,
            namespace=namespace,
            contract=contract,
            canonical_smiles=canonical,
            shard_size=shard_size,
            status="in_progress",
        )
    scoring, BindingDBFeature, UniMolScorer = _import_production(scoring_dir)
    l1_scorer = scoring.Layer1Scorer()
    l3_scorer = scoring.Layer3Scorer()
    mol_featurizer = BindingDBFeature()
    unimol_scorer = UniMolScorer(device=unimol_device, multi_process=multi_process)

    shards = []
    for start in range(0, len(canonical), shard_size):
        expected = canonical[start : start + shard_size]
        shard_path = shard_dir / f"part-{start // shard_size:06d}.parquet"
        if resume and shard_path.is_file():
            frame = pd.read_parquet(shard_path)
            if frame["smiles"].tolist() != expected:
                raise RuntimeError(f"cache shard input mismatch: {shard_path}")
        else:
            frame = _score_shard(
                expected,
                start_index=start,
                l1_scorer=l1_scorer,
                l3_scorer=l3_scorer,
                mol_featurizer=mol_featurizer,
                unimol_scorer=unimol_scorer,
                batch_size=batch_size,
                namespace=namespace,
            )
            if strict and frame["l4_status"].eq("ok").sum() == 0 and len(frame):
                raise RuntimeError(f"Uni-Mol failed for every molecule in shard {start // shard_size}")
            _atomic_parquet(frame, shard_path)
        shards.append(frame)

    if not shards:
        return pd.DataFrame()
    result = pd.concat(shards, ignore_index=True).sort_values("input_index").reset_index(drop=True)
    if len(result) != len(canonical):
        raise RuntimeError("layer cache row count mismatch")
    write_cache_metadata(
        output,
        namespace=namespace,
        contract=contract,
        canonical_smiles=canonical,
        shard_size=shard_size,
        status="complete",
        n_shards=len(shards),
    )
    return result


def seal_existing_cache(
    out_dir: str | Path,
    *,
    scoring_dir: str | Path,
    canonical_smiles: Sequence[str],
    shard_size: int,
) -> dict[str, object]:
    output = Path(out_dir)
    namespace, contract = cache_namespace(scoring_dir)
    expected = list(canonical_smiles)
    paths = sorted((output / "layer_cache_shards").glob("part-*.parquet"))
    rows = 0
    observed = []
    for path in paths:
        frame = pd.read_parquet(path)
        observed.extend(frame["smiles"].astype(str).tolist())
        frame["cache_key"] = [_cache_key(smi, namespace) for smi in frame["smiles"].astype(str)]
        _atomic_parquet(frame, path)
        rows += len(frame)
    if observed != expected:
        raise RuntimeError("existing cache SMILES do not match frozen library order")
    write_cache_metadata(
        output,
        namespace=namespace,
        contract=contract,
        canonical_smiles=expected,
        shard_size=shard_size,
        status="complete",
        n_shards=len(paths),
    )
    return {"namespace": namespace, "n_molecules": rows, "n_shards": len(paths)}
