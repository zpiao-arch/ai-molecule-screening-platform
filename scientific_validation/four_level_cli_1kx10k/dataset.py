from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from sklearn.model_selection import train_test_split

from .contracts import BenchmarkConfig


RDLogger.DisableLog("rdApp.*")


@dataclass(frozen=True)
class SplitResult:
    fit_rows: pd.DataFrame
    heldout_rows: pd.DataFrame
    eligible_heldout_rows: pd.DataFrame
    exact_key_overlap_count: int
    eligible_fit_pair_overlap_count: int


def _canonicalize(smiles: str) -> str | None:
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol, canonical=True) if mol is not None else None
    except Exception:
        return None


def _load_uniprot_to_chembl(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                mapping.setdefault(fields[0].strip(), fields[1].strip())
    return mapping


def _parse_label(value: object) -> int | None:
    try:
        label = int(float(str(value)))
    except (TypeError, ValueError):
        return None
    return label if label in (0, 1) else None


def _target_id_from_bindingdb(row: Mapping[str, str], uniprot_map: Mapping[str, str]) -> str | None:
    text = row.get("target_text", "") or ""
    match = re.search(r"CHEMBL\d+", text)
    if match:
        return match.group(0)
    uniprots = row.get("UniProt (SwissProt) Primary ID of Target Chain 1", "") or ""
    for uniprot in re.split(r"[,;\s]+", uniprots.strip()):
        if uniprot in uniprot_map:
            return uniprot_map[uniprot]
    return None


def _iter_source_rows(config: BenchmarkConfig) -> Iterator[dict[str, object]]:
    uniprot_map = _load_uniprot_to_chembl(config.chembl_uniprot_mapping)
    sources = (
        ("bindingdb", config.bindingdb_examples),
        ("chembl", config.chembl_examples),
    )
    for source, path in sources:
        with path.open(encoding="utf-8", newline="") as handle:
            for source_row, row in enumerate(csv.DictReader(handle)):
                label = _parse_label(row.get("label"))
                smiles = (row.get("canonical_smiles") or "").strip()
                target_text = (row.get("target_text") or "").strip()
                if label is None or not smiles or not target_text:
                    continue
                canonical = _canonicalize(smiles)
                if canonical is None:
                    continue
                if source == "chembl":
                    target_id = (row.get("target_chembl_id") or "").strip() or None
                else:
                    target_id = _target_id_from_bindingdb(row, uniprot_map)
                yield {
                    "source": source,
                    "source_row": source_row,
                    "source_smiles": smiles,
                    "canonical_smiles": canonical,
                    "target_text": target_text,
                    "target_id": target_id,
                    "label": label,
                }


def _pair_keys(frame: pd.DataFrame, smiles_column: str) -> set[tuple[str, str]]:
    if frame.empty:
        return set()
    return set(zip(frame["target_id"].astype(str), frame[smiles_column].astype(str)))


def reconstruct_combined_split(config: BenchmarkConfig) -> SplitResult:
    """Mirror the fitted model's row order and split, then build a clean evaluation view."""
    seen: set[tuple[str, str]] = set()
    rows: list[dict[str, object]] = []
    for row in _iter_source_rows(config):
        exact_key = (str(row["source_smiles"]), str(row["target_text"]))
        if exact_key in seen:
            continue
        seen.add(exact_key)
        rows.append(row)

    all_rows = pd.DataFrame(rows)
    if all_rows.empty:
        raise RuntimeError("No valid BindingDB/ChEMBL rows were reconstructed")

    indices = np.arange(len(all_rows))
    fit_indices, heldout_indices = train_test_split(
        indices,
        test_size=0.2,
        random_state=0,
        stratify=all_rows["label"].to_numpy(),
    )
    fit_rows = all_rows.iloc[fit_indices].reset_index(drop=True)
    heldout_rows = all_rows.iloc[heldout_indices].reset_index(drop=True)

    fit_exact = set(zip(fit_rows["source_smiles"], fit_rows["target_text"]))
    heldout_exact = set(zip(heldout_rows["source_smiles"], heldout_rows["target_text"]))
    exact_overlap = len(fit_exact & heldout_exact)

    fit_pairs = _pair_keys(fit_rows.dropna(subset=["target_id"]), "canonical_smiles")
    eligible = heldout_rows.dropna(subset=["target_id"]).copy()
    eligible["fit_pair"] = list(zip(eligible["target_id"], eligible["canonical_smiles"]))
    eligible = eligible.loc[~eligible["fit_pair"].isin(fit_pairs)].drop(columns="fit_pair")

    label_counts = eligible.groupby(["target_id", "canonical_smiles"])["label"].nunique()
    conflict_keys = set(label_counts[label_counts > 1].index)
    if conflict_keys:
        keys = list(zip(eligible["target_id"], eligible["canonical_smiles"]))
        eligible = eligible.loc[[key not in conflict_keys for key in keys]]
    eligible = eligible.drop_duplicates(["target_id", "canonical_smiles"], keep="first")
    eligible = eligible.reset_index(drop=True)

    eligible_overlap = len(_pair_keys(eligible, "canonical_smiles") & fit_pairs)
    return SplitResult(
        fit_rows=fit_rows,
        heldout_rows=heldout_rows,
        eligible_heldout_rows=eligible,
        exact_key_overlap_count=exact_overlap,
        eligible_fit_pair_overlap_count=eligible_overlap,
    )


def build_target_pool(
    *,
    target_id: str,
    target_index: int,
    fit_rows: pd.DataFrame,
    heldout_rows: pd.DataFrame,
    library: pd.DataFrame,
    pool_size: int,
    seed: int,
) -> pd.DataFrame:
    """Build one deterministic pool without treating background as confirmed negative."""
    fit = fit_rows.loc[fit_rows["target_id"] == target_id]
    heldout = heldout_rows.loc[heldout_rows["target_id"] == target_id].copy()
    fit_smiles = set(fit["canonical_smiles"].dropna().astype(str))
    heldout = heldout.loc[~heldout["canonical_smiles"].isin(fit_smiles)]

    counts = heldout.groupby("canonical_smiles")["label"].nunique()
    conflicts = set(counts[counts > 1].index)
    heldout = heldout.loc[~heldout["canonical_smiles"].isin(conflicts)]
    heldout = heldout.drop_duplicates("canonical_smiles", keep="first")

    lib = library.copy()
    if "canonical_smiles" not in lib:
        raise ValueError("library must contain canonical_smiles")
    lib = lib.dropna(subset=["canonical_smiles"]).drop_duplicates("canonical_smiles", keep="first")
    library_smiles = set(lib["canonical_smiles"].astype(str))
    heldout = heldout.loc[heldout["canonical_smiles"].isin(library_smiles)]
    heldout = heldout.sort_values(["label", "canonical_smiles"], ascending=[False, True])

    labeled_rows = [
        {
            "target_id": target_id,
            "canonical_smiles": row.canonical_smiles,
            "label": int(row.label),
            "label_role": "heldout_positive" if int(row.label) == 1 else "heldout_negative",
        }
        for row in heldout.itertuples(index=False)
    ]
    if len(labeled_rows) > pool_size:
        raise ValueError(f"{target_id} has {len(labeled_rows)} heldout labels for pool size {pool_size}")

    excluded = fit_smiles | {str(row["canonical_smiles"]) for row in labeled_rows}
    background = sorted(library_smiles - excluded)
    needed = pool_size - len(labeled_rows)
    if len(background) < needed:
        raise ValueError(
            f"{target_id} has only {len(background)} eligible background molecules; {needed} required"
        )

    rng = np.random.RandomState(seed + target_index)
    chosen = rng.choice(np.asarray(background, dtype=object), size=needed, replace=False)
    rows = labeled_rows + [
        {
            "target_id": target_id,
            "canonical_smiles": str(smiles),
            "label": None,
            "label_role": "unlabeled_background",
        }
        for smiles in chosen
    ]
    order = rng.permutation(len(rows))
    pool = pd.DataFrame(rows).iloc[order].reset_index(drop=True)
    if len(pool) != pool_size or not pool["canonical_smiles"].is_unique:
        raise AssertionError(f"invalid pool cardinality for {target_id}")
    return pool


def select_targets(
    heldout_rows: pd.DataFrame,
    *,
    n_targets: int,
    fallback_target_ids: Sequence[str] | Iterable[str],
) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    selected_ids: set[str] = set()

    grouped: list[dict[str, object]] = []
    for target_id, rows in heldout_rows.dropna(subset=["target_id"]).groupby("target_id", sort=True):
        n_positive = int(rows.loc[rows["label"] == 1, "canonical_smiles"].nunique())
        n_negative = int(rows.loc[rows["label"] == 0, "canonical_smiles"].nunique())
        if n_positive == 0:
            continue
        tier = "strict_labeled" if n_positive >= 5 and n_negative >= 5 else "limited_positive"
        texts = rows["target_text"].dropna().astype(str)
        grouped.append(
            {
                "target_id": str(target_id),
                "target_text": texts.iloc[0] if not texts.empty else str(target_id),
                "evaluation_tier": tier,
                "n_heldout_positive": n_positive,
                "n_heldout_negative": n_negative,
                "resolution_source": "heldout_label",
            }
        )

    tier_order = {"strict_labeled": 0, "limited_positive": 1}
    grouped.sort(
        key=lambda row: (
            tier_order[str(row["evaluation_tier"])],
            -int(row["n_heldout_positive"]),
            -int(row["n_heldout_negative"]),
            str(row["target_id"]),
        )
    )
    for row in grouped:
        if len(selected) >= n_targets:
            break
        selected.append(row)
        selected_ids.add(str(row["target_id"]))

    for target_id in fallback_target_ids:
        target_id = str(target_id)
        if len(selected) >= n_targets:
            break
        if target_id in selected_ids:
            continue
        selected.append(
            {
                "target_id": target_id,
                "target_text": target_id,
                "evaluation_tier": "unlabeled_fallback",
                "n_heldout_positive": 0,
                "n_heldout_negative": 0,
                "resolution_source": "chembl_fallback",
            }
        )
        selected_ids.add(target_id)

    if len(selected) != n_targets:
        raise ValueError(f"selected {len(selected)} targets, expected {n_targets}")
    return selected
