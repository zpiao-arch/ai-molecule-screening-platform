#!/usr/bin/env python3
"""Download small public-data caches used by the demo workflow.

The project intentionally does not vendor large chemistry/biology databases.
This helper downloads selected public structure files and writes a manifest so
the same inputs can be reproduced without committing bulky caches to Git.
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Iterable, List


DEFAULT_PDB_IDS = ["3TI6", "3TI5", "3TI3", "6FS6", "6FS8", "2KQT"]


def public_sources() -> List[dict]:
    return [
        {
            "id": "rcsb_pdb",
            "name": "RCSB Protein Data Bank",
            "use": "selected receptor/co-crystal structures, ligands, and structure metadata",
            "url": "https://www.rcsb.org/",
            "api": "https://files.rcsb.org/download/{pdb_id}.pdb",
            "runtime_policy": "download_selected_cache",
        },
        {
            "id": "pubchem",
            "name": "PubChem / PUG REST",
            "use": "compound CID, SMILES, synonyms, and lightweight enrichment",
            "url": "https://pubchem.ncbi.nlm.nih.gov/",
            "runtime_policy": "query_on_demand_or_filtered_cache",
        },
        {
            "id": "chembl",
            "name": "ChEMBL",
            "use": "target-filtered activity records and drug-like molecule metadata",
            "url": "https://chembl.gitbook.io/chembl-interface-documentation/downloads",
            "runtime_policy": "filtered_subset_only",
        },
        {
            "id": "bindingdb",
            "name": "BindingDB",
            "use": "experimental protein-ligand affinity records for selected targets",
            "url": "https://www.bindingdb.org/",
            "runtime_policy": "target_filtered_tsv_only",
        },
        {
            "id": "drugcentral",
            "name": "DrugCentral",
            "use": "known drugs, indications, targets, and structures for seed/control panels",
            "url": "https://drugcentral.org/download",
            "runtime_policy": "small_seed_subset",
        },
        {
            "id": "open_targets",
            "name": "Open Targets Platform",
            "use": "disease-target-drug evidence summaries",
            "url": "https://platform.opentargets.org/",
            "runtime_policy": "filtered_subset_only",
        },
    ]


def fetch_bytes(url: str, timeout: int) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "ai-molecule-screening-platform/0.1"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def write_file(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)


def download_rcsb(pdb_ids: Iterable[str], out: Path, timeout: int, include_cif: bool) -> List[dict]:
    records: List[dict] = []
    rcsb_dir = out / "rcsb"
    for pdb_id in pdb_ids:
        normalized = pdb_id.strip().upper()
        if not normalized:
            continue
        entry = {
            "pdb_id": normalized,
            "pdb": None,
            "cif": None,
            "entry_json": None,
            "status": "pending",
            "errors": [],
        }
        endpoints = [
            ("pdb", f"https://files.rcsb.org/download/{normalized}.pdb", rcsb_dir / f"{normalized}.pdb"),
            ("entry_json", f"https://data.rcsb.org/rest/v1/core/entry/{normalized}", rcsb_dir / f"{normalized}.entry.json"),
        ]
        if include_cif:
            endpoints.append(("cif", f"https://files.rcsb.org/download/{normalized}.cif", rcsb_dir / f"{normalized}.cif"))
        for key, url, path in endpoints:
            try:
                write_file(path, fetch_bytes(url, timeout))
                entry[key] = str(path)
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                entry["errors"].append(f"{key}: {exc}")
        entry["status"] = "ok" if not entry["errors"] else "partial"
        records.append(entry)
    return records


def write_manifest(out: Path, mode: str, pdb_ids: List[str], records: List[dict]) -> Path:
    out.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": "0.1",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "mode": mode,
        "pdb_ids": pdb_ids,
        "public_sources": public_sources(),
        "records": records,
        "policy": [
            "The repository stores small seed metadata only.",
            "Downloaded structure/database caches should stay outside Git or under ignored data/cache paths.",
            "DrugBank and commercial software outputs are not redistributed by this project.",
        ],
    }
    path = out / "manifest.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    return path


def parse_args(argv: List[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download selected public data caches for reproducible demos.")
    parser.add_argument("--out", default="data/public_cache", help="Output cache directory.")
    parser.add_argument("--pdb", action="append", default=[], help="PDB ID to download. Can be repeated.")
    parser.add_argument("--include-cif", action="store_true", help="Also download mmCIF files.")
    parser.add_argument("--timeout", type=int, default=30, help="Network timeout in seconds.")
    parser.add_argument("--dry-run", action="store_true", help="Write manifest only; do not download network files.")
    parser.add_argument("--sources-only", action="store_true", help="Print public source policy JSON and exit.")
    return parser.parse_args(argv)


def main(argv: List[str] | None = None) -> int:
    args = parse_args(argv)
    out = Path(args.out).expanduser().resolve()
    pdb_ids = [item.strip().upper() for item in (args.pdb or DEFAULT_PDB_IDS) if item.strip()]
    if args.sources_only:
        print(json.dumps(public_sources(), ensure_ascii=False, indent=2))
        return 0
    records: List[dict] = []
    mode = "dry_run" if args.dry_run else "download"
    if not args.dry_run:
        records = download_rcsb(pdb_ids, out, args.timeout, args.include_cif)
    manifest = write_manifest(out, mode, pdb_ids, records)
    print(f"wrote manifest: {manifest}")
    if not args.dry_run:
        print(f"downloaded/checked {len(records)} RCSB entries under {out / 'rcsb'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
