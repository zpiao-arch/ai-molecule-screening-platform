# Data Sources and Local Storage Policy

The product should be offline-first at runtime. Public databases are used to build small, target-specific caches and evidence summaries. Large public or restricted databases are not vendored in this repository.

## Included Seed Data

`ai_mol_loop/targets/influenza/` contains a small public seed set:

- target catalog for influenza NA, PA endonuclease, and M2 examples;
- known-drug/control metadata, including oseltamivir, zanamivir, peramivir, laninamivir, baloxavir, amantadine, and rimantadine contexts;
- representative PDB IDs such as `3TI6`, `3TI5`, `3TI3`, `6FS6`, `6FS8`, and `2KQT`;
- offline evidence summaries generated from public metadata.

These files are small enough for Git and are meant to bootstrap a reproducible demo.

## Public Sources

| Source | Role | Local policy |
|---|---|---|
| RCSB PDB | receptor structures, co-crystal ligands, structure metadata | download selected PDB/mmCIF/JSON files only |
| PubChem | CID, SMILES, synonyms, lightweight compound enrichment | query on demand or store small filtered caches |
| ChEMBL | target-filtered activity and molecule records | store filtered subsets, not full DB |
| BindingDB | experimental affinity records | store target-filtered TSV extracts only |
| DrugCentral | known drugs, indications, targets, structures | store small seed/control subset |
| Open Targets | disease-target-drug evidence | store filtered summaries only |

## Download Selected Public Structures

```bash
python scripts/download_public_data.py --out data/public_cache --pdb 3TI6 --pdb 6FS6
```

Dry run:

```bash
python scripts/download_public_data.py --out data/public_cache --pdb 3TI6 --dry-run
```

The script writes `data/public_cache/manifest.json`. The `data/` directory is ignored by Git.

## Optional External Source Mirrors

To inspect generation/docking projects without adding them to this repository:

```bash
python scripts/download_ai_molecule_repos_api.py
```

The script writes to `data/external_repos/`, which is ignored by Git.

## Restricted or Commercial Data

DrugBank, commercial docking software outputs, vendor datasets, and non-public negative experimental data are not redistributed. If used locally, keep them outside Git and document only the source, license, and import schema.

## Claims Boundary

Evidence source availability, docking readiness, or known-drug control ranking is not biological validation. The platform may prioritize candidates for follow-up, but it does not prove efficacy, safety, toxicity, dose, or clinical benefit.
