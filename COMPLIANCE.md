# Release Compliance Notice

## Current license status

This source package is **source-available, not open source** under the current `LICENSE`. The notice does not grant copying, modification, redistribution, sublicensing, or commercial-use rights. Only the copyright holder or an explicitly authorized distributor should publish or share the archive. Do not add an OSI license header or call this repository open source until the rights holder chooses and applies one.

## External assets

The companion offline archive contains BindingDB L2 weights, ADMET weights, Uni-Mol weights, FDA reference embeddings, aligned ChEMBL/BindingDB data, library SMILES, and receptor files. Those assets remain governed by their upstream licenses and terms. smina and obabel are platform-specific executables and are intentionally not redistributed. See [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md) for the public archive data boundary.

Before using external assets, run:

```bash
python scripts/verify_assets.py --asset-root . --manifest ASSET_MANIFEST.json
```

The formal runtime also requires executable `smina` and `obabel` paths through `DOCK_BIN_DIR`, `SMINA_BIN`, or `OBABEL_BIN`. `runtime_doctor --strict` checks the same paths and verifies the asset manifest before probing model backends. An explicit `FOUR_LEVEL_ALLOW_UNVERIFIED_ASSETS=1` is unsafe and causes strict doctor failure.

## Distribution exclusions

The GitHub source archive excludes model/data/receptor assets, local benchmark data, 10M candidate partitions, generated outputs, logs, bytecode, archives, host-specific binaries, and machine-readable row-level data under `legacy/`. The release builder checks these exclusions while preserving the compact frozen-run evidence and historical source documentation.

## Verification scope

Default tests are source-only and do not imply that model assets or docking binaries are present. Asset integration and docking integration are explicit pytest markers. Full local verification must run all tests with the companion assets and `DOCK_BIN_DIR` configured, followed by `verify_snapshot` for the GitHub evidence and the strict `verify_run` command for the complete local run.

The 1000 x 10000 result is pair-heldout computational evidence. It is not a cold-start, temporal, scaffold-external, or wet-lab validation claim.
