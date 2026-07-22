# Four-Level CLI Release Hardening V2 Design

## Goal

Close the eight confirmed release defects without changing model weights, scoring formulas, benchmark labels, candidate ordering, or the scientific interpretation of the verified 1000 x 10000 run.

## Release Boundary

The GitHub archive remains source-available under the existing restrictive `LICENSE`. This change does not choose an OSI license on behalf of the copyright holder.

The public archive keeps current CLI code, validation code, tests, documentation, historical Python/shell scripts, aggregate Markdown reports, and figures. Row-level or machine-readable legacy research data (`legacy/**/*.csv`, `legacy/**/*.json`, and `legacy/**/*.parquet`) are excluded until a per-file redistribution review is complete. The local working copy is not deleted; only the deterministic public archive changes.

## Trusted Asset Loading

Add one shared asset-integrity module under `scoring/`. It locates the manifest in this order:

1. `FOUR_LEVEL_ASSET_MANIFEST` when explicitly set;
2. `ASSET_MANIFEST.json` in the merged release root;
3. `assets/ASSET_MANIFEST.json` in the source release root.

Every packaged joblib, pickle, NumPy reference cache, and Uni-Mol weight is verified against the manifest before deserialization. A missing manifest entry, missing file, or SHA-256 mismatch raises an `AssetIntegrityError`. Custom assets require a custom manifest. The only escape hatch is the explicit environment variable `FOUR_LEVEL_ALLOW_UNVERIFIED_ASSETS=1`; using it is reported by runtime doctor and must never be enabled by default or in CI.

Hash results are cached by resolved path, size, modification time, manifest path, and expected digest so batch scoring does not repeatedly hash a model for every molecule.

Runtime doctor uses the same verifier and includes a manifest summary in `overall`. It also resolves `SMINA_BIN`, `OBABEL_BIN`, `DOCK_BIN_DIR`, package-local `bin/`, and `PATH` in the same order as the production docking module.

## Compact Frozen Snapshot

The 1 GB complete run remains untouched and continues to use its complete `MANIFEST.sha256` and strict `verify_run` path.

The GitHub snapshot becomes a separate evidence type:

- its `MANIFEST.sha256` lists only files actually included in `validation/frozen_run/`;
- `snapshot.json` records the compact schema, the complete run manifest digest, expected target/pair counts, and the intentionally omitted artifact classes;
- `verify_snapshot.py` validates the compact manifest and cross-checks `summary.json`, `score_summary.json`, failure evidence, parity evidence, per-target metrics, checkpoints, and docking summaries;
- it never claims to re-check 10 million row-level scores, fit-pair membership, or omitted docking poses.

Documentation must distinguish `verify_snapshot` from the complete-run `verify_run --strict` command.

## Receptor Protonation

Keep `--mcce` for CLI compatibility but describe it accurately as Open Babel pH 7.4 protonation, not a real MCCE run.

Protonation writes to a unique temporary directory. The code removes any pre-existing output before invocation, requires return code zero, requires a non-empty output file, and includes stderr in a bounded error message. The temporary receptor remains alive only for the docking call and is always removed afterward. A failed protonation can never fall back to an adjacent stale `.protonated.pdb` file.

## CLI Contracts And Failure State

All count, batch, worker, timeout, exhaustiveness, CPU, heavy-atom, and top-N options use shared argparse validators. Values that must be positive fail during parsing; seeds are non-negative.

`--dock-mode` and `--cascade-top-n` use `None` parser defaults so the CLI can distinguish omission from explicit use. Any explicit docking-only option is rejected when the effective route is `library`, including an `auto` route that discovers no usable receptor. A valid explicit receptor plus box makes `auto` choose the cascade branch.

CSV input uses `utf-8-sig`, preserving ordinary UTF-8 behavior while accepting an Excel BOM.

Provenance accepts `failed` checkpoints. If a stage raises after writing `in_progress`, the CLI atomically rewrites that checkpoint with the same input fingerprint, error type, bounded message, and failure timestamp, then re-raises. Retrying with matching inputs is permitted; a failed checkpoint never counts as complete.

## Installable CLI And CI

Add a `pyproject.toml` with Python 3.11 metadata and these console scripts:

- `four-level-molecule` -> small-batch scoring CLI;
- `four-level-benchmark` -> staged 1000 x 10000 CLI;
- `four-level-doctor` -> runtime and asset preflight;
- `four-level-verify-snapshot` -> compact evidence verification.

The `scoring` directory becomes an importable package while retaining direct `python scoring/scoring.py` compatibility through dual package/script imports.

GitHub Actions runs the source-only test suite, builds a wheel, installs it without dependencies into the tested environment, exercises all console-script `--help` commands, builds the deterministic source ZIP, and verifies its exclusions. Asset and docking integration remain local because the redistributable GitHub repository intentionally omits those assets and platform binaries.

## Verification Gates

Completion requires fresh evidence for every boundary:

1. Each regression test fails against the old behavior and passes after its focused fix.
2. Default source tests and `pip check` pass.
3. Clean assets pass; a one-byte or expected-hash mismatch fails before deserialization and makes doctor fail.
4. BOM CSV succeeds; invalid numeric arguments fail with argparse exit code 2 before stage work starts.
5. A forced stage exception produces a `failed` checkpoint and failure event.
6. Failed Open Babel protonation cannot reach `DockingReranker` and leaves no receptor artifact.
7. Compact snapshot verification passes while complete-run strict verification still reports 1000 targets, 10,000,000 rows, zero contamination, valid docking, manifest, and provenance.
8. Current release scoring remains numerically identical on the 20-molecule parity fixture and at least one complete 10,000-row target partition.
9. The wheel entry points work, CI syntax is valid, and the rebuilt GitHub ZIP contains no row-level legacy data, caches, logs, host paths, model weights, receptors, or binaries.

