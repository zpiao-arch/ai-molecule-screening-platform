# Four-Level CLI Release Hardening Design

## Goal

Remove the confirmed CLI and release-engineering defects without changing model weights, scoring formulas, benchmark labels, frozen candidate scores, or the scientific interpretation of the 1000 x 10000 run.

## Scope

The change covers six boundaries:

1. CLI option semantics for cascade fusion and forced library mode.
2. Explicit non-strict behavior when the packaged BindingDB L2 asset is unavailable.
3. Ranking output that matches the score used for sorting.
4. Provenance evidence backed by stage events instead of a fabricated report-time log.
5. Portable default tests separated from asset and docking integration tests.
6. A deterministic UTF-8-safe source ZIP and release verification artifacts.

The current restrictive `LICENSE` remains authoritative until the copyright holder explicitly selects an open-source license. Documentation and release manifests must describe the archive as source-available, not open source, while that notice remains in place.

## CLI Contract

`--dock-fusion` uses `None` as the parser default. An omitted value receives the automatic cascade default of `0.30`; an explicit `0` remains zero and disables score fusion. Values outside `[0, 1]` are rejected by argparse.

`--mode library` rejects receptor, box, flexible-residue, protonation, reranking, and explicitly supplied fusion options. This prevents a command from announcing the library branch and then invoking smina.

When BindingDB cannot load in non-strict mode, the scorer retains the `bindingdb` method and uses an explicit unavailable backend result. Each row records a failed L2 status, backend name, and zero L2 contribution without importing the unbundled DeepPurpose implementation. Strict mode continues to fail at initialization.

When `final_score_dock` exists, the main ranking table displays that value and labels the column as the cascade score. The base `final_score` remains preserved in CSV output.

## Provenance

Every call to `write_stage_checkpoint` appends one JSON object to `run.log` after the checkpoint is atomically written. The event records timestamp, stage, status, input fingerprint, and checkpoint-relative path.

The report stage no longer creates a narrative placeholder `run.log`. Strict provenance requires all checkpoint files plus one of two evidence forms:

- current format: at least one valid stage-checkpoint JSONL event matching an archived checkpoint;
- frozen-run compatibility: a nonempty legacy run index plus at least one nonempty archived stage stdout/stderr log.

This keeps the sealed 2026-07-20 run verifiable without allowing a new directory containing only placeholder files to pass.

## Test And Packaging Contract

Default `pytest` runs only tests that can pass from the GitHub source archive. Tests requiring model/data assets use `integration_assets`; tests requiring smina/obabel use `integration_docking`. Full local verification overrides the default marker filter after assets and binaries are configured.

The release builder walks files in sorted order, excludes caches, generated logs, outputs, local data/assets, and archives, then writes a deterministic ZIP. Non-ASCII archive names must carry the ZIP UTF-8 flag and round-trip through Python `zipfile` with their original names.

## Verification

Completion requires all of the following fresh evidence:

- targeted regression tests pass after having failed against the old behavior;
- default source-only tests pass without external assets or docking binaries;
- all tests pass with merged offline assets and `DOCK_BIN_DIR` configured;
- runtime doctor reports all four layers and both docking binaries available;
- actual CLI smoke proves explicit zero fusion and rejects library/docking conflicts;
- strict frozen-run verification still reports 1000 targets, 10,000,000 rows, zero fit-pair contamination, valid docking, manifest, provenance, and byte-identical summary;
- rebuilt ZIP preserves all Chinese filenames, contains no caches/logs/host paths, and passes checksum and extraction checks.
