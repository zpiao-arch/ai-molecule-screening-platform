# Four-Level CLI Release Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the confirmed CLI contract, provenance, portable-test, and UTF-8 archive defects while preserving the frozen scientific run.

**Architecture:** Keep behavior changes local to `scoring.py`, provenance/report verification, pytest configuration, and one release builder. Regression tests exercise public CLI behavior and verifier outputs; integration-only dependencies are marked explicitly instead of weakening their assertions.

**Tech Stack:** Python 3.11, argparse, pytest, pandas, pathlib, JSONL, zipfile.

---

### Task 1: Cascade Option Semantics

**Files:**
- Modify: `tests/test_four_level_cli_contract.py`
- Modify: `scoring/scoring.py`

- [ ] Add a CLI regression test that invokes a mocked cascade with `--dock-fusion 0` and asserts the fusion function is not called.
- [ ] Add a parser regression test asserting `--dock-fusion -0.1` and `--dock-fusion 1.1` exit with an interval error.
- [ ] Add a CLI regression test asserting `--mode library` with a valid receptor and box exits before scorer or reranker construction.
- [ ] Run the three tests and confirm they fail because zero becomes `0.30`, ranges are accepted, and library mode reaches scoring.
- [ ] Implement an argparse unit-interval type, use `None` as the fusion default, apply `0.30` only when omitted on a cascade branch, and reject library/docking argument conflicts.
- [ ] Run the three tests and the existing cascade contract tests until green.

### Task 2: L2 Unavailable State And Ranking Display

**Files:**
- Modify: `tests/test_four_level_cli_contract.py`
- Modify: `scoring/scoring.py`

- [ ] Replace the legacy fallback expectation with a failing test that forbids DeepPurpose construction and asserts an explicit failed BindingDB-unavailable row.
- [ ] Add a failing output test asserting `print_ranking` displays `final_score_dock` and a cascade-score heading when present.
- [ ] Run both tests and confirm the current DeepPurpose fallback and base-score display cause the failures.
- [ ] Add a minimal BindingDB-unavailable scorer implementing the normal BindingDB score signature and keep `l2_method="bindingdb"` after load failure.
- [ ] Select the ranking score key from the result schema and display the selected value.
- [ ] Run the targeted tests and all CLI contract tests.

### Task 3: Real Stage Log Evidence

**Files:**
- Modify: `tests/test_four_level_run_verifier.py`
- Modify: `scientific_validation/four_level_cli_1kx10k/provenance.py`
- Modify: `scientific_validation/four_level_cli_1kx10k/report.py`
- Modify: `scientific_validation/four_level_cli_1kx10k/verify_run.py`

- [ ] Add a failing test asserting `prepare_delivery_artifacts` does not create `run.log` in an otherwise logless run.
- [ ] Add a failing test asserting stage checkpoint writes append parseable JSONL events to `run.log`.
- [ ] Add a failing test asserting placeholder checkpoints and a placeholder run log do not pass strict provenance.
- [ ] Add a compatibility test constructing a legacy run index plus a nonempty stage log and complete checkpoints, then asserting provenance passes.
- [ ] Run the tests and confirm the current placeholder creation/existence-only checks fail them.
- [ ] Append structured events from `write_stage_checkpoint`, remove report-time placeholder generation, and validate checkpoint contents plus current or legacy log evidence.
- [ ] Run verifier tests and the frozen-run strict verifier.

### Task 4: Portable Test Layers

**Files:**
- Create: `pytest.ini`
- Modify: `tests/test_four_level_cli_contract.py`
- Modify: `tests/test_four_level_dataset.py`
- Modify: `README.md`
- Modify: `VALIDATION.md`

- [ ] Mark tests that require external model/data assets with `integration_assets` and tests that require actual docking binaries with `integration_docking`.
- [ ] Configure default pytest to exclude integration markers and register both markers.
- [ ] Document the source-only, asset integration, and full docking verification commands with accurate scope.
- [ ] Run default pytest in the source tree with model/data/receptor assets absent and no `DOCK_BIN_DIR`; require zero failures.
- [ ] Run the complete suite with the offline assets merged and `DOCK_BIN_DIR` configured; require zero failures.

### Task 5: UTF-8-Safe Release Builder

**Files:**
- Create: `scripts/build_source_release.py`
- Modify: `tests/test_four_level_cli_contract.py`
- Modify: `README.md`

- [ ] Add a failing test importing `build_source_release`, building a temporary archive containing a Chinese filename, and asserting the original name and UTF-8 flag round-trip.
- [ ] Run the test and confirm it fails because the builder does not exist.
- [ ] Implement deterministic sorted ZIP construction with fixed timestamps, preserved Unix modes, UTF-8 names, and exclusions for caches, logs, outputs, local assets, large data, and nested archives.
- [ ] Run the builder test and CLI contract suite.
- [ ] Document the release-build command and the exclusions.

### Task 6: Compliance And Final Release Verification

**Files:**
- Create: `COMPLIANCE.md`
- Modify: `README.md`
- Create outside source tree: versioned ZIP, SHA256SUMS, and release manifest under `packages/`

- [ ] Document that the current package is source-available only under the restrictive license and that third-party assets/binaries are excluded and separately governed.
- [ ] Scan the source tree and archive for model weights, receptor files, host absolute paths, caches, logs, generated outputs, and stale field names.
- [ ] Build the new UTF-8-safe source ZIP and versioned checksum/manifest files; retain the separately verified offline asset archive.
- [ ] Extract the new ZIP into a clean directory and run default source-only tests.
- [ ] Merge verified offline assets, configure `DOCK_BIN_DIR`, run all tests, runtime doctor, and README CLI smoke.
- [ ] Run actual cascade smoke cases for explicit zero fusion, automatic default fusion, and rejected library/docking conflict.
- [ ] Run the strict frozen 1000 x 10000 verifier read-only and record its exact cardinality, manifest, contamination, docking, provenance, and summary results.
- [ ] Confirm Python `zipfile` reports the UTF-8 flag on every non-ASCII entry and no archive entry contains traversal, symlinks, or packaging junk.

No Git commit steps are included because this workspace is not a Git repository. Test output and versioned release manifests serve as execution checkpoints.
