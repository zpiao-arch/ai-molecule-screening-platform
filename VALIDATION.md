# Validation Record

The frozen local run was executed on 2026-07-20/21 from the design at `docs/2026-07-17-four-level-cli-1000x10000-design.md`.

## Engineering gates

| Gate | Result |
|---|---:|
| Targets | 1000 |
| Candidates per target | 10000 |
| Candidate score rows | 10000000 |
| Score partitions | 1000 |
| Unique cached molecules | 98688 |
| L1/L2/L3/L4 failures | 0 / 0 / 0 / 0 |
| Fit-pair contamination | 0 |
| Batch/scalar parity max difference | 0 |
| Compact snapshot manifest files | 19 |
| Historical contract/verifier tests | 43 passed |

The complete local run's strict verifier reported `ok=True`, `manifest.ok=True`, `provenance.ok=True`, and `summary.ok=True`. The GitHub archive contains a separate 19-file compact snapshot; run `python -m scientific_validation.four_level_cli_1kx10k.verify_snapshot --snapshot-dir validation/frozen_run` for that evidence type. It intentionally does not claim to verify omitted row-level score partitions.

## 2026-07-21 release hardening verification

The source-only suite and the asset/docking integration suite are run separately because the GitHub archive intentionally omits model weights and platform binaries. The strengthened verifier is used read-only against the same full 1,000-target run: it recomputes 10,000,000 score rows, score ranges/formulas, zero fit-pair contamination, parity, L2 top-300 docking selection, failed-row preservation, successful-row LE fusion, the complete-run manifest, provenance and the byte-identical summary. The compact snapshot verifier checks only the evidence files shipped in this archive. These are verification passes over frozen artifacts, not relabeling or reruns of model scores.

The release hardening regression suite additionally verifies explicit zero fusion, automatic `0.30` cascade fusion, rejection of library/docking conflicts, explicit BindingDB-unavailable rows, fused-score ranking display, real checkpoint JSONL events, placeholder provenance rejection, portable test markers, and UTF-8 archive names.

## Scientific interpretation

## Open Molecule Lab verified-resume validation (2026-07-23)

The source-only API contract passed CSV BOM/duplicate/blank/oversize cases, exact input/result identity checks, fused-score ranked/failed pagination, invalid pagination rejection, manifest-tamper rejection, stage summary/resume status codes, and fail-closed resume preflight. The stage contract check covered canonical fingerprints, runtime code identity, immutable attempt numbering, terminal-state protection, store-computed output hashes, contiguous stage selection and output tamper detection. The worker lifecycle check rejected an unrelated persisted PID, terminated a worker and its spawned child process, proved cascade score/dock command separation, proved library dock skipping, and waited for both attempt and run cancellation persistence.

With the local Python 3.11 environment, companion offline asset manifest (18 files), and explicit smina/obabel paths, `real-run-check` produced these observed chains:

| Case | Stage evidence |
|---|---|
| Two-row library | prepare 1 complete; score 1 complete; dock skipped; report 1 complete |
| Eight-row uninterrupted cascade | prepare/score/dock/report each complete on attempt 1 |
| Eight-row interrupted and resumed cascade | score remained attempt 1 byte-identical; cancelled dock attempt 1 was preserved; dock attempt 2 and report attempt 1 completed |
| Tampered score output | resume returned HTTP 409 `blocked/checkpoint_mismatch` at score; no new dock attempt; no run-owned smina process |

The resumed and uninterrupted cascade `results/scores.csv` files matched for all columns after ID ordering; numeric comparison used absolute tolerance `1e-4` and reported no mismatches. Both exposed `final_score_dock`, retained the base score, and contained real `structure_docking_status=ok` rows. Completed manifests were rechecked, and all public JSON/JSONL/log/Markdown evidence was scanned for source-root, asset-root and run-root path leakage. Desktop 1440×1000 and mobile 390×844 browser checks found four persisted stage rows, no horizontal overflow, no stage overlap and no button-text overflow.

This is an engineering integration check, not a new scientific benchmark or efficacy validation. The smoke fixtures do not support performance claims.

The strict labeled tier contains 167 targets. L2 AUC median is 0.870833 and four-level AUC median is 0.850694. Top-1% recall medians are both 0; means are 0.019064 (L2) and 0.040870 (four-level). These are pair-heldout measurements with an unlabeled background, not confirmed-negative or cold-start estimates.

Only CHEMBL2051 has a registered receptor. Its L2 top-300 docking produced 157 `ok`, 82 `skipped_hac`, 39 `dock_timeout`, 21 `dock_no_score`, and 1 `prep_timeout`. Fusion AUC was 0.844444, equal to the pre-docking L2 AUC. The historical 0.935 value is not a generalization claim for this run.

## Reproduction limits

The source archive does not contain model weights, source databases, receptor files, binaries, or the 10M candidate-level output. Use the companion offline asset package and the upstream licenses/terms before running a full benchmark.
