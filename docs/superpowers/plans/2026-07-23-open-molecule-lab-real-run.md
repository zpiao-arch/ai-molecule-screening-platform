# Open Molecule Lab Real Run Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the prompt-first Open Molecule Lab prototype into a fail-closed local workflow that seals a molecule CSV, checks the exact Python/assets/binaries required by the resolved route, executes the existing four-level CLI in a background worker, and exposes persisted status plus real result summaries.

**Architecture:** Keep the Python CLI as the only scientific computation boundary. The Node server owns immutable molecule sets, RunSpec attachment, preflight, process lifecycle and evidence files; a small Python bridge reuses `scoring.read_molecules_csv` and pandas for CSV contracts instead of duplicating parsing in JavaScript. Source-only verification must end in an explicit `blocked` run when model assets are absent, while a separate local integration check uses the existing offline asset bundle and smina environment to prove a real run.

**Tech Stack:** Python 3.11, existing four-level CLI/RDKit/pandas, Node.js ESM HTTP server, React 18, TypeScript, Vite, JSON/JSONL evidence files, SHA-256 manifests.

---

## File Responsibility Map

- `scoring/asset_paths.py`: resolve an explicit external asset root without symlinks.
- `scoring/scoring.py`: accept `--asset-root` and pass resolved L2/L3/L4 paths into existing scorers.
- `scoring/scripts/unimol_scorer.py`: load references and UniMol weights from an injected model directory.
- `scoring/pipeline_router.py`: resolve registry-relative receptor paths against the external asset root.
- `apps/open-molecule-lab/server/python-bridge.py`: strict molecule-set validation and result summarization using the CLI's Python contracts.
- `apps/open-molecule-lab/server/store.mjs`: immutable molecule-set/run storage, JSONL events and SHA-256 manifests.
- `apps/open-molecule-lab/server/preflight.mjs`: Python, module, asset, route and binary checks.
- `apps/open-molecule-lab/server/worker.mjs`: background CLI process, cancellation, terminal-state persistence and logs.
- `apps/open-molecule-lab/server/server.mjs`: HTTP routing and composition only.
- `apps/open-molecule-lab/server/contract-check.mjs`: source-only API/fail-closed contract.
- `apps/open-molecule-lab/server/real-run-check.mjs`: opt-in local asset-backed execution check.
- `apps/open-molecule-lab/src/types.ts`: API contracts for molecule sets and executable runs.
- `apps/open-molecule-lab/src/App.tsx`: CSV attachment, run launch, polling and result summary.

### Task 1: Add an explicit external asset root to the scientific CLI

**Files:**
- Create: `scoring/asset_paths.py`
- Modify: `scoring/scoring.py`
- Modify: `scoring/scripts/unimol_scorer.py`
- Modify: `scoring/pipeline_router.py`
- Test: `tests/test_four_level_cli_contract.py`

- [ ] **Step 1: Write failing tests for injected model and receptor roots**

Add tests that create `asset_root/scoring/models/...` and `asset_root/scoring/receptors/...`, then assert:

```python
paths = resolve_asset_paths(asset_root)
assert paths.model_root == asset_root / "scoring" / "models"
assert paths.l2_model.name == "l2_model_sklearn_1_7_2.joblib"
assert paths.admet_model_dir == asset_root / "scoring" / "models" / "admet"
assert paths.unimol_model_dir == asset_root / "scoring" / "models"
```

Patch `FOUR_LEVEL_ASSET_ROOT` and verify `pipeline_router.lookup_receptor("CHEMBL2051")` returns the receptor under `asset_root/scoring/receptors`, never the source directory.

- [ ] **Step 2: Run the focused tests and confirm they fail**

Run:

```bash
/Users/ruiny_park/Documents/药物分子证明/.venv_four_level_cli/bin/python -m pytest -q tests/test_four_level_cli_contract.py -k 'asset_root or receptor_root'
```

Expected: failures because `scoring.asset_paths` and `--asset-root` do not yet exist.

- [ ] **Step 3: Implement the asset-path contract**

Create an immutable `AssetPaths` dataclass with `root`, `manifest`, `model_root`, `l2_model`, `admet_model_dir`, `unimol_model_dir` and `receptor_root`. Resolve the explicit argument first, then `FOUR_LEVEL_ASSET_ROOT`, and return `None` only when neither is set. Prefer `l2_model_sklearn_1_7_2.joblib` when present and fall back to `l2_model.joblib`.

- [ ] **Step 4: Thread the paths through existing scorers**

Add `asset_root: str | Path | None = None` to `MoleculeScorer`, pass the injected L2 path to `Layer2BindingDB`, pass `model_dir` to `Layer3Scorer`, and construct `UniMolScorer(model_dir=...)`. Add CLI `--asset-root`; when provided, set `FOUR_LEVEL_ASSET_ROOT` and `FOUR_LEVEL_ASSET_MANIFEST=<root>/ASSET_MANIFEST.json` before the router or any model is loaded.

- [ ] **Step 5: Make UniMol and receptor resolution injectable**

Change `UniMolScorer.__init__` to store `self.model_dir` and derive reference caches, `drugbank_ref.txt` and weight paths from it. In `pipeline_router._resolved_entry`, resolve relative receptor entries against `<FOUR_LEVEL_ASSET_ROOT>/scoring` when that file exists; otherwise preserve the current registry-relative behavior.

- [ ] **Step 6: Run all Python source tests**

Run:

```bash
/Users/ruiny_park/Documents/药物分子证明/.venv_four_level_cli/bin/python -m pytest -q
```

Expected: `92 passed, 10 deselected` or a larger passing count after the new tests.

### Task 2: Implement immutable MoleculeSet ingestion

**Files:**
- Create: `apps/open-molecule-lab/server/python-bridge.py`
- Create: `apps/open-molecule-lab/server/store.mjs`
- Modify: `apps/open-molecule-lab/server/server.mjs`
- Test: `apps/open-molecule-lab/server/contract-check.mjs`

- [ ] **Step 1: Add failing API cases**

Extend the contract check to POST `/api/molecule-sets` with `{name,csvText,license}`. Assert a BOM-prefixed valid CSV returns `201`, `nRows`, columns, `inputSha256` and `moleculeSetId`; missing columns, duplicate IDs, blank values, empty datasets and payloads above the configured byte limit return structured `400`/`413` errors and do not leave a molecule-set directory.

- [ ] **Step 2: Run the contract check and confirm the endpoint is missing**

Run:

```bash
npm run contract-check
```

Expected: HTTP 404 for `/api/molecule-sets`.

- [ ] **Step 3: Implement the Python CSV bridge**

The `validate-molecule-set` subcommand imports `read_molecules_csv` from the source tree, rejects zero rows or more than 100,000 rows, and prints one JSON object containing `ok`, `nRows`, `columns:["id","smiles"]` and a short sample. It must write errors to JSON without printing model output.

- [ ] **Step 4: Implement immutable storage**

`store.mjs` writes uploads to a staging file, invokes the bridge, derives `molset_<first 16 sha256 chars>`, and atomically renames the validated directory to `molecule_sets/<id>`. Store `input.csv`, `metadata.json` and `MANIFEST.sha256`; a repeated byte-identical upload returns the existing ID, while no API may overwrite content under an existing ID.

- [ ] **Step 5: Add the endpoint and body limits**

Keep JSON upload for the local MVP: the browser reads CSV text and sends it as JSON. Set a configurable maximum (default 8 MiB), return `413 payload_too_large` before JSON parsing, and expose `POST /api/molecule-sets` plus `GET /api/molecule-sets/:id`.

- [ ] **Step 6: Re-run the contract check**

Expected: all valid/invalid/BOM/idempotency cases pass and invalid uploads leave no artifacts.

### Task 3: Add strict preflight and immutable executable RunSpec creation

**Files:**
- Create: `apps/open-molecule-lab/server/preflight.mjs`
- Modify: `apps/open-molecule-lab/server/store.mjs`
- Modify: `apps/open-molecule-lab/server/server.mjs`
- Test: `apps/open-molecule-lab/server/contract-check.mjs`

- [ ] **Step 1: Add a fail-closed run contract test**

Create a plan, attach a valid molecule set, POST `/api/runs`, and assert source-only execution creates a new `run_<timestamp>_<suffix>` with status `blocked`, a persisted attached RunSpec, `preflight.json`, `events.jsonl`, `run.json` and `MANIFEST.sha256`. Assert no `scores.csv` exists and the response names each missing backend/asset.

- [ ] **Step 2: Implement environment discovery**

Resolve `OPEN_MOLECULE_PYTHON` first, then `process.execPath` is not valid for Python, then `python3`. Require Python `3.11.x`; probe imports for `rdkit`, `numpy`, `pandas`, `pyarrow`, `sklearn`, `torch` and `unimol_tools`. Resolve `OPEN_MOLECULE_ASSET_ROOT`, `SMINA_BIN` and `OBABEL_BIN` without searching arbitrary source directories.

- [ ] **Step 3: Verify the external asset manifest**

Invoke `scripts/verify_assets.py --asset-root <root> --manifest <root>/ASSET_MANIFEST.json`. Parse its JSON even on non-zero exit. A missing manifest, missing listed file or SHA mismatch produces a named failed check and blocks strict execution.

- [ ] **Step 4: Make route-aware checks**

Library runs require Python plus L1-L4 assets. Cascade runs additionally require a registered receptor resolved under the asset root and executable smina/obabel. `auto` uses the plan's already persisted resolved branch; preflight must not silently switch a requested/recorded branch.

- [ ] **Step 5: Persist the attached RunSpec before execution**

Read the plan bundle and molecule metadata, require `nRows == expectedCandidateCount`, replace `moleculeSet.attached:false` with the immutable ID/hash/count, set `mode:"execute"`, and hash the complete spec. If preflight fails, emit `run_created` and `preflight_blocked` JSONL events and finish in `blocked` without spawning a process.

- [ ] **Step 6: Verify source-only behavior**

Run `npm run contract-check` with no asset-root variables. Expected: the API is healthy, the executable run is persisted as blocked, and no fake scores are generated.

### Task 4: Execute and supervise the real CLI worker

**Files:**
- Create: `apps/open-molecule-lab/server/worker.mjs`
- Modify: `apps/open-molecule-lab/server/server.mjs`
- Modify: `apps/open-molecule-lab/server/store.mjs`
- Test: `apps/open-molecule-lab/server/real-run-check.mjs`

- [ ] **Step 1: Define the exact worker command**

Build an argument array, never a shell string:

```text
<python> <sourceRoot>/scoring/scoring.py
  --input <runRoot>/inputs/molecules.csv
  --output <runRoot>/results/scores.csv
  --target <targetId>
  --mode <resolvedBranch>
  --strict-backends
  --asset-root <assetRoot>
```

Add `--cascade-top-n` only for cascade runs. Pass `FOUR_LEVEL_ASSET_ROOT`, `FOUR_LEVEL_ASSET_MANIFEST`, `FOUR_LEVEL_RECEPTOR_REGISTRY`, `SMINA_BIN` and `OBABEL_BIN` in a controlled child environment.

- [ ] **Step 2: Persist lifecycle evidence**

Before spawn, write status `queued`; on spawn write `running`, `startedAt`, PID and `worker_started`; pipe stdout/stderr to `logs/stdout.log` and `logs/stderr.log`. On exit 0 require a non-empty scores CSV, summarize it through the Python bridge, then write `complete`; on non-zero exit write `failed` with exit code and a bounded stderr tail. Every transition appends one JSONL event and refreshes `MANIFEST.sha256`.

- [ ] **Step 3: Add status, results and cancellation endpoints**

Implement `GET /api/runs/:id`, `GET /api/runs/:id/results` and `POST /api/runs/:id/cancel`. Cancellation sends SIGTERM only to an in-memory process owned by the server, emits `cancel_requested`, and ends as `cancelled`; completed/failed/blocked runs are immutable and return HTTP 409 on cancellation.

- [ ] **Step 4: Add restart recovery**

At server start, scan runs whose persisted state is `queued` or `running`. If their PID is not an owned live child of this server, mark them `failed` with `worker_interrupted` and retain logs/checkpoints. Do not auto-resume scientific execution until checkpoint fingerprints are implemented.

- [ ] **Step 5: Add the opt-in real-run checker**

`real-run-check.mjs` uploads a two-row CSV, plans `CHEMBL2051` in `library` mode, launches a run, polls terminal state with a bounded timeout, asserts `complete`, asserts every result row preserves `layer1_status` through `layer4_status`, and verifies the run manifest. The script refuses to start unless `OPEN_MOLECULE_PYTHON` and `OPEN_MOLECULE_ASSET_ROOT` are explicitly set.

- [ ] **Step 6: Run the local asset-backed check**

Run:

```bash
OPEN_MOLECULE_PYTHON=/Users/ruiny_park/Documents/药物分子证明/.venv_four_level_cli/bin/python \
OPEN_MOLECULE_ASSET_ROOT=/Users/ruiny_park/Documents/药物分子证明/open_source_release/four-level-molecule-cli-offline-assets \
SMINA_BIN=/Users/ruiny_park/Documents/药物分子证明/ai_drug_eval_tools/micromamba-root/envs/smina-local/bin/smina \
OBABEL_BIN=/Users/ruiny_park/Documents/药物分子证明/ai_drug_eval_tools/micromamba-root/envs/smina-local/bin/obabel \
npm run real-run-check
```

Expected: terminal status `complete`, two real scored rows, no placeholder values, manifest verification true.

### Task 5: Connect the frontend to molecule sets and executable runs

**Files:**
- Modify: `apps/open-molecule-lab/src/types.ts`
- Modify: `apps/open-molecule-lab/src/App.tsx`
- Modify: `apps/open-molecule-lab/src/styles.css`

- [ ] **Step 1: Add typed API contracts**

Define `MoleculeSet`, `ExecutionRun`, `PreflightCheck`, `RunEvent`, `ResultSummary` and terminal statuses `complete|failed|blocked|cancelled`. Preserve `plan_only` types for the planning stage rather than widening every field to `string`.

- [ ] **Step 2: Add CSV attachment to the prompt workflow**

Use a file input limited to `.csv`; read text in the browser, POST `/api/molecule-sets`, and show immutable ID, SHA prefix, row count and validation error. Do not infer scientific validity from the filename or row count.

- [ ] **Step 3: Add a route-aware run action**

Enable the run button only when a plan and molecule set exist and their row counts match. Label it `运行四级 CLI`; a blocked cascade plan remains non-runnable until its route assets are available. Keep `生成执行计划` as a distinct action.

- [ ] **Step 4: Poll and render persisted status**

While `queued` or `running`, poll `GET /api/runs/:id` every second. Render a compact stage/status timeline, preflight failures, elapsed time and a cancel icon button with a tooltip. Never turn `blocked` or `failed` into an empty-success view.

- [ ] **Step 5: Render real result summaries**

On `complete`, fetch `/results` and show ranked valid rows separately from failed/non-evaluable rows. Display ID, four layer statuses, layer scores, final score, gate status and gate reason. Update the scientific boundary from `plan_only` to `real CLI result` only for that completed run.

- [ ] **Step 6: Build and inspect responsive output**

Run `npm run build`, start the local server and inspect desktop/mobile layouts. Expected: no overlapping controls, file labels wrap, stage rows remain stable, and no nested decorative cards are introduced.

### Task 6: Final verification, documentation, archive and GitHub synchronization

**Files:**
- Modify: `apps/open-molecule-lab/README.md`
- Modify: `apps/open-molecule-lab/DESIGN.md`
- Modify: `README.md`
- Modify: `FULL_SOURCE_INVENTORY.md`
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Document operational boundaries**

Describe the JSON CSV upload limit, immutable storage layout, required environment variables, strict blocked behavior, local real-run command, cancellation semantics and the fact that automatic checkpoint resume remains outside this increment.

- [ ] **Step 2: Run all fresh verification commands**

Run Python tests, `npm run build`, `npm run contract-check`, `npm run real-run-check`, `python scripts/verify_assets.py` against the offline root, and `four-level-verify-snapshot` against the shipped compact snapshot. Record exact pass counts and terminal run ID.

- [ ] **Step 3: Rebuild the deterministic source archive**

Use `scripts/build_source_release.py` to create a new dated ZIP under `packages/`. Confirm it excludes `node_modules`, `dist`, `runs`, offline weights, receptors and local binaries; confirm its uncompressed source size and SHA-256.

- [ ] **Step 4: Replace the GitHub main tree**

Use the existing GitHub API channel to create blobs/tree/commit from the clean source archive and update `refs/heads/main` without preserving obsolete source paths. Verify `ai_mol_loop/` count is zero and compare every archive blob path/hash to the remote tree.

- [ ] **Step 5: Publish the archive asset**

Upload the new ZIP, manifest and SHA256SUMS to a dated Release. Download the Release asset again and require byte-identical comparison with the local ZIP.

- [ ] **Step 6: Handle the protected workflow path explicitly**

Before claiming exact tree parity, require `gh auth status` to show `workflow`. If it is missing, report `.github/workflows/ci.yml` as the only protected-path blocker and do not claim the remote tree is byte-for-byte complete. The ordinary source replacement and Release upload may still proceed, but the gap remains visible until `gh auth refresh -h github.com -s workflow` succeeds.

## Plan Self-Review

- The plan covers the selected stage-1 increment: molecule-set attachment, strict preflight, real local worker, cancellation, persisted status, result summary and evidence files.
- It deliberately does not claim checkpoint resume, run comparison, RO-Crate export or 10M interactive execution; those remain later increments from the approved product specification.
- Source-only CI and local scientific integration are separated so that missing proprietary/external assets cannot be replaced with mock scores.
- Public archives continue to omit weights, receptor files, binaries, `node_modules`, build output and run artifacts.
