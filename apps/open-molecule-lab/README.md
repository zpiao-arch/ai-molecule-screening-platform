# Open Molecule Lab

Prompt-first local workbench for the four-level molecule CLI.

This local-first workbench converts a natural-language research request into an auditable `RunSpec`, seals an `id,smiles` MoleculeSet, performs strict preflight, and executes the existing four-level CLI as immutable `prepare -> score -> dock -> report` stage attempts. Missing assets remain explicitly `blocked`; the server never manufactures scores.

## Origin

The interaction shell and visual design were copied from the local `web/evidence-workbench` project. The legacy HelixGuard task runner, generated assets, logs, task state and `node_modules` were not copied. See `NOTICE_THIRD_PARTY.md`.

## Run

```bash
npm ci
npm run build
npm run serve
```

Open `http://127.0.0.1:4173`.

For frontend development, run the API server and Vite separately:

```bash
npm run serve
npm run dev
```

Vite proxies `/api` to the local server on port 4173.

## Verify

```bash
npm run build
npm run contract-check
npm run stage-contract-check
npm run worker-lifecycle-check
```

The contract check verifies health, prompt validation, RunSpec creation, stage/resume HTTP behavior, paged ranked/failed result retrieval, SHA-256 manifest presence and absolute-path exclusion. The stage contract check verifies canonical fingerprints, runtime code identity, immutable attempt history, output tamper detection and contiguous recovery. The worker lifecycle check verifies score/dock separation and process-group cancellation, including a spawned child process.

## Runtime configuration

The source release intentionally omits model weights, receptors and platform binaries. A real run requires an external asset bundle and Python 3.11:

```bash
export OPEN_MOLECULE_PYTHON=/path/to/python3.11
export OPEN_MOLECULE_ASSET_ROOT=/path/to/four-level-molecule-cli-offline-assets
export SMINA_BIN=/path/to/smina
export OBABEL_BIN=/path/to/obabel
npm run serve
```

Without `OPEN_MOLECULE_ASSET_ROOT`, a run is persisted as `blocked` during strict preflight. `queued`, `running`, `complete`, `failed`, `blocked` and `cancelled` are persisted in `run.json`; cancellation and server shutdown signal only the active detached stage process group so docking children do not remain active. Server restart terminates a recorded owned worker before marking its active attempt `failed/worker_interrupted`.

Resume is explicit, not automatic. It revalidates the run manifest, RunSpec, sealed MoleculeSet, asset manifest, runtime code identity, stage fingerprints and output hashes. It then appends a new attempt at the first incomplete stage. A mismatch becomes `blocked/checkpoint_mismatch`; changing code, inputs, assets or policy requires a new run. This contract resumes only at complete stage boundaries and does not restore Python, UniMol, RDKit or smina process memory.

Stage evidence is stored under:

```text
runs/<run_id>/stages/<prepare|score|dock|report>/attempt-0001/
  checkpoint.json
  command.json
  logs/
  outputs/
```

Completed attempts are immutable. Retry/resume creates `attempt-0002` and preserves prior bytes.

## Current API

- `GET /api/health`: local source/CLI discovery and `local_execution` capability.
- `GET /api/model-status`: L1/L2/L3/L4/docking asset availability.
- `POST /api/prompt-plan`: validate a research prompt and persist an auditable plan bundle.
- `POST /api/molecule-sets`: validate and content-address an `id,smiles` CSV (BOM supported).
- `GET /api/molecule-sets/:id`: retrieve immutable MoleculeSet metadata.
- `POST /api/runs`: attach a MoleculeSet, run strict preflight and queue/blocked the CLI execution.
- `GET /api/runs/:id`: retrieve persisted status, preflight and result counts.
- `GET /api/runs/:id/stages`: retrieve persisted attempt counts, states, timestamps and resumability.
- `GET /api/runs/:id/results?view=ranked|failed|all&offset=0&limit=50`: retrieve one bounded result page after completion (`limit` is capped at 200).
- `POST /api/runs/:id/cancel`: request cancellation of a queued/running local worker.
- `POST /api/runs/:id/resume`: revalidate and resume a failed/cancelled run from its next verified stage.

Example request:

```json
{
  "prompt": "Plan an auditable four-level screen for influenza neuraminidase.",
  "target": "CHEMBL2051",
  "candidatePool": 1000,
  "finalSelectionCount": 10,
  "routeMode": "auto"
}
```

## Local verification

Source-only contract, including invalid CSVs and fail-closed preflight:

```bash
npm run contract-check
npm run stage-contract-check
npm run worker-lifecycle-check
```

Two-row real L1-L4 library and smina cascade executions with the offline asset bundle:

```bash
OPEN_MOLECULE_PYTHON=/path/to/python3.11 \
OPEN_MOLECULE_ASSET_ROOT=/path/to/four-level-molecule-cli-offline-assets \
SMINA_BIN=/path/to/smina \
OBABEL_BIN=/path/to/obabel \
npm run real-run-check
```

The real-run checker requires all four layer statuses to be `ok`, runs a library chain and an eight-molecule cascade with real smina rows, interrupts a second cascade after score, and resumes it without changing score attempt 1. The resumed and uninterrupted result CSVs must match column-by-column within absolute tolerance `1e-4`. A separate score-output tamper case must return HTTP 409 `blocked/checkpoint_mismatch`, create no new dock attempt and start no smina process. It also checks complete manifests and rejects host-absolute paths in public evidence.

## Next Boundary

The next product increment is run comparison, evidence bundle download/reproduce, and the 10,000-candidate background scheduler. Those features must reuse the same immutable RunSpec, StageAttempt and manifest contracts; they must not weaken verified resume.
