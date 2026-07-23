# Open Molecule Lab

Prompt-first local workbench for the four-level molecule CLI.

This local-first workbench converts a natural-language research request into an auditable `RunSpec`, seals an `id,smiles` MoleculeSet, performs strict preflight, and can execute the existing four-level CLI in a background worker. Missing assets remain explicitly `blocked`; the server never manufactures scores.

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
```

The contract check verifies health, prompt validation, RunSpec creation, plan bundle contents, SHA-256 manifest presence and absolute-path exclusion.

## Runtime configuration

The source release intentionally omits model weights, receptors and platform binaries. A real run requires an external asset bundle and Python 3.11:

```bash
export OPEN_MOLECULE_PYTHON=/path/to/python3.11
export OPEN_MOLECULE_ASSET_ROOT=/path/to/four-level-molecule-cli-offline-assets
export SMINA_BIN=/path/to/smina
export OBABEL_BIN=/path/to/obabel
npm run serve
```

Without `OPEN_MOLECULE_ASSET_ROOT`, a run is persisted as `blocked` during strict preflight. `queued`, `running`, `complete`, `failed`, `blocked` and `cancelled` are persisted in `run.json`; server restart marks orphaned active runs as `failed/worker_interrupted`. Automatic checkpoint resume is not yet part of this increment.

## Current API

- `GET /api/health`: local source/CLI discovery and plan-only capability.
- `GET /api/model-status`: L1/L2/L3/L4/docking asset availability.
- `POST /api/prompt-plan`: validate a research prompt and persist an auditable plan bundle.
- `POST /api/molecule-sets`: validate and content-address an `id,smiles` CSV (BOM supported).
- `GET /api/molecule-sets/:id`: retrieve immutable MoleculeSet metadata.
- `POST /api/runs`: attach a MoleculeSet, run strict preflight and queue/blocked the CLI execution.
- `GET /api/runs/:id`: retrieve persisted status, preflight and result counts.
- `GET /api/runs/:id/results`: retrieve ranked rows and preserved failed rows after completion.
- `POST /api/runs/:id/cancel`: request cancellation of a queued/running local worker.

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
```

Two-row real L1-L4 execution with the offline asset bundle:

```bash
OPEN_MOLECULE_PYTHON=/path/to/python3.11 \
OPEN_MOLECULE_ASSET_ROOT=/path/to/four-level-molecule-cli-offline-assets \
SMINA_BIN=/path/to/smina \
OBABEL_BIN=/path/to/obabel \
npm run real-run-check
```

The real-run checker requires all four layer statuses to be `ok`, verifies the result summary, checks the complete evidence manifest and rejects host-absolute paths in the public bundle files.

## Next Boundary

The next implementation unit is checkpoint fingerprinting and resume from a verified stage boundary, followed by result comparison and evidence publish metadata. It must preserve the same RunSpec and bundle contracts.
