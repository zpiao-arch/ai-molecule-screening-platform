# Open Molecule Lab

Prompt-first local workbench for the four-level molecule CLI.

This initial prototype converts a natural-language research request into an auditable, persisted `RunSpec`, route decision, stage plan and SHA-256 manifest. It is intentionally `plan_only`: it does not produce molecular scores or docking evidence until a validated molecule-set adapter and strict execution worker are connected.

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

## Current API

- `GET /api/health`: local source/CLI discovery and plan-only capability.
- `GET /api/model-status`: L1/L2/L3/L4/docking asset availability.
- `POST /api/prompt-plan`: validate a research prompt and persist an auditable plan bundle.

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

## Next Boundary

The next implementation unit is a molecule-set adapter plus a strict CLI execution worker. It must validate CSV IDs/SMILES, seal the input hash, run the existing four-level CLI, preserve per-layer statuses, and write real checkpoints. It must not generate placeholder scores.
