# Open Molecule Lab Initial Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Copy the local evidence-workbench interaction shell into the four-level CLI source tree and create a prompt-first Open Molecule Lab prototype that turns a research request into an auditable RunSpec and execution plan.

**Architecture:** Reuse the existing local React/Vite evidence-workbench visual shell and its warm evidence-ledger design, but keep a new Open Molecule Lab app directory with no dependency on the legacy HelixGuard task runner. A small Node server will serve the built frontend, expose health/model status/prompt-plan APIs, validate prompt input deterministically, and persist plan-only run bundles under a local runs directory. The server will not claim to execute four-level scoring until a molecule-set input and verified assets are supplied.

**Tech Stack:** React 18, TypeScript, Vite, lucide-react, Node.js ESM HTTP server, JSON RunSpec files, existing Python four-level CLI discovered from the app directory.

---

### Task 1: Copy the evidence-workbench source boundary

**Files:**
- Create: `apps/open-molecule-lab/package.json`
- Create: `apps/open-molecule-lab/package-lock.json`
- Create: `apps/open-molecule-lab/index.html`
- Create: `apps/open-molecule-lab/DESIGN.md`
- Create: `apps/open-molecule-lab/src/App.tsx`
- Create: `apps/open-molecule-lab/src/main.tsx`
- Create: `apps/open-molecule-lab/src/styles.css`
- Create: `apps/open-molecule-lab/src/types.ts`
- Create: `apps/open-molecule-lab/vite.config.ts`
- Create: `apps/open-molecule-lab/tsconfig.json`
- Create: `apps/open-molecule-lab/tsconfig.node.json`
- Create: `apps/open-molecule-lab/NOTICE_THIRD_PARTY.md`

- [ ] Copy only source, configuration, lockfile and design documentation from `web/evidence-workbench`; exclude `node_modules/`, `dist/`, `.task-runner.log`, `.task-runner.pid`, `task_runner_state/` and build caches.
- [ ] Preserve the local source provenance in `NOTICE_THIRD_PARTY.md`, including the original path and the fact that the current package remains source-available.
- [ ] Confirm the copied tree contains no host-specific task-runner path or legacy HelixGuard server file.

Run:
```bash
find apps/open-molecule-lab -type f | sort
```
Expected: source/config files only; no `node_modules`, `dist`, PID, log or task state entries.

### Task 2: Rename and reduce the frontend to a prompt-first Open Molecule Lab shell

**Files:**
- Modify: `apps/open-molecule-lab/package.json`
- Modify: `apps/open-molecule-lab/index.html`
- Modify: `apps/open-molecule-lab/src/App.tsx`
- Modify: `apps/open-molecule-lab/src/styles.css`
- Modify: `apps/open-molecule-lab/src/types.ts`

- [ ] Rename the product-facing strings from HelixGuard/customer screening to Open Molecule Lab/research request.
- [ ] Replace the legacy disease form with a single research prompt textarea plus explicit target, candidate-pool and output-count controls.
- [ ] Keep RDKit as a required baseline indicator, but show L2/L3/L4 and docking as asset-dependent status rows.
- [ ] Make the primary interaction submit a prompt to `/api/prompt-plan`; render the returned normalized intent, route decision, stages, boundaries and equivalent CLI preview.
- [ ] Keep the evidence ledger and audit drawer, but label the result as `plan_only` until a molecule set is attached and a real run is executed.
- [ ] Preserve responsive layout and accessible landmarks from the copied workbench; remove unused imports and legacy result types so TypeScript remains strict.

Run:
```bash
npm run build
```
Expected: Vite and TypeScript build successfully with no HelixGuard strings in `src/`.

### Task 3: Add the local prompt-plan server

**Files:**
- Create: `apps/open-molecule-lab/server/server.mjs`
- Modify: `apps/open-molecule-lab/package.json`
- Create: `apps/open-molecule-lab/README.md`

- [ ] Serve `dist/` and expose `GET /api/health` with product name, source root, CLI path, asset status and `plan_only` capability.
- [ ] Expose `GET /api/model-status` with deterministic rows for RDKit, L2 BindingDB, L3 ADMET, L4 UniMol and docking; status must be `available`, `not_found` or `planned`, never inferred from a fake score.
- [ ] Expose `POST /api/prompt-plan` accepting `{prompt, target, candidatePool, finalSelectionCount, routeMode}`.
- [ ] Validate non-empty prompt, integer bounds (`1..100000` candidates, `1..100` final selections), and route mode (`auto`, `library`, `cascade`). Return HTTP 400 with a structured error for invalid input.
- [ ] Normalize target text, infer a safe target slug only when the target field is empty, and choose `library` unless an explicit cascade route and a registered receptor are both available.
- [ ] Persist a plan-only run bundle at `runs/<run_id>/run.json`, `runs/<run_id>/prompt.txt`, `runs/<run_id>/DESIGN.md` and `runs/<run_id>/MANIFEST.sha256` using relative paths and SHA-256 values.
- [ ] Return `run_id`, normalized request, route rationale, stage list, execution boundary, asset requirements and an equivalent CLI command. Do not invoke scoring or docking in this task.

Run:
```bash
node --check server/server.mjs
```
Expected: exit code 0.

### Task 4: Add prompt API contract tests and smoke scripts

**Files:**
- Create: `apps/open-molecule-lab/server/contract-check.mjs`
- Modify: `apps/open-molecule-lab/package.json`

- [ ] Start the server on an ephemeral local port in the contract script.
- [ ] Assert health reports `open-molecule-lab` and `plan_only`.
- [ ] Assert a valid prompt produces a run ID, a `RunSpec`, a non-empty stage list, and a bundle manifest.
- [ ] Assert an empty prompt, out-of-range candidate count, and invalid route each return HTTP 400 without creating a run.
- [ ] Assert the generated bundle contains no absolute workspace path.

Run:
```bash
npm run contract-check
```
Expected: JSON output with `ok: true` and a generated plan-run ID.

### Task 5: Verify the copied app and document the next integration boundary

**Files:**
- Modify: `apps/open-molecule-lab/README.md`
- Modify: `docs/superpowers/specs/2026-07-22-open-molecule-lab-design.md` only if the initial prototype changes a public contract.

- [ ] Run `npm run build` and `npm run contract-check` from `apps/open-molecule-lab`.
- [ ] Start the preview server and use `curl` to verify `/api/health`, `/api/model-status` and `/api/prompt-plan`.
- [ ] Verify the existing four-level source tests remain unchanged and still pass; do not add model-dependent tests for a plan-only prototype.
- [ ] Document that the next task is a molecule-set adapter and strict CLI execution worker, not a mock score generator.

Run:
```bash
npm run build
npm run contract-check
```
Expected: both commands pass; no claim of real scoring is made by the prototype.
