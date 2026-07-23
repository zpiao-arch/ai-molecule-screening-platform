# Open Molecule Lab Verified Stage Resume Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split interactive Open Molecule Lab runs into fingerprinted prepare/score/dock/report StageAttempts and allow failed, cancelled or interrupted runs to resume only from verified immutable stage outputs.

**Architecture:** Keep `scoring.py` as the scientific boundary, adding a strict `--base-scores` path so cascade docking can consume verified L1-L4 output without recomputing it. Add a pure stage-contract module, a filesystem StageAttempt store and a sequential worker orchestrator; the Run API exposes stage history and explicit resume, while the React UI renders persisted attempts rather than inferred progress.

**Tech Stack:** Python 3.11, RDKit/pandas/NumPy, Node.js ESM, React 18, TypeScript, Vite, JSON/JSONL, SHA-256, existing smina/Open Babel integration.

---

The source directory is not a Git checkout. Local commit steps are intentionally omitted; every task ends with a named red/green verification command, and the final task rebuilds the deterministic archive and verifies the GitHub tree.

## File Responsibility Map

- `scoring/scoring.py`: load and validate immutable base scores, skip L1-L4 when `--base-scores` is present, and prove docking does not mutate base columns.
- `tests/test_four_level_cli_contract.py`: Python contract and CLI compatibility tests for staged scoring.
- `apps/open-molecule-lab/server/stage-contract.mjs`: canonical JSON, code identity, stage ordering and input fingerprint pure functions.
- `apps/open-molecule-lab/server/stage-store.mjs`: atomic StageAttempt creation/transitions, output hashing and contiguous checkpoint verification.
- `apps/open-molecule-lab/server/python-bridge.py`: stage result identity, ranking field and complete structure-docking success counts.
- `apps/open-molecule-lab/server/stage-contract-check.mjs`: source-only tests for fingerprints, immutable attempts and resume verification.
- `apps/open-molecule-lab/server/store.mjs`: execution-run metadata and manifest integration only.
- `apps/open-molecule-lab/server/worker.mjs`: sequential score/dock/report orchestration and process-group lifecycle.
- `apps/open-molecule-lab/server/server.mjs`: stage/resume endpoints, fresh preflight and restart recovery composition.
- `apps/open-molecule-lab/server/contract-check.mjs`: source-only HTTP contracts, including non-resumable blocked runs.
- `apps/open-molecule-lab/server/real-run-check.mjs`: asset-backed stage-chain, interruption, resume, parity and tamper checks.
- `apps/open-molecule-lab/src/types.ts`: typed StageAttempt and resumable run contracts.
- `apps/open-molecule-lab/src/App.tsx`: persisted stage timeline and explicit resume action.
- `apps/open-molecule-lab/src/styles.css`: stable stage timeline and attempt metadata layout.

## Task 1: Add a Strict Base-Score Scientific Boundary

**Files:**
- Modify: `scoring/scoring.py`
- Modify: `tests/test_four_level_cli_contract.py`

- [ ] **Step 1: Write failing base-score identity and mutation tests**

Add tests that create an input CSV plus base score CSV and assert exact identity, type coercion and base-column immutability:

```python
def test_read_base_scores_requires_exact_input_identity(tmp_path):
    scoring = _load_module_from_file("staged_scoring_identity", SCORING_DIR / "scoring.py")
    expected = [("mol-1", "CCO"), ("mol-2", "CCN")]
    base = tmp_path / "base.csv"
    base.write_text(
        "id,smiles,layer1_status,layer2_status,layer3_status,layer4_status,"
        "docking_normalized,final_score\n"
        "mol-1,CCO,ok,ok,ok,ok,0.2,0.4\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="identity mismatch"):
        scoring.read_base_scores_csv(base, expected)


def test_validate_base_columns_rejects_docking_mutation():
    scoring = _load_module_from_file("staged_scoring_mutation", SCORING_DIR / "scoring.py")
    before = [{"id": "mol-1", "smiles": "CCO", "final_score": 0.4}]
    after = [{"id": "mol-1", "smiles": "CCO", "final_score": 0.5,
              "final_score_dock": 0.6}]
    with pytest.raises(ValueError, match="base score mutated"):
        scoring.validate_base_columns_unchanged(before, after)
```

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
<workspace>/.venv_four_level_cli/bin/python \
  -m pytest -q tests/test_four_level_cli_contract.py \
  -k 'base_scores or base_columns'
```

Expected: fail because `read_base_scores_csv` and `validate_base_columns_unchanged` do not exist.

- [ ] **Step 3: Implement base-score loading and validation**

Add these public helpers near `read_molecules_csv`:

```python
BASE_REQUIRED_COLUMNS = {
    "id", "smiles",
    "layer1_status", "layer2_status", "layer3_status", "layer4_status",
    "docking_normalized", "final_score",
}

DOCKING_APPEND_COLUMNS = {
    "docking_affinity_kcal_mol", "heavy_atoms", "ligand_efficiency",
    "dock_rerank_rank", "structure_docking_status", "dock_le_norm",
    "final_score_dock",
}


def read_base_scores_csv(path: str | Path,
                         expected: List[Tuple[str, str]]) -> List[Dict]:
    import pandas as pd

    frame = pd.read_csv(path, dtype={"id": "string", "smiles": "string"})
    missing = sorted(BASE_REQUIRED_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"base scores missing required columns: {', '.join(missing)}")
    if frame["id"].isna().any() or frame["id"].duplicated().any():
        raise ValueError("base scores contain blank or duplicate ids")
    expected_map = dict(expected)
    actual_map = dict(zip(frame["id"].astype(str), frame["smiles"].astype(str), strict=True))
    if actual_map != expected_map or len(frame) != len(expected):
        raise ValueError("base score identity mismatch")
    rows = frame.where(frame.notna(), None).to_dict(orient="records")
    return [{str(key): value for key, value in row.items()} for row in rows]


def validate_base_columns_unchanged(before: List[Dict], after: List[Dict]) -> None:
    before_by_id = {str(row["id"]): row for row in before}
    after_by_id = {str(row["id"]): row for row in after}
    if before_by_id.keys() != after_by_id.keys():
        raise ValueError("base score mutated: molecule ids changed")
    for mol_id, original in before_by_id.items():
        current = after_by_id[mol_id]
        for key, value in original.items():
            if key in DOCKING_APPEND_COLUMNS:
                continue
            if current.get(key) != value:
                raise ValueError(f"base score mutated: {mol_id}.{key}")
```

- [ ] **Step 4: Add `--base-scores` without changing default CLI behavior**

Add:

```python
parser.add_argument(
    "--base-scores",
    default=None,
    help="已验证的 L1-L4 基础结果 CSV；仅 cascade 阶段使用并跳过 L1-L4 重算",
)
```

Reject it outside cascade, snapshot the loaded rows, and skip scorer construction:

```python
if args.base_scores and args.mode != "cascade":
    parser.error("--base-scores 仅允许与 --mode cascade 一起使用")

def load_or_score_results(args, rows, decision, default_tt):
    if args.base_scores:
        loaded = read_base_scores_csv(args.base_scores, rows)
        return loaded, [dict(row) for row in loaded]
    scorer = MoleculeScorer(
        deeppurpose_target=args.target,
        deeppurpose_target_seq=args.target_seq,
        l2_method=args.l2,
        default_target_text=default_tt,
        l2_model_path=decision.get("l2_model_path"),
        strict_backends=args.strict_backends,
        asset_root=args.asset_root,
    )
    return scorer.score_batch(rows), None


results, base_snapshot = load_or_score_results(args, rows, decision, default_tt)
```

Use `MoleculeScorer.save_csv`/`print_ranking` statically. Before writing a `--base-scores` result, call `validate_base_columns_unchanged(base_snapshot, results)`.

- [ ] **Step 5: Prove L1-L4 initialization is skipped**

Test the helper directly with a complete base fixture:

```python
def test_load_or_score_results_skips_scorer_for_base_scores(tmp_path, monkeypatch):
    from types import SimpleNamespace

    scoring = _load_module_from_file("staged_scoring_skip", SCORING_DIR / "scoring.py")
    base = tmp_path / "base.csv"
    base.write_text(
        "id,smiles,layer1_status,layer2_status,layer3_status,layer4_status,"
        "docking_normalized,final_score\n"
        "mol-1,CCO,ok,ok,ok,ok,0.2,0.4\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        scoring,
        "MoleculeScorer",
        type("ForbiddenScorer", (), {"__init__": lambda self, **kwargs: (_ for _ in ()).throw(AssertionError("scorer initialized"))}),
    )
    args = SimpleNamespace(
        base_scores=str(base), target="CHEMBL2051", target_seq=None,
        l2="bindingdb", strict_backends=True, asset_root=None,
    )
    rows, snapshot = scoring.load_or_score_results(
        args, [("mol-1", "CCO")], {"l2_model_path": None}, "Neuraminidase"
    )
    assert rows == snapshot
```

Also assert `python scoring/scoring.py --help` contains `--base-scores`.

- [ ] **Step 6: Run Python regression suite**

Run:

```bash
<workspace>/.venv_four_level_cli/bin/python -m pytest -q
```

Expected: all source-only tests pass; existing direct library/cascade CLI contracts remain unchanged.

## Task 2: Implement Deterministic Stage Fingerprints

**Files:**
- Create: `apps/open-molecule-lab/server/stage-contract.mjs`
- Create: `apps/open-molecule-lab/server/stage-contract-check.mjs`
- Modify: `apps/open-molecule-lab/package.json`

- [ ] **Step 1: Write failing pure contract checks**

The check must assert stable canonicalization, changed fingerprints, library/cascade ordering and runtime code identity:

```js
import {
  canonicalJson,
  computeCodeIdentity,
  computeStageFingerprint,
  stageSequence,
} from "./stage-contract.mjs";

assert(canonicalJson({ b: 2, a: 1 }) === '{"a":1,"b":2}', "canonical keys changed");
assert(JSON.stringify(stageSequence("library")) === JSON.stringify(["prepare", "score", "report"]), "library order");
assert(JSON.stringify(stageSequence("cascade")) === JSON.stringify(["prepare", "score", "dock", "report"]), "cascade order");
const left = computeStageFingerprint({ stage: "score", runSpecSha256: "a", moleculeSetSha256: "b" });
const right = computeStageFingerprint({ stage: "score", runSpecSha256: "a", moleculeSetSha256: "c" });
assert(left !== right, "molecule hash did not affect fingerprint");
const identity = await computeCodeIdentity(sourceRoot);
assert(identity.sha256.length === 64 && identity.files.some((row) => row.path === "scoring/scoring.py"), "code identity incomplete");
```

- [ ] **Step 2: Run and verify RED**

Run: `node server/stage-contract-check.mjs`  
Expected: module-not-found for `stage-contract.mjs`.

- [ ] **Step 3: Implement canonical JSON and stage ordering**

Export:

```js
export const STAGE_SCHEMA = "open-molecule-lab.stage-attempt.v0.1";

export function canonicalJson(value) {
  if (Array.isArray(value)) return `[${value.map(canonicalJson).join(",")}]`;
  if (value && typeof value === "object") {
    return `{${Object.keys(value).sort().map((key) => `${JSON.stringify(key)}:${canonicalJson(value[key])}`).join(",")}}`;
  }
  return JSON.stringify(value);
}

export function stageSequence(branch) {
  if (branch === "library") return ["prepare", "score", "report"];
  if (branch === "cascade") return ["prepare", "score", "dock", "report"];
  throw new Error(`unsupported branch: ${branch}`);
}

export function computeStageFingerprint(input) {
  return createHash("sha256").update(canonicalJson(input), "utf8").digest("hex");
}
```

- [ ] **Step 4: Implement runtime code identity**

Walk only these inputs: `scoring/**/*.py`, `apps/open-molecule-lab/server/*.{mjs,py}`, `requirements.lock.txt`, `requirements-runtime.txt`, and `apps/open-molecule-lab/package-lock.json`. Return sorted `{path,sha256}` rows plus the canonical list hash. Reject unreadable required roots; never include absolute paths in the returned object.

- [ ] **Step 5: Register and run the check**

Add to `package.json`:

```json
"stage-contract-check": "node server/stage-contract-check.mjs"
```

Run: `npm run stage-contract-check`  
Expected: JSON `{ "ok": true, ... }` and exit 0.

## Task 3: Add an Atomic StageAttempt Store

**Files:**
- Create: `apps/open-molecule-lab/server/stage-store.mjs`
- Modify: `apps/open-molecule-lab/server/stage-contract-check.mjs`
- Modify: `apps/open-molecule-lab/server/store.mjs`

- [ ] **Step 1: Write RED tests for immutable attempts**

Using a temporary run root, assert:

```js
const first = await stages.createAttempt(runId, "score", fingerprint, command);
assert(first.attempt === 1 && first.status === "queued", "first attempt incorrect");
await stages.transition(runId, "score", 1, "running");
await fs.writeFile(path.join(first.root, "outputs", "scores.csv"), "id,smiles\na,CCO\n");
const complete = await stages.complete(runId, "score", 1, { "scores.csv": "outputs/scores.csv" });
assert(complete.status === "complete", "attempt did not complete");
const second = await stages.createAttempt(runId, "score", fingerprint, command);
assert(second.attempt === 2, "attempt number was not incremented");
await assert.rejects(() => stages.transition(runId, "score", 1, "running"), /terminal attempt/);
```

Tamper with attempt 1 output and assert `verifyAttempt` returns `output_mismatch`.

- [ ] **Step 2: Run and verify RED**

Run: `npm run stage-contract-check`  
Expected: import failure for `stage-store.mjs`.

- [ ] **Step 3: Implement atomic StageAttempt operations**

Create `createStageStore({ runsRoot, refreshManifest })` with:

```js
return {
  createAttempt,
  transition,
  complete,
  fail,
  listAttempts,
  summarize,
  verifyAttempt,
  verifyContiguous,
};
```

Write JSON via `checkpoint.json.tmp-<pid>-<random>` followed by `rename`. `createAttempt` creates `attempt-0001`, `logs/`, and `outputs/`; it must fail if the destination already exists. Terminal attempts are immutable. `complete` hashes declared output files itself rather than accepting caller-provided digests.

- [ ] **Step 4: Implement contiguous resume verification**

`verifyContiguous(run, currentInputs)` must iterate `stageSequence(run.route.branch)`, select the latest complete attempt for each stage, recompute its fingerprint and outputs, and return exactly one of:

```js
{ ok: true, nextStage: "dock", reused: [{ stage: "prepare", attempt: 1 }, { stage: "score", attempt: 1 }] }
{ ok: true, nextStage: null, reused: [...] }
{ ok: false, code: "checkpoint_mismatch", stage: "score", resource: "outputs/scores.csv" }
{ ok: false, code: "checkpoint_invalid", stage: "score", resource: "checkpoint.json" }
```

- [ ] **Step 5: Integrate manifest refresh**

Expose `atomicWriteJson` from `store.mjs` or pass `artifactStore.refreshManifest` into the stage store. Every checkpoint transition must finish its rename before appending an event and refreshing `MANIFEST.sha256`.

- [ ] **Step 6: Run stage contract checks**

Run: `npm run stage-contract-check`  
Expected: stable fingerprint, immutable attempt, tamper detection and contiguous verification cases all pass.

## Task 4: Replace the Monolithic Worker with Stage Orchestration

**Files:**
- Modify: `apps/open-molecule-lab/server/worker.mjs`
- Modify: `apps/open-molecule-lab/server/store.mjs`
- Modify: `apps/open-molecule-lab/server/server.mjs`
- Modify: `apps/open-molecule-lab/server/python-bridge.py`
- Modify: `apps/open-molecule-lab/server/worker-lifecycle-check.mjs`

- [ ] **Step 1: Extend worker lifecycle RED fixtures**

Make the fake `scoring.py` distinguish base and dock calls. Assert a cascade run invokes score without `--base-scores`, then dock with the exact prior score path; a library run never invokes dock. Assert cancellation writes the active StageAttempt terminal state before `shutdown()` resolves.

- [ ] **Step 2: Run and verify RED**

Run: `npm run worker-lifecycle-check`  
Expected: current worker has no StageAttempt orchestration and fails the invocation assertions.

- [ ] **Step 3: Create the prepare attempt during execution-run creation**

After strict preflight succeeds and the attached RunSpec is persisted, compute code identity and prepare fingerprint. Create and complete prepare attempt 1 with logical outputs:

```js
{
  "run-spec.json": "run-spec.json",
  "preflight.json": "preflight.json",
  "molecules.csv": "inputs/molecules.csv",
}
```

Blocked preflight runs retain existing blocked behavior and do not create a fake complete prepare attempt.

- [ ] **Step 4: Build exact stage commands**

Score:

```js
[
  scoringEntry,
  "--input", inputPath,
  "--output", scoreOutput,
  "--target", spec.target.id,
  "--mode", "library",
  "--strict-backends",
  "--asset-root", assetRoot,
]
```

Dock:

```js
[
  scoringEntry,
  "--input", inputPath,
  "--base-scores", verifiedScoreOutput,
  "--output", dockOutput,
  "--target", spec.target.id,
  "--mode", "cascade",
  "--strict-backends",
  "--asset-root", assetRoot,
  "--cascade-top-n", String(spec.policy?.dockTopN || 300),
]
```

No shell string is allowed. Attempt-local `command.json`, `stdout.log` and `stderr.log` use logical placeholders and path redaction.

- [ ] **Step 5: Execute score and dock sequentially**

Refactor `start(run, { startStage } = {})` into an async loop over remaining stages. For score/dock, create attempt -> running -> spawn detached group -> validate exit/output -> complete. Store the active `{child,stage,attempt,cancelRequested,donePromise}` in the process map.

Extend `summarize-results` so its complete scan returns:

```python
structure_ok = (
    int(frame["structure_docking_status"].eq("ok").sum())
    if "structure_docking_status" in frame.columns else 0
)
summary["structureDockingOk"] = structure_ok
```

Score completion must invoke `python-bridge.py summarize-results --expected-input`. Dock completion additionally requires `rankingScoreField === "final_score_dock"` and `structureDockingOk > 0`; do not infer success from a bounded result page.

- [ ] **Step 6: Implement report as a real stage**

Create report attempt, select verified score output for library or dock output for cascade, copy it atomically to `results/scores.csv`, run the summary bridge, write `results/summary.json`, complete report, then mark run complete. If copy, summary, identity, redaction or manifest verification fails, report attempt and run fail with `evidence_incomplete`.

- [ ] **Step 7: Preserve process-group lifecycle semantics**

Cancellation and shutdown terminate only the active attempt group, then wait for attempt checkpoint and run state persistence. Restart recovery must mark the active attempt `failed/worker_interrupted` before marking the run failed. Existing unrelated-PID rejection remains unchanged.

- [ ] **Step 8: Run worker and source-only contract checks**

Run:

```bash
npm run worker-lifecycle-check
npm run stage-contract-check
npm run contract-check
```

Expected: all pass; source-only execution remains fail-closed without assets.

## Task 5: Add Explicit Resume and Stage APIs

**Files:**
- Modify: `apps/open-molecule-lab/server/server.mjs`
- Modify: `apps/open-molecule-lab/server/store.mjs`
- Modify: `apps/open-molecule-lab/server/stage-store.mjs`
- Modify: `apps/open-molecule-lab/server/contract-check.mjs`

- [ ] **Step 1: Add failing HTTP cases**

Add helpers for:

```js
POST /api/runs/:id/resume
GET /api/runs/:id/stages
```

Assert blocked and complete runs return 409 on resume, unknown runs return 404, and stage GET returns a versioned empty/blocked summary without host paths. Unit fixtures must cover failed/cancelled resume, mismatch blocking and attempt increment.

- [ ] **Step 2: Run and verify RED**

Run: `npm run contract-check`  
Expected: resume/stages endpoints return 404.

- [ ] **Step 3: Implement fresh resume preflight**

Resume reloads immutable RunSpec and MoleculeSet metadata, invokes `runPreflight` with the persisted route, and requires all checks to pass. It must not re-resolve or switch the route branch.

- [ ] **Step 4: Implement resume verification and start**

For a failed/cancelled/interrupted run:

```js
const verification = await stageStore.verifyContiguous(run, await buildStageInputs(run));
if (!verification.ok) {
  const blocked = await artifactStore.updateRun(runId, {
    status: "blocked",
    finishedAt: now,
    error: { code: verification.code, message: "Verified resume checkpoint mismatch", stage: verification.stage },
  }, { event: "resume_blocked", status: "blocked", stage: verification.stage });
  return jsonResponse(response, 409, { ok: true, ...blocked });
}
if (!verification.nextStage) return jsonResponse(response, 409, errorPayload("run_already_complete", "No stage remains"));
await artifactStore.updateRun(runId, {
  status: "queued", finishedAt: null, error: null,
  resumeCount: Number(run.resumeCount || 0) + 1,
}, { event: "resume_started", status: "queued", stage: verification.nextStage });
await workerManager.start(await artifactStore.getRun(runId), { startStage: verification.nextStage });
```

- [ ] **Step 5: Return persisted stage summaries**

`GET /api/runs/:id/stages` returns:

```json
{
  "ok": true,
  "schemaVersion": "open-molecule-lab.stage-summary.v0.1",
  "runId": "run_...",
  "resumable": true,
  "stages": [
    { "stage": "prepare", "status": "complete", "attempts": 1 },
    { "stage": "score", "status": "complete", "attempts": 1 },
    { "stage": "dock", "status": "failed", "attempts": 1 },
    { "stage": "report", "status": "waiting", "attempts": 0 }
  ]
}
```

Do not return absolute attempt roots or physical asset paths.

- [ ] **Step 6: Run API and stage checks**

Run: `npm run contract-check && npm run stage-contract-check`  
Expected: all endpoint, mismatch and immutable-history cases pass.

## Task 6: Render Real Stage History and Resume Controls

**Files:**
- Modify: `apps/open-molecule-lab/src/types.ts`
- Modify: `apps/open-molecule-lab/src/App.tsx`
- Modify: `apps/open-molecule-lab/src/styles.css`

- [ ] **Step 1: Add exact TypeScript contracts**

```ts
export type AttemptStatus = "queued" | "running" | "complete" | "failed" | "cancelled" | "blocked";

export interface StageSummaryRow {
  stage: "prepare" | "score" | "dock" | "report";
  status: AttemptStatus | "waiting" | "skipped";
  attempts: number;
  startedAt?: string | null;
  finishedAt?: string | null;
  errorCode?: string | null;
}

export interface StageSummary {
  ok: true;
  schemaVersion: "open-molecule-lab.stage-summary.v0.1";
  runId: string;
  resumable: boolean;
  stages: StageSummaryRow[];
}
```

- [ ] **Step 2: Fetch stages with the run poll**

When a run exists, fetch `/api/runs/:id/stages` after launch and on every poll. Use the existing request-generation guard so a previous run cannot overwrite the current stage timeline.

- [ ] **Step 3: Replace planned-stage rendering after launch**

Before launch, retain plan stages. After launch, render persisted StageSummary rows with fixed dimensions, attempt count and terminal error code. Do not estimate percentages. Map `planned/waiting` to informational/idle, not green.

- [ ] **Step 4: Add explicit resume command**

Use `RotateCcw` from lucide-react. Show the icon+text command only when `stageSummary.resumable` is true and the run is not active. POST resume, replace run state with the response and restart polling. Keep completed/blocked-mismatch runs non-resumable.

- [ ] **Step 5: Add responsive styles**

Use a single unframed timeline list; each row has stable grid tracks for index, label, state and attempt metadata. At `max-width: 560px`, wrap metadata below the label. Do not add nested cards or viewport-scaled fonts.

- [ ] **Step 6: Build and inspect**

Run: `npm run build`  
Expected: TypeScript and Vite pass. Start the server and inspect desktop/mobile widths for overlap, stable stage rows and readable resume state.

## Task 7: Prove Real Interruption, Resume and Parity

**Files:**
- Modify: `apps/open-molecule-lab/server/real-run-check.mjs`
- Modify: `apps/open-molecule-lab/server/worker-lifecycle-check.mjs`

- [ ] **Step 1: Refactor the checker to restart the same server state**

Replace the single const child with `startServer()`/`stopServer(signal)` helpers that reuse the same `dataRoot`, `runsRoot` and port. `stopServer("SIGTERM")` must wait for exit.

- [ ] **Step 2: Assert uninterrupted stage chains**

For library require `prepare/score/report == complete` and dock skipped. For cascade require `prepare/score/dock/report == complete`, one attempt each, real structure docking rows and `rankingScoreField == final_score_dock`.

- [ ] **Step 3: Interrupt after score checkpoint**

Launch a cascade run, poll `stages/score/attempt-0001/checkpoint.json` until complete, record the score file SHA-256, then stop the server while dock/report is incomplete. Restart the server with identical environment and assert the run becomes failed or cancelled with score attempt 1 still complete.

- [ ] **Step 4: Resume and verify immutable reuse**

POST resume, wait complete, and assert:

```js
assert(await sha256(scorePath) === scoreShaBefore, "resume overwrote score output");
assert(stageSummary.stages.find((row) => row.stage === "score").attempts === 1, "score reran");
assert(stageSummary.stages.find((row) => row.stage === "dock").attempts >= 1, "dock did not resume");
```

Compare scientific columns from resumed and uninterrupted cascade with absolute tolerance `1e-4`; exclude timestamps, PIDs, paths and attempt metadata.

- [ ] **Step 5: Prove tamper fail-closed**

Create a second interrupted run, modify one base `final_score`, refresh only the outer run manifest to isolate checkpoint verification, restart and POST resume. Expect HTTP 409, run status blocked/checkpoint_mismatch, no new dock attempt and no smina process.

- [ ] **Step 6: Run the asset-backed checker**

Run:

```bash
OPEN_MOLECULE_PYTHON=<workspace>/.venv_four_level_cli/bin/python \
OPEN_MOLECULE_ASSET_ROOT=<workspace>/open_source_release/four-level-molecule-cli-offline-assets \
SMINA_BIN=<workspace>/ai_drug_eval_tools/micromamba-root/envs/smina-local/bin/smina \
OBABEL_BIN=<workspace>/ai_drug_eval_tools/micromamba-root/envs/smina-local/bin/obabel \
npm run real-run-check
```

Expected: uninterrupted library/cascade, interrupted-resumed cascade parity and tamper-block cases all pass.

## Task 8: Final Verification, Documentation and Publication

**Files:**
- Modify: `apps/open-molecule-lab/README.md`
- Modify: `apps/open-molecule-lab/DESIGN.md`
- Modify: `README.md`
- Modify: `VALIDATION.md`
- Modify: `FULL_SOURCE_INVENTORY.md`
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Document exact resume boundaries**

Document StageAttempt paths, fingerprints, `GET /stages`, `POST /resume`, immutable attempt history, mismatch blocking and the explicit non-goal of mid-process memory resume.

- [ ] **Step 2: Add source-only CI checks**

Run these in the Open Molecule Lab CI job after build:

```yaml
- run: npm run contract-check
- run: npm run stage-contract-check
- run: npm run worker-lifecycle-check
```

Do not add the asset-backed real-run checker to public CI.

- [ ] **Step 3: Run the complete fresh verification matrix**

Run:

```bash
<workspace>/.venv_four_level_cli/bin/python -m pytest -q
<workspace>/.venv_four_level_cli/bin/python \
  -m scientific_validation.four_level_cli_1kx10k.verify_snapshot \
  --snapshot-dir validation/frozen_run
<workspace>/.venv_four_level_cli/bin/python \
  scripts/verify_assets.py \
  --asset-root <workspace>/open_source_release/four-level-molecule-cli-offline-assets
npm --prefix apps/open-molecule-lab run build
npm --prefix apps/open-molecule-lab run contract-check
npm --prefix apps/open-molecule-lab run stage-contract-check
npm --prefix apps/open-molecule-lab run worker-lifecycle-check
```

Then run the Task 7 asset-backed command. Every command must exit 0 before release claims are updated.

- [ ] **Step 4: Restart and verify the local product**

Restart `npm run serve` with explicit Python/assets/smina/obabel variables. Verify `/api/health` is `local_execution`, all model rows are available, the built asset hash is current, and no scoring/smina process remains after a cancellation test.

- [ ] **Step 5: Rebuild the deterministic source archive**

Use `scripts/build_source_release.py`, run `unzip -t`, verify excluded paths/suffixes, record source file count, expanded bytes, compressed bytes and SHA-256, then update the release MANIFEST/SHA256SUMS.

- [ ] **Step 6: Replace GitHub source and Release assets**

Sync from the clean ZIP extraction so obsolete files are deleted. Compare every remote Git blob SHA to the archive, require `ai_mol_loop/ == 0`, upload Release assets with clobber, download the ZIP and require byte-identical comparison.

The current token lacks `workflow`; report `.github/workflows/ci.yml` separately and do not claim exact remote parity until `gh auth refresh -h github.com -s workflow` succeeds.

## Plan Self-Review

- The plan covers every requirement in the verified-resume design: real score/dock separation, immutable attempts, deterministic fingerprints, explicit resume, restart recovery, UI state, tamper blocking and real parity.
- It does not substitute batch CLI checkpoints for interactive uploaded MoleculeSets.
- It preserves default `scoring.py` behavior and does not change scientific formulas.
- It does not claim run comparison, evidence export or 10,000-candidate scheduling; those remain subsequent product increments under the approved high-level specification.
- All paths, function names, endpoint names, commands and terminal states are fixed; there are no implementation placeholders.
