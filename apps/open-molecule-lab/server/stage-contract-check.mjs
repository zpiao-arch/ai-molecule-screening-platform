import { createHash } from "node:crypto";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";

import {
  canonicalJson,
  computeCodeIdentity,
  computeStageFingerprint,
  stageSequence,
} from "./stage-contract.mjs";
import { createStageStore } from "./stage-store.mjs";


const __dirname = path.dirname(fileURLToPath(import.meta.url));
const sourceRoot = path.resolve(__dirname, "..", "..", "..");

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

async function assertRejects(operation, pattern, message) {
  try {
    await operation();
  } catch (error) {
    assert(pattern.test(String(error?.message || error)), message);
    return;
  }
  throw new Error(message);
}

function digest(value) {
  return createHash("sha256").update(value).digest("hex");
}

function stageInput(stage, overrides = {}) {
  return {
    schemaVersion: "open-molecule-lab.stage-input.v0.1",
    stage,
    runSpecSha256: "1".repeat(64),
    moleculeSetSha256: "2".repeat(64),
    assetManifestSha256: "3".repeat(64),
    codeIdentity: "4".repeat(64),
    routeBranch: "cascade",
    stagePolicy: {},
    upstreamCheckpointSha256: null,
    upstreamOutputSha256: null,
    ...overrides,
  };
}

async function writeAttemptOutput(attempt, name, content) {
  const relativePath = path.posix.join("outputs", name);
  await fs.writeFile(path.join(attempt.root, relativePath), content, "utf8");
  return relativePath;
}

async function completeFixtureStage(stages, runId, stage, input, outputName = `${stage}.json`) {
  const command = { executable: "OPEN_MOLECULE_PYTHON", args: [stage] };
  const attempt = await stages.createAttempt(
    runId,
    stage,
    computeStageFingerprint(input),
    command,
    input,
  );
  await stages.transition(runId, stage, attempt.attempt, "running");
  const outputPath = await writeAttemptOutput(attempt, outputName, `${stage}\n`);
  return stages.complete(runId, stage, attempt.attempt, { [outputName]: outputPath });
}

async function checkStageStore() {
  const temporaryRoot = await fs.mkdtemp(path.join(os.tmpdir(), "open-molecule-stage-store-"));
  const runsRoot = path.join(temporaryRoot, "runs");
  const events = [];
  const refreshes = [];
  const stages = createStageStore({
    runsRoot,
    appendEvent: async (runId, event) => events.push({ runId, ...event }),
    refreshManifest: async (runRoot) => refreshes.push(runRoot),
  });
  try {
    const runId = "run_stage_store";
    await fs.mkdir(path.join(runsRoot, runId), { recursive: true });
    const scoreInput = stageInput("score");
    const fingerprint = computeStageFingerprint(scoreInput);
    const command = { executable: "OPEN_MOLECULE_PYTHON", args: ["scoring/scoring.py"] };

    const first = await stages.createAttempt(runId, "score", fingerprint, command, scoreInput);
    assert(first.attempt === 1 && first.status === "queued", "first attempt was not queued");
    await stages.transition(runId, "score", 1, "running");
    const outputPath = await writeAttemptOutput(first, "scores.csv", "id,smiles\na,CCO\n");
    const complete = await stages.complete(runId, "score", 1, { "scores.csv": outputPath });
    assert(complete.status === "complete", "attempt did not complete");
    assert(
      complete.outputs["scores.csv"].sha256 === digest("id,smiles\na,CCO\n"),
      "stage store accepted a caller digest instead of hashing the output",
    );
    assert(
      complete.outputs["scores.csv"].path === "outputs/scores.csv",
      "stage output path was not persisted logically",
    );

    const firstCheckpointPath = path.join(first.root, "checkpoint.json");
    const firstBytes = await fs.readFile(firstCheckpointPath);
    const second = await stages.createAttempt(runId, "score", fingerprint, command, scoreInput);
    assert(second.attempt === 2 && second.status === "queued", "attempt number was not incremented");
    assert(
      Buffer.compare(firstBytes, await fs.readFile(firstCheckpointPath)) === 0,
      "creating attempt 2 rewrote attempt 1",
    );
    await assertRejects(
      () => stages.transition(runId, "score", 1, "running"),
      /terminal attempt/i,
      "terminal attempt was mutable",
    );
    assert(events.map((event) => event.event).slice(0, 4).join(",")
      === "stage_attempt_queued,stage_attempt_started,stage_attempt_complete,stage_attempt_queued",
    "stage event order was not persisted");
    assert(refreshes.length === 4, "manifest was not refreshed after each checkpoint transition");
    assert(!firstBytes.toString("utf8").includes(temporaryRoot), "checkpoint exposed an absolute path");

    await fs.writeFile(path.join(first.root, outputPath), "id,smiles\na,CCC\n", "utf8");
    const tampered = await stages.verifyAttempt(runId, "score", 1, fingerprint);
    assert(!tampered.ok && tampered.code === "output_mismatch", "tampered output was accepted");

    const cascadeRunId = "run_cascade_resume";
    await fs.mkdir(path.join(runsRoot, cascadeRunId), { recursive: true });
    const cascadeInputs = {
      prepare: stageInput("prepare"),
      score: stageInput("score"),
      dock: stageInput("dock"),
      report: stageInput("report"),
    };
    await completeFixtureStage(stages, cascadeRunId, "prepare", cascadeInputs.prepare);
    const cascadeScore = await completeFixtureStage(
      stages,
      cascadeRunId,
      "score",
      cascadeInputs.score,
      "scores.csv",
    );
    const cascadeRun = { runId: cascadeRunId, route: { branch: "cascade" } };
    const resumable = await stages.verifyContiguous(cascadeRun, cascadeInputs);
    assert(resumable.ok && resumable.nextStage === "dock", "cascade did not resume at dock");
    assert(
      JSON.stringify(resumable.reused) === JSON.stringify([
        { stage: "prepare", attempt: 1 },
        { stage: "score", attempt: 1 },
      ]),
      "contiguous reuse history was incorrect",
    );
    const cascadeScorePath = path.join(
      runsRoot,
      cascadeRunId,
      "stages",
      "score",
      "attempt-0001",
      cascadeScore.outputs["scores.csv"].path,
    );
    await fs.appendFile(cascadeScorePath, "tampered,CCC\n", "utf8");
    const blocked = await stages.verifyContiguous(cascadeRun, cascadeInputs);
    assert(
      !blocked.ok
        && blocked.code === "checkpoint_mismatch"
        && blocked.stage === "score"
        && blocked.resource === "outputs/scores.csv",
      "contiguous verification did not fail closed on output tampering",
    );

    const libraryRunId = "run_library_resume";
    await fs.mkdir(path.join(runsRoot, libraryRunId), { recursive: true });
    const libraryInputs = {
      prepare: stageInput("prepare", { routeBranch: "library" }),
      score: stageInput("score", { routeBranch: "library" }),
      report: stageInput("report", { routeBranch: "library" }),
    };
    for (const stage of stageSequence("library")) {
      await completeFixtureStage(stages, libraryRunId, stage, libraryInputs[stage]);
    }
    const libraryRun = { runId: libraryRunId, route: { branch: "library" } };
    const libraryComplete = await stages.verifyContiguous(libraryRun, libraryInputs);
    assert(libraryComplete.ok && libraryComplete.nextStage === null, "complete library run remained resumable");
    const summary = await stages.summarize(libraryRun);
    assert(
      summary.stages.map((row) => `${row.stage}:${row.status}`).join(",")
        === "prepare:complete,score:complete,dock:skipped,report:complete",
      "library summary did not mark dock as skipped",
    );
  } finally {
    await fs.rm(temporaryRoot, { recursive: true, force: true });
  }
}

assert(
  canonicalJson({ b: 2, nested: { z: true, a: null }, a: 1 })
    === '{"a":1,"b":2,"nested":{"a":null,"z":true}}',
  "canonical JSON did not sort object keys recursively",
);
assert(
  JSON.stringify(stageSequence("library")) === JSON.stringify(["prepare", "score", "report"]),
  "library stage order changed",
);
assert(
  JSON.stringify(stageSequence("cascade")) === JSON.stringify(["prepare", "score", "dock", "report"]),
  "cascade stage order changed",
);

const fingerprintInput = {
  schemaVersion: "open-molecule-lab.stage-input.v0.1",
  stage: "score",
  runSpecSha256: "a".repeat(64),
  moleculeSetSha256: "b".repeat(64),
  assetManifestSha256: "c".repeat(64),
  codeIdentity: "d".repeat(64),
  routeBranch: "cascade",
  stagePolicy: { role: "base_four_level_score" },
  upstreamCheckpointSha256: "e".repeat(64),
  upstreamOutputSha256: "f".repeat(64),
};
const first = computeStageFingerprint(fingerprintInput);
const repeated = computeStageFingerprint({ ...fingerprintInput, stagePolicy: { role: "base_four_level_score" } });
const changed = computeStageFingerprint({ ...fingerprintInput, moleculeSetSha256: "9".repeat(64) });
assert(first.length === 64 && first === repeated, "stage fingerprint is not deterministic");
assert(first !== changed, "molecule-set hash did not affect stage fingerprint");

const identity = await computeCodeIdentity(sourceRoot);
assert(identity.schemaVersion === "open-molecule-lab.code-identity.v0.1", "code identity schema missing");
assert(identity.sha256?.length === 64, "code identity digest missing");
assert(identity.files.some((row) => row.path === "scoring/scoring.py"), "scientific entry missing from identity");
assert(
  identity.files.some((row) => row.path === "apps/open-molecule-lab/server/worker.mjs"),
  "worker entry missing from identity",
);
assert(identity.files.every((row) => !path.isAbsolute(row.path)), "code identity exposed an absolute path");
assert(
  JSON.stringify([...identity.files].sort((left, right) => left.path.localeCompare(right.path)))
    === JSON.stringify(identity.files),
  "code identity files are not sorted",
);

await checkStageStore();

process.stdout.write(`${JSON.stringify({
  ok: true,
  fingerprint: first,
  codeIdentity: identity.sha256,
  files: identity.files.length,
}, null, 2)}\n`);
