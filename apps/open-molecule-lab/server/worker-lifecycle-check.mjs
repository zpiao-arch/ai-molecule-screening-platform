import { spawn } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";

import { computeStageFingerprint, stageSequence } from "./stage-contract.mjs";
import { createStageStore } from "./stage-store.mjs";
import { createWorkerManager, terminateOwnedPid, terminateProcessGroup } from "./worker.mjs";


function assert(condition, message) {
  if (!condition) throw new Error(message);
}

function processExists(pid) {
  try {
    process.kill(pid, 0);
    return true;
  } catch (error) {
    return error?.code === "EPERM";
  }
}

async function waitForExit(pid) {
  for (let attempt = 0; attempt < 100; attempt += 1) {
    if (!processExists(pid)) return;
    await new Promise((resolve) => setTimeout(resolve, 25));
  }
  throw new Error(`process ${pid} survived group termination`);
}

async function waitForRun(store, statuses) {
  const expected = new Set(Array.isArray(statuses) ? statuses : [statuses]);
  for (let attempt = 0; attempt < 200; attempt += 1) {
    if (expected.has(store.run.status)) return store.run;
    await new Promise((resolve) => setTimeout(resolve, 25));
  }
  throw new Error(`run did not reach ${[...expected].join("/")}; current=${store.run.status}`);
}

async function waitForStage(stageStore, run, stage, status) {
  for (let attempt = 0; attempt < 200; attempt += 1) {
    const summary = await stageStore.summarize(run);
    if (summary.stages.find((row) => row.stage === stage)?.status === status) return summary;
    await new Promise((resolve) => setTimeout(resolve, 25));
  }
  throw new Error(`${stage} did not reach ${status}`);
}

function argumentValue(args, flag) {
  const index = args.indexOf(flag);
  return index < 0 ? null : args[index + 1];
}

function fixtureStageInput(stage, branch) {
  return {
    schemaVersion: "open-molecule-lab.stage-input.v0.1",
    stage,
    runSpecSha256: "1".repeat(64),
    moleculeSetSha256: "2".repeat(64),
    assetManifestSha256: "3".repeat(64),
    codeIdentity: "4".repeat(64),
    routeBranch: branch,
    stagePolicy: { role: stage },
    upstreamCheckpointSha256: null,
    upstreamOutputSha256: null,
  };
}

const fakeScoring = String.raw`
const fs = require("node:fs");
const path = require("node:path");
const { spawn } = require("node:child_process");
const args = process.argv.slice(2);
const value = (flag) => args[args.indexOf(flag) + 1];
fs.appendFileSync(path.join(__dirname, "..", "invocations.jsonl"), JSON.stringify(args) + "\n");
if (process.env.OPEN_MOLECULE_FIXTURE_HANG === "1") {
  spawn(process.execPath, ["-e", "setInterval(() => {}, 1000)"], { stdio: "ignore" });
  setInterval(() => {}, 1000);
} else {
  const input = value("--input");
  const output = value("--output");
  const base = args.includes("--base-scores") ? value("--base-scores") : "";
  const cascade = value("--mode") === "cascade";
  const inputRows = fs.readFileSync(input, "utf8").trim().split(/\r?\n/).slice(1);
  const baseHeader = "id,smiles,layer1_status,layer2_status,layer3_status,layer4_status,docking_normalized,final_score";
  const baseRows = inputRows.map((row, index) => row + ",ok,ok,ok,ok,0.5," + (0.8 - index * 0.1));
  let lines = [baseHeader, ...baseRows];
  if (cascade) {
    if (base) lines = fs.readFileSync(base, "utf8").trim().split(/\r?\n/);
    lines[0] += ",docking_affinity_kcal_mol,heavy_atoms,ligand_efficiency,dock_rerank_rank,structure_docking_status,dock_le_norm,final_score_dock";
    lines = [lines[0], ...lines.slice(1).map((row, index) => row + ",-7.5,3,-2.5," + (index + 1) + ",ok,0.9," + (0.95 - index * 0.1))];
  }
  fs.mkdirSync(path.dirname(output), { recursive: true });
  fs.writeFileSync(output, lines.join("\n") + "\n", "utf8");
}
`;

const fakeBridge = String.raw`
const fs = require("node:fs");
const args = process.argv.slice(2);
if (args[0] !== "summarize-results") process.exit(2);
const input = args[args.indexOf("--input") + 1];
const lines = fs.readFileSync(input, "utf8").trim().split(/\r?\n/);
const columns = lines[0].split(",");
const docked = columns.includes("final_score_dock");
process.stdout.write(JSON.stringify({
  ok: true,
  nRows: lines.length - 1,
  nRanked: lines.length - 1,
  nFailed: 0,
  columns,
  rankingScoreField: docked ? "final_score_dock" : "final_score",
  structureDockingOk: docked ? lines.length - 1 : 0,
}) + "\n");
`;

async function createRunFixture(branch, { cancelledWriteDelay = 0 } = {}) {
  const fixtureRoot = await fs.mkdtemp(path.join(os.tmpdir(), `open-molecule-worker-${branch}-`));
  const sourceRoot = path.join(fixtureRoot, "source");
  const runsRoot = path.join(fixtureRoot, "runs");
  const runId = `run_20260723000000_${branch}`;
  const runRoot = path.join(runsRoot, runId);
  await fs.mkdir(path.join(sourceRoot, "scoring"), { recursive: true });
  await fs.mkdir(path.join(runRoot, "inputs"), { recursive: true });
  await fs.mkdir(path.join(runRoot, "results"), { recursive: true });
  await fs.mkdir(path.join(runRoot, "logs"), { recursive: true });
  await fs.writeFile(path.join(sourceRoot, "scoring", "scoring.py"), fakeScoring, "utf8");
  await fs.writeFile(path.join(sourceRoot, "bridge.py"), fakeBridge, "utf8");
  await fs.writeFile(path.join(runRoot, "inputs", "molecules.csv"), "id,smiles\nmol-1,CCO\nmol-2,CCN\n", "utf8");
  const spec = {
    target: { id: "CHEMBL2051" },
    moleculeSet: { inputSha256: "2".repeat(64) },
    execution: { strictBackends: true },
    policy: { dockTopN: 2 },
  };
  await fs.writeFile(path.join(runRoot, "run-spec.json"), `${JSON.stringify(spec, null, 2)}\n`, "utf8");
  const run = {
    runId,
    status: "queued",
    route: { branch },
    specSha256: "1".repeat(64),
    startedAt: null,
    finishedAt: null,
    error: null,
  };
  const events = [];
  const store = {
    runsRoot,
    run,
    async refreshManifest() {},
    async verifyManifest() { return { ok: true }; },
    async appendEvent(_runId, event) { events.push(event); },
    async getRun() { return this.run; },
    async updateRun(_runId, patch, event) {
      if (patch.status === "cancelled" && cancelledWriteDelay) {
        await new Promise((resolve) => setTimeout(resolve, cancelledWriteDelay));
      }
      this.run = { ...this.run, ...patch };
      if (event) events.push(event);
      return this.run;
    },
  };
  const stageStore = createStageStore({
    runsRoot,
    appendEvent: (eventRunId, event) => store.appendEvent(eventRunId, event),
    refreshManifest: (eventRunRoot) => store.refreshManifest(eventRunRoot),
  });
  const stageInputs = Object.fromEntries(
    stageSequence(branch).map((stage) => [stage, fixtureStageInput(stage, branch)]),
  );
  const prepare = await stageStore.createAttempt(
    runId,
    "prepare",
    computeStageFingerprint(stageInputs.prepare),
    { executable: "OPEN_MOLECULE_SERVER", args: ["prepare"] },
    stageInputs.prepare,
  );
  await stageStore.transition(runId, "prepare", prepare.attempt, "running");
  await fs.writeFile(path.join(prepare.root, "outputs", "run-spec.json"), `${JSON.stringify(spec)}\n`, "utf8");
  await stageStore.complete(runId, "prepare", prepare.attempt, {
    "run-spec.json": "outputs/run-spec.json",
  });

  const manager = createWorkerManager({
    store,
    stageStore,
    buildStageInputs: async () => stageInputs,
    sourceRoot,
    pythonPath: process.execPath,
    bridgePath: path.join(sourceRoot, "bridge.py"),
    assetRoot: fixtureRoot,
    sminaBin: "fixture-smina",
    obabelBin: "fixture-obabel",
  });
  return {
    fixtureRoot,
    run,
    store,
    stageStore,
    manager,
    invocationPath: path.join(sourceRoot, "invocations.jsonl"),
    async invocations() {
      const text = await fs.readFile(this.invocationPath, "utf8").catch(() => "");
      return text.trim().split("\n").filter(Boolean).map((line) => JSON.parse(line));
    },
    async cleanup() {
      await manager.shutdown();
      await fs.rm(fixtureRoot, { recursive: true, force: true });
    },
  };
}

const child = spawn(
  process.execPath,
  [
    "-e",
    [
      "const { spawn } = require('node:child_process');",
      "const grandchild = spawn(process.execPath, ['-e', 'setInterval(() => {}, 1000)'], { stdio: 'ignore' });",
      "process.stdout.write(String(grandchild.pid) + '\\n');",
      "setInterval(() => {}, 1000);",
    ].join(""),
  ],
  { detached: true, stdio: ["ignore", "pipe", "inherit"] },
);

const grandchildPid = Number(await new Promise((resolve, reject) => {
  child.stdout.once("data", (chunk) => resolve(String(chunk).trim()));
  child.once("error", reject);
}));

try {
  assert(Number.isInteger(child.pid), "worker PID missing");
  assert(Number.isInteger(grandchildPid), "grandchild PID missing");
  assert(terminateProcessGroup(child, "SIGTERM"), "group termination was not attempted");
  await Promise.all([waitForExit(child.pid), waitForExit(grandchildPid)]);

  const unrelated = spawn(process.execPath, ["-e", "setInterval(() => {}, 1000)"], {
    detached: true,
    stdio: "ignore",
  });
  try {
    assert(
      terminateOwnedPid(unrelated.pid, { sourceRoot: "/tmp/not-this-worker", pythonPath: "/tmp/not-python" }) === false,
      "orphan recovery accepted an unrelated reused PID",
    );
    assert(processExists(unrelated.pid), "orphan recovery terminated an unrelated process");
  } finally {
    terminateProcessGroup(unrelated, "SIGKILL");
  }

  const cascade = await createRunFixture("cascade");
  try {
    await cascade.manager.start(cascade.run);
    await waitForRun(cascade.store, "complete");
    const invocations = await cascade.invocations();
    assert(invocations.length === 2, `cascade invoked scoring ${invocations.length} times instead of 2`);
    assert(argumentValue(invocations[0], "--mode") === "library", "cascade score stage did not use library mode");
    assert(!invocations[0].includes("--base-scores"), "score stage unexpectedly consumed base scores");
    assert(argumentValue(invocations[1], "--mode") === "cascade", "dock stage did not use cascade mode");
    const expectedBase = path.join(
      cascade.store.runsRoot,
      cascade.run.runId,
      "stages",
      "score",
      "attempt-0001",
      "outputs",
      "scores.csv",
    );
    assert(argumentValue(invocations[1], "--base-scores") === expectedBase, "dock did not consume exact score output");
    assert(argumentValue(invocations[1], "--cascade-top-n") === "2", "dock policy was not preserved");
    const summary = await cascade.stageStore.summarize(cascade.store.run);
    assert(
      summary.stages.map((row) => `${row.stage}:${row.status}`).join(",")
        === "prepare:complete,score:complete,dock:complete,report:complete",
      "cascade stage chain was incomplete",
    );
  } finally {
    await cascade.cleanup();
  }

  const library = await createRunFixture("library");
  try {
    await library.manager.start(library.run);
    await waitForRun(library.store, "complete");
    const invocations = await library.invocations();
    assert(invocations.length === 1, `library invoked scoring ${invocations.length} times instead of 1`);
    assert(argumentValue(invocations[0], "--mode") === "library", "library score mode changed");
    assert(!invocations[0].includes("--base-scores"), "library attempted docking base-score mode");
    const summary = await library.stageStore.summarize(library.store.run);
    assert(summary.stages.find((row) => row.stage === "dock")?.status === "skipped", "library dock was not skipped");
  } finally {
    await library.cleanup();
  }

  const cancellation = await createRunFixture("library", { cancelledWriteDelay: 800 });
  const previousHang = process.env.OPEN_MOLECULE_FIXTURE_HANG;
  process.env.OPEN_MOLECULE_FIXTURE_HANG = "1";
  try {
    await cancellation.manager.start(cancellation.run);
    await waitForStage(cancellation.stageStore, cancellation.store.run, "score", "running");
    await cancellation.manager.shutdown();
    assert(cancellation.store.run.status === "cancelled", "shutdown returned before run cancellation persisted");
    const summary = await cancellation.stageStore.summarize(cancellation.store.run);
    assert(
      summary.stages.find((row) => row.stage === "score")?.status === "cancelled",
      "shutdown returned before active StageAttempt cancellation persisted",
    );
  } finally {
    if (previousHang === undefined) delete process.env.OPEN_MOLECULE_FIXTURE_HANG;
    else process.env.OPEN_MOLECULE_FIXTURE_HANG = previousHang;
    await cancellation.cleanup();
  }

  process.stdout.write(`${JSON.stringify({ ok: true, pid: child.pid, grandchildPid }, null, 2)}\n`);
} finally {
  if (processExists(child.pid)) process.kill(child.pid, "SIGKILL");
  if (processExists(grandchildPid)) process.kill(grandchildPid, "SIGKILL");
}
