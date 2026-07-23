import { createHash } from "node:crypto";
import { spawn, spawnSync } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { runPreflight } from "./preflight.mjs";


const __dirname = path.dirname(fileURLToPath(import.meta.url));
const appRoot = path.resolve(__dirname, "..");
const sourceRoot = path.resolve(appRoot, "..", "..");
const pythonPath = process.env.OPEN_MOLECULE_PYTHON;
const assetRoot = process.env.OPEN_MOLECULE_ASSET_ROOT;
if (!pythonPath || !assetRoot) {
  throw new Error("OPEN_MOLECULE_PYTHON and OPEN_MOLECULE_ASSET_ROOT are required");
}

const dataRoot = await fs.mkdtemp(path.join(os.tmpdir(), "open-molecule-lab-real-data-"));
const runsRoot = await fs.mkdtemp(path.join(os.tmpdir(), "open-molecule-lab-real-runs-"));
const port = 36000 + (process.pid % 10000);
const baseUrl = `http://127.0.0.1:${port}`;
let serverChild = null;
let childOutput = "";

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

function wait(milliseconds) {
  return new Promise((resolve) => setTimeout(resolve, milliseconds));
}

async function requestJson(url, options) {
  const response = await fetch(url, options);
  const payload = await response.json();
  return { response, payload };
}

async function post(pathname, body = {}) {
  return requestJson(`${baseUrl}${pathname}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
}

async function waitForHealth() {
  for (let attempt = 0; attempt < 240; attempt += 1) {
    if (serverChild?.exitCode !== null) {
      throw new Error(`server exited during startup: ${childOutput}`);
    }
    try {
      const result = await requestJson(`${baseUrl}/api/health`);
      if (result.response.ok && result.payload.ok) return result.payload;
    } catch {
      // Startup race.
    }
    await wait(500);
  }
  throw new Error(`server did not become healthy: ${childOutput}`);
}

async function startServer() {
  assert(!serverChild, "server is already running");
  childOutput = "";
  const child = spawn(process.execPath, [path.join(__dirname, "server.mjs")], {
    cwd: appRoot,
    env: {
      ...process.env,
      OPEN_MOLECULE_HOST: "127.0.0.1",
      OPEN_MOLECULE_PORT: String(port),
      OPEN_MOLECULE_DATA_DIR: dataRoot,
      OPEN_MOLECULE_RUNS_DIR: runsRoot,
      OPEN_MOLECULE_PYTHON: pythonPath,
      OPEN_MOLECULE_ASSET_ROOT: assetRoot,
    },
    stdio: ["ignore", "pipe", "pipe"],
  });
  child.stdout.on("data", (chunk) => { childOutput += chunk.toString(); });
  child.stderr.on("data", (chunk) => { childOutput += chunk.toString(); });
  serverChild = child;
  return waitForHealth();
}

async function waitForChildExit(child, timeoutMs = 60_000) {
  if (child.exitCode !== null || child.signalCode !== null) return;
  let timeout;
  try {
    await Promise.race([
      new Promise((resolve) => child.once("exit", resolve)),
      new Promise((_, reject) => {
        timeout = setTimeout(() => reject(new Error("server did not stop")), timeoutMs);
      }),
    ]);
  } finally {
    if (timeout) clearTimeout(timeout);
  }
}

async function stopServer(signal = "SIGTERM") {
  const child = serverChild;
  if (!child) return;
  child.kill(signal);
  try {
    await waitForChildExit(child);
  } catch (error) {
    child.kill("SIGKILL");
    await waitForChildExit(child, 10_000).catch(() => {});
    throw error;
  } finally {
    serverChild = null;
  }
}

async function waitForTerminal(runId, timeoutMinutes = 15) {
  const attempts = timeoutMinutes * 60;
  for (let attempt = 0; attempt < attempts; attempt += 1) {
    const result = await requestJson(`${baseUrl}/api/runs/${runId}`);
    if (["complete", "failed", "blocked", "cancelled"].includes(result.payload.status)) return result.payload;
    await wait(1000);
  }
  throw new Error(`run did not finish within ${timeoutMinutes} minutes: ${runId}`);
}

async function waitForStageCheckpoint(runId, stage, status = "complete", timeoutMinutes = 10) {
  const checkpointPath = path.join(
    runsRoot,
    runId,
    "stages",
    stage,
    "attempt-0001",
    "checkpoint.json",
  );
  const attempts = timeoutMinutes * 60 * 4;
  for (let attempt = 0; attempt < attempts; attempt += 1) {
    try {
      const checkpoint = JSON.parse(await fs.readFile(checkpointPath, "utf8"));
      if (checkpoint.status === status) return { checkpoint, checkpointPath };
      if (["failed", "blocked", "cancelled"].includes(checkpoint.status) && checkpoint.status !== status) {
        throw new Error(`${stage} reached ${checkpoint.status}: ${JSON.stringify(checkpoint.error)}`);
      }
    } catch (error) {
      if (error?.code !== "ENOENT" && !String(error?.message || "").includes("Unexpected end")) throw error;
    }
    await wait(250);
  }
  throw new Error(`${stage} did not reach ${status}: ${runId}`);
}

async function getStages(runId) {
  const result = await requestJson(`${baseUrl}/api/runs/${runId}/stages`);
  assert(result.response.ok, `stage summary failed: ${JSON.stringify(result.payload)}`);
  return result.payload;
}

function assertStageChain(summary, expected) {
  const actual = summary.stages.map((row) => `${row.stage}:${row.status}`).join(",");
  assert(actual === expected, `stage chain mismatch: ${actual}`);
}

async function sha256File(filePath) {
  return createHash("sha256").update(await fs.readFile(filePath)).digest("hex");
}

async function collectFiles(root, relative = "") {
  const entries = await fs.readdir(path.join(root, relative), { withFileTypes: true });
  const names = [];
  for (const entry of entries) {
    const child = path.posix.join(relative.split(path.sep).join("/"), entry.name);
    if (child === "MANIFEST.sha256") continue;
    if (entry.isDirectory()) names.push(...await collectFiles(root, child));
    else if (entry.isFile()) names.push(child);
  }
  return names;
}

async function refreshManifest(runRoot) {
  const lines = [];
  for (const relative of (await collectFiles(runRoot)).sort()) {
    lines.push(`${await sha256File(path.join(runRoot, relative))}  ${relative}`);
  }
  await fs.writeFile(path.join(runRoot, "MANIFEST.sha256"), `${lines.join("\n")}\n`, "utf8");
}

async function verifyManifest(runRoot) {
  const manifest = await fs.readFile(path.join(runRoot, "MANIFEST.sha256"), "utf8");
  for (const line of manifest.trim().split("\n")) {
    const [expected, ...nameParts] = line.split("  ");
    const relative = nameParts.join("  ");
    const actual = await sha256File(path.join(runRoot, relative));
    assert(actual === expected, `manifest mismatch: ${relative}`);
  }
}

async function verifyNoHostPaths(runRoot) {
  const textSuffixes = new Set([".json", ".jsonl", ".log", ".md", ".txt"]);
  for (const name of await collectFiles(runRoot)) {
    if (!textSuffixes.has(path.extname(name))) continue;
    const text = await fs.readFile(path.join(runRoot, name), "utf8");
    assert(!text.includes(sourceRoot), `evidence contains source root: ${name}`);
    assert(!text.includes(assetRoot), `evidence contains asset root: ${name}`);
    assert(!text.includes(runRoot), `evidence contains run root: ${name}`);
  }
}

async function attemptCount(runId, stage) {
  const stageRoot = path.join(runsRoot, runId, "stages", stage);
  const entries = await fs.readdir(stageRoot, { withFileTypes: true }).catch(() => []);
  return entries.filter((entry) => entry.isDirectory() && /^attempt-\d{4}$/.test(entry.name)).length;
}

async function createPlanAndRun({ routeMode, moleculeSetId, candidatePool, label }) {
  const plan = await post("/api/prompt-plan", {
    prompt: label,
    target: "CHEMBL2051",
    candidatePool,
    finalSelectionCount: 1,
    routeMode,
  });
  assert(plan.response.ok, `plan failed: ${JSON.stringify(plan.payload)}`);
  const launch = await post("/api/runs", {
    planRunId: plan.payload.runId,
    moleculeSetId,
  });
  assert(launch.response.ok, `run launch failed: ${JSON.stringify(launch.payload)}`);
  return { plan: plan.payload, launch: launch.payload };
}

function compareResultCsv(leftPath, rightPath) {
  const script = String.raw`
import json
import numpy as np
import pandas as pd
import sys

left = pd.read_csv(sys.argv[1]).sort_values("id").reset_index(drop=True)
right = pd.read_csv(sys.argv[2]).sort_values("id").reset_index(drop=True)
mismatches = []
if list(left.columns) != list(right.columns):
    mismatches.append({"columns": [list(left.columns), list(right.columns)]})
elif list(left["id"].astype(str)) != list(right["id"].astype(str)):
    mismatches.append({"ids": [list(left["id"].astype(str)), list(right["id"].astype(str))]})
else:
    for column in left.columns:
        if pd.api.types.is_numeric_dtype(left[column]) and pd.api.types.is_numeric_dtype(right[column]):
            if not np.allclose(left[column].to_numpy(float), right[column].to_numpy(float), atol=1e-4, rtol=0, equal_nan=True):
                mismatches.append(column)
        elif list(left[column].fillna("").astype(str)) != list(right[column].fillna("").astype(str)):
            mismatches.append(column)
print(json.dumps({"ok": not mismatches, "mismatches": mismatches}))
`;
  const result = spawnSync(pythonPath, ["-c", script, leftPath, rightPath], {
    cwd: sourceRoot,
    encoding: "utf8",
    timeout: 60_000,
  });
  assert(result.status === 0, `result comparison failed: ${result.stderr}`);
  return JSON.parse(String(result.stdout || "{}").trim());
}

function tamperFinalScore(scorePath) {
  const script = [
    "import pandas as pd,sys",
    "p=sys.argv[1]",
    "f=pd.read_csv(p)",
    "f.loc[0,'final_score']=float(f.loc[0,'final_score'])+0.01",
    "f.to_csv(p,index=False)",
  ].join(";");
  const result = spawnSync(pythonPath, ["-c", script, scorePath], {
    cwd: sourceRoot,
    encoding: "utf8",
    timeout: 30_000,
  });
  assert(result.status === 0, `score tamper fixture failed: ${result.stderr}`);
}

function assertNoRunSmina(runId) {
  const processList = spawnSync("ps", ["-axo", "command="], { encoding: "utf8", timeout: 5000 });
  assert(processList.status === 0, "could not inspect local processes");
  const runRoot = path.join(runsRoot, runId);
  const active = String(processList.stdout || "")
    .split("\n")
    .filter((line) => line.includes(runRoot) && /smina/i.test(line));
  assert(active.length === 0, `tampered resume started smina: ${active.join("\n")}`);
}

try {
  await startServer();
  const aliasPreflight = await runPreflight({
    sourceRoot,
    pythonPath,
    assetRoot,
    routeBranch: "cascade",
    targetId: "neuraminidase",
    expectedCandidateCount: 2,
    actualCandidateCount: 2,
    sminaBin: process.env.SMINA_BIN || "",
    obabelBin: process.env.OBABEL_BIN || "",
  });
  assert(aliasPreflight.ok, `alias cascade preflight failed: ${JSON.stringify(aliasPreflight.checks)}`);
  assert(
    aliasPreflight.checks.find((check) => check.id === "registered-receptor")?.status === "passed",
    "alias receptor was not resolved",
  );
  assert(
    aliasPreflight.checks.find((check) => check.id === "runtime-versions")?.status === "passed",
    "manifest runtime versions were not verified",
  );
  const modelStatus = await requestJson(`${baseUrl}/api/model-status`);
  assert(modelStatus.response.ok, "model status request failed");
  for (const layer of ["L1", "L2", "L3", "L4", "Dock"]) {
    const row = modelStatus.payload.models.find((item) => item.layer === layer);
    assert(row?.status === "available", `${layer} status is not available: ${JSON.stringify(row)}`);
  }

  const libraryCsv = "id,smiles\nethanol,CCO\naspirin,CC(=O)Oc1ccccc1C(=O)O\n";
  const librarySet = await post("/api/molecule-sets", {
    name: "real-library-smoke",
    csvText: libraryCsv,
    license: "test-only",
  });
  assert(librarySet.response.ok, `library molecule set failed: ${JSON.stringify(librarySet.payload)}`);
  const libraryRun = await createPlanAndRun({
    routeMode: "library",
    moleculeSetId: librarySet.payload.moleculeSetId,
    candidatePool: 2,
    label: "Run a strict real two-molecule library screen for neuraminidase.",
  });
  const libraryTerminal = await waitForTerminal(libraryRun.launch.runId);
  assert(
    libraryTerminal.status === "complete",
    `library run ended as ${libraryTerminal.status}: ${JSON.stringify(libraryTerminal.error)}`,
  );
  const libraryStages = await getStages(libraryTerminal.runId);
  assertStageChain(libraryStages, "prepare:complete,score:complete,dock:skipped,report:complete");
  assert(libraryStages.stages.every((row) => row.stage === "dock" || row.attempts === 1), "library attempts changed");
  const libraryResults = await requestJson(`${baseUrl}/api/runs/${libraryTerminal.runId}/results`);
  assert(libraryResults.response.ok && libraryResults.payload.nRows === 2, "library result rows missing");
  await verifyManifest(path.join(runsRoot, libraryTerminal.runId));
  await verifyNoHostPaths(path.join(runsRoot, libraryTerminal.runId));

  const cascadeCsv = [
    "id,smiles",
    "mol-1,CCO",
    "mol-2,CCN",
    "mol-3,CCC",
    "mol-4,CCCl",
    "mol-5,CCBr",
    "mol-6,CCOC",
    "mol-7,CC(=O)O",
    "mol-8,c1ccccc1",
  ].join("\n") + "\n";
  const cascadeSet = await post("/api/molecule-sets", {
    name: "real-cascade-resume",
    csvText: cascadeCsv,
    license: "test-only",
  });
  assert(cascadeSet.response.ok, `cascade molecule set failed: ${JSON.stringify(cascadeSet.payload)}`);

  const baseline = await createPlanAndRun({
    routeMode: "cascade",
    moleculeSetId: cascadeSet.payload.moleculeSetId,
    candidatePool: 8,
    label: "Run an uninterrupted eight-molecule docking cascade for neuraminidase.",
  });
  const baselineTerminal = await waitForTerminal(baseline.launch.runId);
  assert(
    baselineTerminal.status === "complete",
    `baseline cascade ended as ${baselineTerminal.status}: ${JSON.stringify(baselineTerminal.error)}`,
  );
  const baselineStages = await getStages(baselineTerminal.runId);
  assertStageChain(baselineStages, "prepare:complete,score:complete,dock:complete,report:complete");
  assert(baselineStages.stages.every((row) => row.attempts === 1), "baseline cascade attempt count changed");
  const baselineResults = await requestJson(`${baseUrl}/api/runs/${baselineTerminal.runId}/results`);
  assert(baselineResults.response.ok, `baseline results failed: ${JSON.stringify(baselineResults.payload)}`);
  assert(baselineResults.payload.rankingScoreField === "final_score_dock", "baseline ignored fused score");
  assert(
    baselineResults.payload.rows.some((row) => row.structure_docking_status === "ok"),
    "baseline cascade had no successful docking row",
  );

  const interrupted = await createPlanAndRun({
    routeMode: "cascade",
    moleculeSetId: cascadeSet.payload.moleculeSetId,
    candidatePool: 8,
    label: "Interrupt and resume an eight-molecule docking cascade for neuraminidase.",
  });
  const interruptedScore = await waitForStageCheckpoint(interrupted.launch.runId, "score");
  const scoreOutput = interruptedScore.checkpoint.outputs["scores.csv"];
  const scorePath = path.join(path.dirname(interruptedScore.checkpointPath), scoreOutput.path);
  const scoreShaBefore = await sha256File(scorePath);
  await stopServer("SIGTERM");
  await startServer();
  const interruptedRun = await requestJson(`${baseUrl}/api/runs/${interrupted.launch.runId}`);
  assert(
    ["failed", "cancelled"].includes(interruptedRun.payload.status),
    `interrupted run became ${interruptedRun.payload.status}`,
  );
  const interruptedBeforeResume = await getStages(interrupted.launch.runId);
  assert(
    interruptedBeforeResume.stages.find((row) => row.stage === "score")?.status === "complete",
    "restart lost complete score attempt",
  );
  const resume = await post(`/api/runs/${interrupted.launch.runId}/resume`);
  assert(resume.response.status === 202, `resume returned HTTP ${resume.response.status}: ${JSON.stringify(resume.payload)}`);
  const resumedTerminal = await waitForTerminal(interrupted.launch.runId);
  assert(
    resumedTerminal.status === "complete",
    `resumed cascade ended as ${resumedTerminal.status}: ${JSON.stringify(resumedTerminal.error)}`,
  );
  assert(await sha256File(scorePath) === scoreShaBefore, "resume overwrote verified score output");
  const resumedStages = await getStages(interrupted.launch.runId);
  assertStageChain(resumedStages, "prepare:complete,score:complete,dock:complete,report:complete");
  assert(resumedStages.stages.find((row) => row.stage === "score")?.attempts === 1, "resume reran score");
  assert(resumedStages.stages.find((row) => row.stage === "dock")?.attempts >= 1, "resume did not run dock");
  const parity = compareResultCsv(
    path.join(runsRoot, baselineTerminal.runId, "results", "scores.csv"),
    path.join(runsRoot, resumedTerminal.runId, "results", "scores.csv"),
  );
  assert(parity.ok, `resumed result parity failed: ${JSON.stringify(parity.mismatches)}`);
  await verifyManifest(path.join(runsRoot, resumedTerminal.runId));
  await verifyNoHostPaths(path.join(runsRoot, resumedTerminal.runId));

  const tampered = await createPlanAndRun({
    routeMode: "cascade",
    moleculeSetId: cascadeSet.payload.moleculeSetId,
    candidatePool: 8,
    label: "Tamper with an interrupted score checkpoint and require fail-closed resume.",
  });
  const tamperedScore = await waitForStageCheckpoint(tampered.launch.runId, "score");
  const tamperedOutput = tamperedScore.checkpoint.outputs["scores.csv"];
  const tamperedScorePath = path.join(path.dirname(tamperedScore.checkpointPath), tamperedOutput.path);
  await stopServer("SIGTERM");
  tamperFinalScore(tamperedScorePath);
  const tamperedRunRoot = path.join(runsRoot, tampered.launch.runId);
  await refreshManifest(tamperedRunRoot);
  const dockAttemptsBefore = await attemptCount(tampered.launch.runId, "dock");
  await startServer();
  const blockedResume = await post(`/api/runs/${tampered.launch.runId}/resume`);
  assert(blockedResume.response.status === 409, `tampered resume returned HTTP ${blockedResume.response.status}`);
  assert(blockedResume.payload.status === "blocked", "tampered resume was not blocked");
  assert(blockedResume.payload.error?.code === "checkpoint_mismatch", "tampered resume error code changed");
  assert(blockedResume.payload.error?.stage === "score", "tampered resume did not identify score stage");
  assert(
    await attemptCount(tampered.launch.runId, "dock") === dockAttemptsBefore,
    "tampered resume created a new dock attempt",
  );
  assertNoRunSmina(tampered.launch.runId);

  process.stdout.write(`${JSON.stringify({
    ok: true,
    roots: { dataRoot, runsRoot },
    library: { runId: libraryTerminal.runId, stages: libraryStages.stages },
    baseline: { runId: baselineTerminal.runId, stages: baselineStages.stages },
    resumed: { runId: resumedTerminal.runId, stages: resumedStages.stages, parity },
    tampered: { runId: tampered.launch.runId, status: blockedResume.payload.status },
  }, null, 2)}\n`);
} finally {
  await stopServer("SIGTERM").catch(() => {});
  if (process.env.OPEN_MOLECULE_KEEP_REAL_RUN !== "1") {
    await fs.rm(dataRoot, { recursive: true, force: true });
    await fs.rm(runsRoot, { recursive: true, force: true });
  } else {
    process.stderr.write(`kept real run roots: ${dataRoot} ${runsRoot}\n`);
  }
}
