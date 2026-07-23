import { createHash } from "node:crypto";
import { spawn } from "node:child_process";
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

let childOutput = "";
child.stdout.on("data", (chunk) => { childOutput += chunk.toString(); });
child.stderr.on("data", (chunk) => { childOutput += chunk.toString(); });

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

async function requestJson(url, options) {
  const response = await fetch(url, options);
  const payload = await response.json();
  return { response, payload };
}

async function post(pathname, body) {
  return requestJson(`${baseUrl}${pathname}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
}

async function waitForHealth() {
  for (let attempt = 0; attempt < 100; attempt += 1) {
    try {
      const result = await requestJson(`${baseUrl}/api/health`);
      if (result.response.ok && result.payload.ok) return;
    } catch {
      // Startup race.
    }
    await new Promise((resolve) => setTimeout(resolve, 100));
  }
  throw new Error(`server did not become healthy: ${childOutput}`);
}

async function waitForTerminal(runId) {
  for (let attempt = 0; attempt < 600; attempt += 1) {
    const result = await requestJson(`${baseUrl}/api/runs/${runId}`);
    if (["complete", "failed", "blocked", "cancelled"].includes(result.payload.status)) return result.payload;
    await new Promise((resolve) => setTimeout(resolve, 1000));
  }
  throw new Error(`run did not finish within 10 minutes: ${runId}`);
}

async function verifyManifest(runRoot) {
  const manifest = await fs.readFile(path.join(runRoot, "MANIFEST.sha256"), "utf8");
  for (const line of manifest.trim().split("\n")) {
    const [expected, ...nameParts] = line.split("  ");
    const relative = nameParts.join("  ");
    const bytes = await fs.readFile(path.join(runRoot, relative));
    const actual = createHash("sha256").update(bytes).digest("hex");
    assert(actual === expected, `manifest mismatch: ${relative}`);
  }
}

async function verifyNoHostPaths(runRoot) {
  const names = [
    "run.json",
    "run-spec.json",
    "preflight.json",
    "command.json",
    "events.jsonl",
    "logs/stdout.log",
    "logs/stderr.log",
  ];
  for (const name of names) {
    const text = await fs.readFile(path.join(runRoot, name), "utf8");
    assert(!text.includes(sourceRoot), `evidence contains source root: ${name}`);
    assert(!text.includes(assetRoot), `evidence contains asset root: ${name}`);
    assert(!text.includes(runRoot), `evidence contains run root: ${name}`);
  }
}

try {
  await waitForHealth();
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
  assert(aliasPreflight.checks.find((check) => check.id === "registered-receptor")?.status === "passed", "alias receptor was not resolved");
  assert(aliasPreflight.checks.find((check) => check.id === "runtime-versions")?.status === "passed", "manifest runtime versions were not verified");
  assert(aliasPreflight.checks.find((check) => check.id === "python-modules")?.details?.checked?.includes("joblib"), "joblib was not included in strict module preflight");
  const modelStatus = await requestJson(`${baseUrl}/api/model-status`);
  assert(modelStatus.response.ok, "model status request failed");
  for (const layer of ["L1", "L2", "L3", "L4", "Dock"]) {
    const row = modelStatus.payload.models.find((item) => item.layer === layer);
    assert(row?.status === "available", `${layer} status is not available: ${JSON.stringify(row)}`);
  }
  const csvText = "id,smiles\nethanol,CCO\naspirin,CC(=O)Oc1ccccc1C(=O)O\n";
  const moleculeSet = await post("/api/molecule-sets", { name: "real-run-smoke", csvText, license: "test-only" });
  assert(moleculeSet.response.ok, `molecule set failed: ${JSON.stringify(moleculeSet.payload)}`);
  const autoPlan = await post("/api/prompt-plan", {
    prompt: "Resolve the available route for neuraminidase.",
    target: "CHEMBL2051",
    candidatePool: 2,
    finalSelectionCount: 1,
    routeMode: "auto",
  });
  assert(autoPlan.payload.route?.branch === "cascade", `auto route ignored external receptor: ${JSON.stringify(autoPlan.payload.route)}`);
  const plan = await post("/api/prompt-plan", {
    prompt: "Run a strict real two-molecule library screen for neuraminidase.",
    target: "CHEMBL2051",
    candidatePool: 2,
    finalSelectionCount: 1,
    routeMode: "library",
  });
  assert(plan.response.ok, `plan failed: ${JSON.stringify(plan.payload)}`);
  const launch = await post("/api/runs", {
    planRunId: plan.payload.runId,
    moleculeSetId: moleculeSet.payload.moleculeSetId,
  });
  assert(launch.response.ok, `run launch failed: ${JSON.stringify(launch.payload)}`);
  const terminal = await waitForTerminal(launch.payload.runId);
  assert(terminal.status === "complete", `real run ended as ${terminal.status}: ${JSON.stringify(terminal.error)}`);
  const results = await requestJson(`${baseUrl}/api/runs/${terminal.runId}/results`);
  assert(results.response.ok, `results failed: ${JSON.stringify(results.payload)}`);
  assert(results.payload.nRows === 2, "real result row count mismatch");
  assert(results.payload.rows.length === 2, "real result rows missing");
  for (const row of results.payload.rows) {
    for (let layer = 1; layer <= 4; layer += 1) {
      assert(row[`layer${layer}_status`] === "ok", `${row.id} L${layer} did not complete`);
    }
  }
  await verifyManifest(path.join(runsRoot, terminal.runId));
  await verifyNoHostPaths(path.join(runsRoot, terminal.runId));

  const cascadePlan = await post("/api/prompt-plan", {
    prompt: "Run a strict real two-molecule docking cascade for neuraminidase.",
    target: "CHEMBL2051",
    candidatePool: 2,
    finalSelectionCount: 1,
    routeMode: "cascade",
  });
  assert(cascadePlan.response.ok, `cascade plan failed: ${JSON.stringify(cascadePlan.payload)}`);
  assert(cascadePlan.payload.route?.branch === "cascade", "cascade plan resolved to the wrong branch");
  const cascadeLaunch = await post("/api/runs", {
    planRunId: cascadePlan.payload.runId,
    moleculeSetId: moleculeSet.payload.moleculeSetId,
  });
  assert(cascadeLaunch.response.ok, `cascade launch failed: ${JSON.stringify(cascadeLaunch.payload)}`);
  const cascadeTerminal = await waitForTerminal(cascadeLaunch.payload.runId);
  assert(
    cascadeTerminal.status === "complete",
    `cascade run ended as ${cascadeTerminal.status}: ${JSON.stringify(cascadeTerminal.error)}`,
  );
  const cascadeResults = await requestJson(`${baseUrl}/api/runs/${cascadeTerminal.runId}/results`);
  assert(cascadeResults.response.ok, `cascade results failed: ${JSON.stringify(cascadeResults.payload)}`);
  assert(
    cascadeResults.payload.rankingScoreField === "final_score_dock",
    `cascade results ignored fused score: ${JSON.stringify(cascadeResults.payload)}`,
  );
  assert(
    cascadeResults.payload.rows.every((row) => typeof row.final_score_dock === "number"),
    "cascade results omitted fused scores",
  );
  assert(
    cascadeResults.payload.rows.some((row) => row.structure_docking_status === "ok"),
    "cascade run had no successful structure docking row",
  );
  await verifyManifest(path.join(runsRoot, cascadeTerminal.runId));
  await verifyNoHostPaths(path.join(runsRoot, cascadeTerminal.runId));
  process.stdout.write(`${JSON.stringify({
    ok: true,
    library: { runId: terminal.runId, status: terminal.status, results: results.payload },
    cascade: { runId: cascadeTerminal.runId, status: cascadeTerminal.status, results: cascadeResults.payload },
  }, null, 2)}\n`);
} finally {
  child.kill("SIGTERM");
  if (process.env.OPEN_MOLECULE_KEEP_REAL_RUN !== "1") {
    await fs.rm(dataRoot, { recursive: true, force: true });
    await fs.rm(runsRoot, { recursive: true, force: true });
  } else {
    process.stderr.write(`kept real run roots: ${dataRoot} ${runsRoot}\n`);
  }
}
