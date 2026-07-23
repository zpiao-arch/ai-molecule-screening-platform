import { spawn } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const appRoot = path.resolve(__dirname, "..");
const sourceRoot = path.resolve(appRoot, "..", "..");
const workspacePython = path.resolve(sourceRoot, "..", "..", ".venv_four_level_cli", "bin", "python");
const runsRoot = await fs.mkdtemp(path.join(os.tmpdir(), "open-molecule-lab-contract-"));
const dataRoot = await fs.mkdtemp(path.join(os.tmpdir(), "open-molecule-lab-data-"));
const moleculeSetsRoot = path.join(dataRoot, "molecule_sets");
const port = 32000 + (process.pid % 10000);
const baseUrl = `http://127.0.0.1:${port}`;
const child = spawn(process.execPath, [path.join(__dirname, "server.mjs")], {
  cwd: appRoot,
  env: {
    ...process.env,
    OPEN_MOLECULE_HOST: "127.0.0.1",
    OPEN_MOLECULE_PORT: String(port),
    OPEN_MOLECULE_DATA_DIR: dataRoot,
    OPEN_MOLECULE_RUNS_DIR: runsRoot,
    OPEN_MOLECULE_PYTHON: process.env.OPEN_MOLECULE_PYTHON || workspacePython,
    OPEN_MOLECULE_MAX_BODY_BYTES: "4096",
  },
  stdio: ["ignore", "pipe", "pipe"],
});

let childOutput = "";
child.stdout.on("data", (chunk) => {
  childOutput += chunk.toString();
});
child.stderr.on("data", (chunk) => {
  childOutput += chunk.toString();
});

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

async function requestJson(url, options) {
  const response = await fetch(url, options);
  const payload = await response.json();
  return { response, payload };
}

async function waitForHealth() {
  for (let attempt = 0; attempt < 60; attempt += 1) {
    try {
      const { response, payload } = await requestJson(`${baseUrl}/api/health`);
      if (response.ok && payload.ok) return payload;
    } catch {
      // Server startup race.
    }
    await new Promise((resolve) => setTimeout(resolve, 100));
  }
  throw new Error(`server did not become healthy: ${childOutput}`);
}

async function postPrompt(body) {
  return requestJson(`${baseUrl}/api/prompt-plan`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
}

async function postMoleculeSet(body) {
  return requestJson(`${baseUrl}/api/molecule-sets`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
}

async function postRun(body) {
  return requestJson(`${baseUrl}/api/runs`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
}

async function postCancel(runId) {
  return requestJson(`${baseUrl}/api/runs/${runId}/cancel`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: "{}",
  });
}

try {
  const health = await waitForHealth();
  assert(health.product === "open-molecule-lab", `unexpected product: ${health.product}`);
  assert(health.mode === "plan_only", `unexpected mode: ${health.mode}`);
  assert(health.cli?.available === true, "four-level CLI source entries must be discoverable");

  const { response, payload } = await postPrompt({
    prompt: "Plan an auditable four-level screen for influenza neuraminidase.",
    target: "CHEMBL2051",
    candidatePool: 1000,
    finalSelectionCount: 10,
    routeMode: "auto",
  });
  assert(response.status === 201, `valid prompt returned HTTP ${response.status}`);
  assert(payload.ok === true && payload.mode === "plan_only", "valid prompt must return a plan-only response");
  assert(payload.runId && payload.spec?.schemaVersion === "open-molecule-lab.run-spec.v0.1", "RunSpec missing");
  assert(Array.isArray(payload.stages) && payload.stages.length === 5, "expected five planned stages");
  assert(payload.bundle?.files?.includes("MANIFEST.sha256"), "bundle manifest missing");

  const runRoot = path.join(runsRoot, payload.runId);
  const expectedFiles = ["DESIGN.md", "MANIFEST.sha256", "prompt.txt", "run.json"];
  assert(JSON.stringify((await fs.readdir(runRoot)).sort()) === JSON.stringify(expectedFiles), "unexpected bundle files");
  const bundleText = (
    await Promise.all(expectedFiles.map((name) => fs.readFile(path.join(runRoot, name), "utf8")))
  ).join("\n");
  assert(!bundleText.includes(sourceRoot), "bundle must not contain the absolute workspace path");

  const validCsv = "\ufeffid,smiles\nmol-1,CCO\nmol-2,CCN\n";
  const moleculeSet = await postMoleculeSet({
    name: "contract-fixture",
    csvText: validCsv,
    license: "CC0-test",
  });
  assert(moleculeSet.response.status === 201, `valid molecule set returned HTTP ${moleculeSet.response.status}`);
  assert(moleculeSet.payload.ok === true, "valid molecule set must return ok");
  assert(moleculeSet.payload.moleculeSetId?.startsWith("molset_"), "molecule set ID missing");
  assert(moleculeSet.payload.nRows === 2, "molecule set row count mismatch");
  assert(moleculeSet.payload.inputSha256?.length === 64, "molecule set hash missing");
  const moleculeSetId = moleculeSet.payload.moleculeSetId;
  const moleculeSetRoot = path.join(moleculeSetsRoot, moleculeSetId);
  assert(JSON.stringify((await fs.readdir(moleculeSetRoot)).sort()) === JSON.stringify(["MANIFEST.sha256", "input.csv", "metadata.json"]), "molecule set files mismatch");
  const fetchedSet = await requestJson(`${baseUrl}/api/molecule-sets/${moleculeSetId}`);
  assert(fetchedSet.response.status === 200 && fetchedSet.payload.nRows === 2, "molecule set GET failed");
  const repeated = await postMoleculeSet({ name: "same-bytes", csvText: validCsv, license: "CC0-test" });
  assert(repeated.response.status === 200, "identical molecule set should be idempotent");
  assert(repeated.payload.moleculeSetId === moleculeSetId, "identical molecule set changed ID");

  const beforeInvalidMoleculeSets = (await fs.readdir(moleculeSetsRoot)).length;
  const invalidMoleculeSets = [
    { csvText: "name,smiles\na,CCO\n", field: "columns" },
    { csvText: "id,smiles\na,CCO\na,CCN\n", field: "id" },
    { csvText: "id,smiles\na,\n", field: "smiles" },
    { csvText: "id,smiles\n", field: "rows" },
  ];
  for (const invalid of invalidMoleculeSets) {
    const result = await postMoleculeSet({ name: `invalid-${invalid.field}`, csvText: invalid.csvText });
    assert(result.response.status === 400, `${invalid.field} molecule set returned HTTP ${result.response.status}`);
    assert(result.payload.error?.field === invalid.field, `${invalid.field} error field mismatch`);
  }
  assert((await fs.readdir(moleculeSetsRoot)).length === beforeInvalidMoleculeSets, "invalid molecule sets left artifacts");
  const oversized = await postMoleculeSet({
    name: "oversized",
    csvText: `id,smiles\nmol-1,${"C".repeat(5000)}\n`,
  });
  assert(oversized.response.status === 413, `oversized molecule set returned HTTP ${oversized.response.status}`);
  assert(oversized.payload.error?.code === "payload_too_large", "oversized error code mismatch");

  const executablePlan = await postPrompt({
    prompt: "Run a strict two-molecule library screen for influenza neuraminidase.",
    target: "CHEMBL2051",
    candidatePool: 2,
    finalSelectionCount: 1,
    routeMode: "library",
  });
  assert(executablePlan.response.status === 201, "executable plan creation failed");
  const execution = await postRun({
    planRunId: executablePlan.payload.runId,
    moleculeSetId,
  });
  assert(execution.response.status === 201, `blocked execution returned HTTP ${execution.response.status}`);
  assert(execution.payload.ok === true, "blocked execution must still persist successfully");
  assert(execution.payload.status === "blocked", `source-only execution must block, got ${execution.payload.status}`);
  assert(execution.payload.runId?.startsWith("run_"), "execution run ID missing");
  assert(execution.payload.preflight?.ok === false, "source-only preflight unexpectedly passed");
  assert(
    execution.payload.preflight.checks.some((check) => check.required && check.status === "failed"),
    "blocked preflight must name a failed required check",
  );
  const executionRoot = path.join(runsRoot, execution.payload.runId);
  const executionFiles = (await fs.readdir(executionRoot)).sort();
  assert(executionFiles.includes("run.json"), "execution run.json missing");
  assert(executionFiles.includes("run-spec.json"), "attached RunSpec missing");
  assert(executionFiles.includes("preflight.json"), "preflight evidence missing");
  assert(executionFiles.includes("events.jsonl"), "event log missing");
  assert(executionFiles.includes("MANIFEST.sha256"), "execution manifest missing");
  assert((await fs.readFile(path.join(executionRoot, "inputs", "molecules.csv"), "utf8")) === validCsv, "sealed input changed");
  try {
    await fs.access(path.join(executionRoot, "results", "scores.csv"));
    throw new Error("blocked execution produced scores.csv");
  } catch (error) {
    if (error?.message === "blocked execution produced scores.csv") throw error;
  }
  const fetchedRun = await requestJson(`${baseUrl}/api/runs/${execution.payload.runId}`);
  assert(fetchedRun.response.status === 200 && fetchedRun.payload.status === "blocked", "blocked run GET failed");
  const blockedResults = await requestJson(`${baseUrl}/api/runs/${execution.payload.runId}/results`);
  assert(blockedResults.response.status === 409, "blocked run results must return HTTP 409");
  const blockedCancel = await postCancel(execution.payload.runId);
  assert(blockedCancel.response.status === 409, "blocked run cancellation must return HTTP 409");

  const beforeInvalid = (await fs.readdir(runsRoot)).length;
  const invalidCases = [
    { body: { prompt: "", candidatePool: 1000, finalSelectionCount: 10, routeMode: "auto" }, field: "prompt" },
    { body: { prompt: "x", candidatePool: 0, finalSelectionCount: 1, routeMode: "auto" }, field: "candidatePool" },
    { body: { prompt: "x", candidatePool: 10, finalSelectionCount: 1, routeMode: "magic" }, field: "routeMode" },
  ];
  for (const invalid of invalidCases) {
    const result = await postPrompt(invalid.body);
    assert(result.response.status === 400, `${invalid.field} case returned HTTP ${result.response.status}`);
    assert(result.payload.error?.field === invalid.field, `${invalid.field} error did not preserve field`);
  }
  const afterInvalid = (await fs.readdir(runsRoot)).length;
  assert(beforeInvalid === afterInvalid, "invalid prompts must not create run directories");

  process.stdout.write(
    `${JSON.stringify({ ok: true, mode: health.mode, runId: payload.runId, files: expectedFiles, route: payload.route }, null, 2)}\n`,
  );
} finally {
  child.kill("SIGTERM");
  await fs.rm(runsRoot, { recursive: true, force: true });
  await fs.rm(dataRoot, { recursive: true, force: true });
}
