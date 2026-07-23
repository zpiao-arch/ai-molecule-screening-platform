import { createHash } from "node:crypto";
import { spawn, spawnSync } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const appRoot = path.resolve(__dirname, "..");
const sourceRoot = path.resolve(appRoot, "..", "..");
const workspacePython = path.resolve(sourceRoot, "..", "..", ".venv_four_level_cli", "bin", "python");
const runtimePython = process.env.OPEN_MOLECULE_PYTHON || workspacePython;
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
    OPEN_MOLECULE_PYTHON: runtimePython,
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

async function collectManifestFiles(root, relative = "") {
  const entries = await fs.readdir(path.join(root, relative), { withFileTypes: true });
  const files = [];
  for (const entry of entries) {
    const child = path.posix.join(relative.split(path.sep).join("/"), entry.name);
    if (child === "MANIFEST.sha256") continue;
    if (entry.isDirectory()) files.push(...await collectManifestFiles(root, child));
    else if (entry.isFile()) files.push(child);
  }
  return files;
}

async function refreshManifest(root) {
  const lines = [];
  for (const relative of (await collectManifestFiles(root)).sort()) {
    const digest = createHash("sha256").update(await fs.readFile(path.join(root, relative))).digest("hex");
    lines.push(`${digest}  ${relative}`);
  }
  await fs.writeFile(path.join(root, "MANIFEST.sha256"), `${lines.join("\n")}\n`, "utf8");
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

async function postResume(runId) {
  return requestJson(`${baseUrl}/api/runs/${runId}/resume`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: "{}",
  });
}

async function getStages(runId) {
  return requestJson(`${baseUrl}/api/runs/${runId}/stages`);
}

try {
  const health = await waitForHealth();
  assert(health.product === "open-molecule-lab", `unexpected product: ${health.product}`);
  assert(health.mode === "local_execution", `unexpected mode: ${health.mode}`);
  assert(health.cli?.available === true, "four-level CLI source entries must be discoverable");
  assert(!JSON.stringify(health).includes(sourceRoot), "health response must not expose the host source root");

  const bridgeFixture = path.join(dataRoot, "bridge-results.csv");
  const bridgeInput = path.join(dataRoot, "bridge-input.csv");
  const bridgeRows = ["id,smiles,layer1_status,layer2_status,layer3_status,layer4_status,final_score,final_score_dock"];
  const bridgeInputRows = ["id,smiles"];
  for (let index = 0; index < 250; index += 1) {
    bridgeRows.push(`mol-${index},CCO,ok,ok,ok,ok,${(index / 250).toFixed(4)},${((249 - index) / 250).toFixed(4)}`);
    bridgeInputRows.push(`mol-${index},CCO`);
  }
  await fs.writeFile(bridgeFixture, `${bridgeRows.join("\n")}\n`, "utf8");
  await fs.writeFile(bridgeInput, `${bridgeInputRows.join("\n")}\n`, "utf8");
  const summaryBridge = spawnSync(runtimePython, [path.join(__dirname, "python-bridge.py"), "summarize-results", "--input", bridgeFixture, "--expected-input", bridgeInput], { encoding: "utf8" });
  assert(summaryBridge.status === 0, `summary bridge failed: ${summaryBridge.stderr}`);
  const summaryFixture = JSON.parse(summaryBridge.stdout);
  assert(summaryFixture.nRows === 250 && !Object.hasOwn(summaryFixture, "rows"), "summary bridge must not materialize all result rows");
  const missingRowFixture = path.join(dataRoot, "bridge-results-missing-row.csv");
  await fs.writeFile(missingRowFixture, `${bridgeRows.slice(0, -1).join("\n")}\n`, "utf8");
  const missingRowBridge = spawnSync(runtimePython, [path.join(__dirname, "python-bridge.py"), "summarize-results", "--input", missingRowFixture, "--expected-input", bridgeInput], { encoding: "utf8" });
  assert(missingRowBridge.status !== 0, "summary bridge accepted a result file with a missing input row");
  const pageBridge = spawnSync(runtimePython, [path.join(__dirname, "python-bridge.py"), "page-results", "--input", bridgeFixture, "--offset", "20", "--limit", "25", "--view", "ranked"], { encoding: "utf8" });
  assert(pageBridge.status === 0, `page bridge failed: ${pageBridge.stderr}`);
  const pageFixture = JSON.parse(pageBridge.stdout);
  assert(pageFixture.rows.length === 25 && pageFixture.offset === 20 && pageFixture.limit === 25, "page bridge did not honor bounds");
  assert(pageFixture.rankingScoreField === "final_score_dock", "page bridge ignored cascade fused scores");
  assert(pageFixture.rows[0]?.id === "mol-20", "ranked page was not ordered by cascade fused score");

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
  assert(Array.isArray(payload.stages) && payload.stages.length === 4, "expected four planned stages");
  assert(
    payload.stages.map((stage) => stage.id).join(",") === "prepare,score,dock,report",
    "planned stage order must match the persisted StageAttempt chain",
  );
  assert(payload.bundle?.files?.includes("MANIFEST.sha256"), "bundle manifest missing");
  assert(payload.equivalentCli.includes("--mode library"), "equivalent CLI omitted resolved mode");
  assert(payload.equivalentCli.includes("--asset-root"), "equivalent CLI omitted asset root");

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
  const blockedStages = await getStages(execution.payload.runId);
  assert(blockedStages.response.status === 200, `blocked stage summary returned HTTP ${blockedStages.response.status}`);
  assert(
    blockedStages.payload.schemaVersion === "open-molecule-lab.stage-summary.v0.1",
    "stage summary schema missing",
  );
  assert(blockedStages.payload.resumable === false, "blocked run was marked resumable");
  assert(
    blockedStages.payload.stages.map((row) => `${row.stage}:${row.status}`).join(",")
      === "prepare:waiting,score:waiting,dock:skipped,report:waiting",
    "blocked library stage summary was inconsistent",
  );
  assert(!JSON.stringify(blockedStages.payload).includes(runsRoot), "stage summary exposed a host path");
  const blockedResume = await postResume(execution.payload.runId);
  assert(blockedResume.response.status === 409, "blocked run resume must return HTTP 409");
  assert(blockedResume.payload.error?.code === "run_not_resumable", "blocked resume error code mismatch");
  const missingStages = await getStages("run_deadbeef_missing");
  assert(missingStages.response.status === 404, "unknown stage summary must return HTTP 404");
  const missingResume = await postResume("run_deadbeef_missing");
  assert(missingResume.response.status === 404, "unknown run resume must return HTTP 404");

  const failedRunPath = path.join(executionRoot, "run.json");
  const failedRun = JSON.parse(await fs.readFile(failedRunPath, "utf8"));
  await fs.writeFile(
    failedRunPath,
    `${JSON.stringify({
      ...failedRun,
      status: "failed",
      finishedAt: new Date().toISOString(),
      error: { code: "worker_interrupted", message: "fixture interruption" },
    }, null, 2)}\n`,
    "utf8",
  );
  await refreshManifest(executionRoot);
  const failedResume = await postResume(execution.payload.runId);
  assert(failedResume.response.status === 409, "failed source-only resume must return HTTP 409");
  assert(failedResume.payload.status === "blocked", "failed fresh preflight did not block resume");
  assert(failedResume.payload.error?.code === "resume_preflight_failed", "resume preflight error code mismatch");

  const pagedRunRoot = path.join(runsRoot, execution.payload.runId);
  const pagedRunPath = path.join(pagedRunRoot, "run.json");
  const pagedRun = JSON.parse(await fs.readFile(pagedRunPath, "utf8"));
  await fs.writeFile(
    pagedRunPath,
    `${JSON.stringify({
      ...pagedRun,
      status: "complete",
      finishedAt: new Date().toISOString(),
      error: null,
      resultSummary: { nRows: 250, nRanked: 249, nFailed: 1 },
    }, null, 2)}\n`,
    "utf8",
  );
  const pagedRows = [...bridgeRows];
  pagedRows[1] = "mol-0,CCO,failed,ok,ok,ok,,";
  await fs.writeFile(path.join(pagedRunRoot, "results", "scores.csv"), `${pagedRows.join("\n")}\n`, "utf8");
  await fs.writeFile(
    path.join(pagedRunRoot, "results", "summary.json"),
    `${JSON.stringify({ ok: true, nRows: 250, nRanked: 249, nFailed: 1, columns: bridgeRows[0].split(","), rankingScoreField: "final_score_dock" }, null, 2)}\n`,
    "utf8",
  );
  const tamperedPage = await requestJson(
    `${baseUrl}/api/runs/${execution.payload.runId}/results?view=ranked&offset=0&limit=25`,
  );
  assert(tamperedPage.response.status === 409, "results API accepted files that do not match MANIFEST.sha256");
  await refreshManifest(pagedRunRoot);
  const completeResume = await postResume(execution.payload.runId);
  assert(completeResume.response.status === 409, "complete run resume must return HTTP 409");
  assert(completeResume.payload.error?.code === "run_not_resumable", "complete resume error code mismatch");
  const rankedPage = await requestJson(
    `${baseUrl}/api/runs/${execution.payload.runId}/results?view=ranked&offset=20&limit=25`,
  );
  assert(rankedPage.response.status === 200, `ranked results page returned HTTP ${rankedPage.response.status}`);
  assert(rankedPage.payload.view === "ranked", "ranked results page omitted view");
  assert(rankedPage.payload.offset === 20 && rankedPage.payload.limit === 25, "ranked results page bounds mismatch");
  assert(rankedPage.payload.total === 249 && rankedPage.payload.rows.length === 25, "ranked results page total mismatch");
  assert(!Object.hasOwn(rankedPage.payload, "rankedRows"), "paged API must not materialize rankedRows");
  const failedPage = await requestJson(
    `${baseUrl}/api/runs/${execution.payload.runId}/results?view=failed&offset=0&limit=25`,
  );
  assert(failedPage.response.status === 200, `failed results page returned HTTP ${failedPage.response.status}`);
  assert(failedPage.payload.view === "failed" && failedPage.payload.total === 1, "failed results page total mismatch");
  assert(failedPage.payload.rows[0]?.id === "mol-0", "failed results page returned the wrong row");
  const invalidPage = await requestJson(
    `${baseUrl}/api/runs/${execution.payload.runId}/results?view=ranked&offset=-1&limit=25`,
  );
  assert(invalidPage.response.status === 400, "invalid results pagination must return HTTP 400");

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
