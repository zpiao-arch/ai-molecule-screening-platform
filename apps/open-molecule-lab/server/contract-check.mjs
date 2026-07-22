import { spawn } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const appRoot = path.resolve(__dirname, "..");
const sourceRoot = path.resolve(appRoot, "..", "..");
const runsRoot = await fs.mkdtemp(path.join(os.tmpdir(), "open-molecule-lab-contract-"));
const port = 32000 + (process.pid % 10000);
const baseUrl = `http://127.0.0.1:${port}`;
const child = spawn(process.execPath, [path.join(__dirname, "server.mjs")], {
  cwd: appRoot,
  env: {
    ...process.env,
    OPEN_MOLECULE_HOST: "127.0.0.1",
    OPEN_MOLECULE_PORT: String(port),
    OPEN_MOLECULE_RUNS_DIR: runsRoot,
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
}
