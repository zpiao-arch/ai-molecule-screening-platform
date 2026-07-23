import { createHash, randomBytes } from "node:crypto";
import { createServer } from "node:http";
import { spawnSync } from "node:child_process";
import { existsSync } from "node:fs";
import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

import { computeStageFingerprint } from "./stage-contract.mjs";
import { createStageStore } from "./stage-store.mjs";
import { atomicWriteJson, createArtifactStore, StoreError } from "./store.mjs";
import { runPreflight } from "./preflight.mjs";
import { createStageInputBuilder, createWorkerManager } from "./worker.mjs";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const appRoot = path.resolve(__dirname, "..");
const sourceRoot = path.resolve(appRoot, "..", "..");
const distRoot = path.join(appRoot, "dist");
const dataRoot = path.resolve(process.env.OPEN_MOLECULE_DATA_DIR || path.join(appRoot, "data"));
const runsRoot = path.resolve(process.env.OPEN_MOLECULE_RUNS_DIR || path.join(appRoot, "runs"));
const pythonPath = process.env.OPEN_MOLECULE_PYTHON || "python3";
const bridgePath = path.join(__dirname, "python-bridge.py");
const externalAssetRoot = process.env.OPEN_MOLECULE_ASSET_ROOT
  ? path.resolve(process.env.OPEN_MOLECULE_ASSET_ROOT)
  : "";
const scoringEntry = path.join(sourceRoot, "scoring", "scoring.py");
const benchmarkEntry = path.join(
  sourceRoot,
  "scientific_validation",
  "four_level_cli_1kx10k",
  "batch_cli.py",
);
const assetManifest = path.join(sourceRoot, "assets", "ASSET_MANIFEST.json");
const receptorRegistry = path.join(sourceRoot, "scoring", "receptor_registry.json");
const host = process.env.OPEN_MOLECULE_HOST || process.env.HOST || "127.0.0.1";
const port = Number(process.env.OPEN_MOLECULE_PORT || process.env.PORT || 4173);
const maxBodyBytes = Number(process.env.OPEN_MOLECULE_MAX_BODY_BYTES || 8 * 1024 * 1024);

function externalAssetBundleVerified() {
  if (!externalAssetRoot) return false;
  const manifestPath = path.join(externalAssetRoot, "ASSET_MANIFEST.json");
  if (!existsSync(manifestPath)) return false;
  const verify = spawnSync(
    pythonPath,
    [
      path.join(sourceRoot, "scripts", "verify_assets.py"),
      "--asset-root",
      externalAssetRoot,
      "--manifest",
      manifestPath,
    ],
    { cwd: sourceRoot, encoding: "utf8", timeout: 120_000 },
  );
  try {
    return verify.status === 0 && JSON.parse(String(verify.stdout || "{}")).ok === true;
  } catch {
    return false;
  }
}

const assetBundleVerified = externalAssetBundleVerified();
const artifactStore = createArtifactStore({ dataRoot, runsRoot, sourceRoot, pythonPath, bridgePath });
const stageStore = createStageStore({
  runsRoot,
  appendEvent: artifactStore.appendEvent,
  refreshManifest: artifactStore.refreshManifest,
});
const buildStageInputs = createStageInputBuilder({
  store: artifactStore,
  stageStore,
  sourceRoot,
  assetRoot: externalAssetRoot,
});
const workerManager = createWorkerManager({
  store: artifactStore,
  stageStore,
  buildStageInputs,
  sourceRoot,
  pythonPath,
  bridgePath,
  assetRoot: externalAssetRoot,
  sminaBin: process.env.SMINA_BIN || "",
  obabelBin: process.env.OBABEL_BIN || "",
});

async function createPrepareAttempt(run) {
  const runRoot = path.join(runsRoot, run.runId);
  const stageInputs = await buildStageInputs(run);
  const stageInput = stageInputs.prepare;
  if (!stageInput) throw new Error("prepare stage input is unavailable");
  const command = { executable: "OPEN_MOLECULE_SERVER", args: ["prepare"] };
  let attempt = null;
  try {
    attempt = await stageStore.createAttempt(
      run.runId,
      "prepare",
      computeStageFingerprint(stageInput),
      command,
      stageInput,
    );
    await atomicWriteJson(path.join(attempt.root, "command.json"), command);
    await stageStore.transition(run.runId, "prepare", attempt.attempt, "running");
    await fs.copyFile(path.join(runRoot, "run-spec.json"), path.join(attempt.root, "outputs", "run-spec.json"));
    await fs.copyFile(path.join(runRoot, "preflight.json"), path.join(attempt.root, "outputs", "preflight.json"));
    await fs.copyFile(
      path.join(runRoot, "inputs", "molecules.csv"),
      path.join(attempt.root, "outputs", "molecules.csv"),
    );
    await stageStore.complete(run.runId, "prepare", attempt.attempt, {
      "run-spec.json": "outputs/run-spec.json",
      "preflight.json": "outputs/preflight.json",
      "molecules.csv": "outputs/molecules.csv",
    });
  } catch (error) {
    if (attempt) {
      await stageStore.fail(
        run.runId,
        "prepare",
        attempt.attempt,
        { code: "prepare_failed", message: "Failed to seal prepare-stage evidence" },
      ).catch(() => {});
    }
    throw error;
  }
}

async function recoverInterruptedRuns() {
  for (const run of await artifactStore.listRuns()) {
    if (!["queued", "running"].includes(run.status)) continue;
    if (run.pid) await workerManager.terminateOrphan(run.pid, "SIGTERM");
    if (run.currentStage && run.currentAttempt) {
      const attempts = await stageStore.listAttempts(run.runId, run.currentStage);
      const current = attempts.find((attempt) => attempt.attempt === run.currentAttempt);
      if (current && ["queued", "running"].includes(current.status)) {
        await stageStore.fail(
          run.runId,
          run.currentStage,
          run.currentAttempt,
          { code: "worker_interrupted", message: "Server stopped before this stage attempt completed" },
        );
      }
    }
    await artifactStore.updateRun(
      run.runId,
      {
        status: "failed",
        finishedAt: new Date().toISOString(),
        pid: null,
        error: { code: "worker_interrupted", message: "server restarted before the worker reached a terminal state" },
      },
      { event: "worker_interrupted", status: "failed" },
    );
  }
}

await recoverInterruptedRuns();

const scientificBoundary =
  "plan_only: no molecule set is attached and no L1/L2/L3/L4 scoring or docking has been executed.";

function sha256(value) {
  return createHash("sha256").update(value).digest("hex");
}

function jsonResponse(response, status, payload) {
  const body = JSON.stringify(payload, null, 2);
  response.writeHead(status, {
    "Content-Type": "application/json; charset=utf-8",
    "Content-Length": Buffer.byteLength(body),
    "Cache-Control": "no-store",
  });
  response.end(body);
}

function errorPayload(code, message, field) {
  return {
    ok: false,
    error: {
      code,
      message,
      ...(field ? { field } : {}),
    },
  };
}

async function readJsonBody(request) {
  const chunks = [];
  let total = 0;
  for await (const chunk of request) {
    total += chunk.length;
    if (total > maxBodyBytes) {
      const error = new Error(`request body exceeds ${maxBodyBytes} bytes`);
      error.status = 413;
      error.code = "payload_too_large";
      error.field = "body";
      throw error;
    }
    chunks.push(chunk);
  }
  const text = Buffer.concat(chunks).toString("utf8");
  if (!text) return {};
  return JSON.parse(text);
}

function pythonModuleAvailable(moduleName) {
  const probe = spawnSync(
    pythonPath,
    ["-c", `import importlib.util,sys;sys.exit(0 if importlib.util.find_spec(${JSON.stringify(moduleName)}) else 1)`],
    { stdio: "ignore", timeout: 3000 },
  );
  return probe.status === 0;
}

function executablePath(names) {
  const explicit = names
    .map((name) => process.env[`${String(name).toUpperCase()}_BIN`])
    .filter(Boolean);
  for (const candidate of [...explicit, ...names]) {
    const value = String(candidate || "").trim();
    if (!value) continue;
    if (value.includes(path.sep) && existsSync(value)) return value;
    const probe = spawnSync("which", [value], { encoding: "utf8", timeout: 2000 });
    if (probe.status === 0 && String(probe.stdout || "").trim()) return String(probe.stdout).trim();
  }
  return "";
}

function modelStatusRows() {
  const scoringRoot = externalAssetRoot
    ? path.join(externalAssetRoot, "scoring")
    : path.join(sourceRoot, "scoring");
  const modelRoot = path.join(scoringRoot, "models");
  const l2Available = assetBundleVerified && [
    path.join(modelRoot, "bindingdb_l2", "l2_model_sklearn_1_7_2.joblib"),
    path.join(modelRoot, "bindingdb_l2", "l2_model.joblib"),
  ].some(existsSync);
  const l3Available = assetBundleVerified && ["tox21.pkl", "bbbp.pkl", "clintox.pkl", "sider.pkl"].every((name) =>
    existsSync(path.join(modelRoot, "admet", name)),
  );
  const l4Available = assetBundleVerified && [
    path.join(modelRoot, "ref_embeddings.npz"),
    path.join(modelRoot, "ref_smiles.pkl"),
  ].every(existsSync);
  const dockingAvailable = Boolean(executablePath(["smina"])) && Boolean(executablePath(["obabel"]));
  return [
    {
      id: "l1-rdkit",
      label: "RDKit molecular quality",
      layer: "L1",
      status: pythonModuleAvailable("rdkit") ? "available" : "not_found",
      role: "分子质量、描述符和结构合法性",
      requirement: "Python RDKit runtime",
    },
    {
      id: "l2-bindingdb",
      label: "BindingDB target binding",
      layer: "L2",
      status: l2Available ? "available" : "not_found",
      role: "靶点感知的结合概率粗筛",
      requirement: "BindingDB L2 model asset",
    },
    {
      id: "l3-admet",
      label: "ADMET safety panel",
      layer: "L3",
      status: l3Available ? "available" : "not_found",
      role: "ADMET 与毒性风险证据",
      requirement: "four ADMET model assets",
    },
    {
      id: "l4-unimol",
      label: "UniMol reference similarity",
      layer: "L4",
      status: l4Available ? "available" : "not_found",
      role: "三维表征与参考药物相似度",
      requirement: "UniMol weights and reference embeddings",
    },
    {
      id: "dock-smina",
      label: "smina docking cascade",
      layer: "Dock",
      status: dockingAvailable ? "available" : "not_found",
      role: "有受体时对 L2 头部候选进行精排",
      requirement: "registered receptor, smina and obabel",
    },
  ];
}

function normalizeInteger(value, field, min, max) {
  const parsed = Number(value);
  if (!Number.isInteger(parsed) || parsed < min || parsed > max) {
    const error = new Error(`${field} must be an integer from ${min} to ${max}`);
    error.field = field;
    throw error;
  }
  return parsed;
}

function safeSlug(value, fallback = "research_target") {
  const normalized = String(value || "")
    .normalize("NFKC")
    .toLowerCase()
    .replace(/[^\p{L}\p{N}]+/gu, "_")
    .replace(/^_+|_+$/g, "")
    .slice(0, 96);
  return normalized || fallback;
}

function inferTarget(prompt, explicitTarget) {
  const target = String(explicitTarget || "").trim();
  if (target) return target.slice(0, 160);
  const chembl = String(prompt || "").match(/CHEMBL\d+/i)?.[0];
  if (chembl) return chembl.toUpperCase();
  const text = String(prompt || "").toLowerCase();
  if (/neuraminidase|神经氨酸酶|influenza|流感/.test(text)) return "CHEMBL2051";
  if (/egfr/.test(text)) return "EGFR";
  return safeSlug(text);
}

function normalizeTargetName(value) {
  return String(value || "").trim().toLowerCase().replace(/[_\s-]+/g, "");
}

async function receptorForTarget(targetId) {
  try {
    const registry = JSON.parse(await fs.readFile(receptorRegistry, "utf8"));
    const entries = registry.entries || {};
    const exact = entries[targetId];
    const matched =
      exact ||
      Object.entries(entries).find(([key, entry]) => {
        const target = normalizeTargetName(targetId);
        return normalizeTargetName(key) === target || normalizeTargetName(entry.target_name) === target;
      })?.[1];
    if (!matched) return { registered: false, available: false, resourceId: "" };
    const sourcePath = path.isAbsolute(matched.pdbqt)
      ? matched.pdbqt
      : path.resolve(path.dirname(receptorRegistry), matched.pdbqt);
    const externalPath = externalAssetRoot && !path.isAbsolute(matched.pdbqt)
      ? path.resolve(externalAssetRoot, "scoring", matched.pdbqt)
      : "";
    const receptorPath = externalAssetRoot && !path.isAbsolute(matched.pdbqt)
      ? externalPath
      : sourcePath;
    return {
      registered: true,
      available: existsSync(receptorPath),
      resourceId: safeSlug(path.basename(receptorPath), "registered_receptor"),
    };
  } catch {
    return { registered: false, available: false, resourceId: "" };
  }
}

function resolveRoute(requestedRoute, receptor) {
  if (requestedRoute === "cascade") {
    if (receptor.available) {
      return {
        branch: "cascade",
        status: "ready",
        receptorAvailable: true,
        rationale: "Explicit cascade route accepted because the registered receptor asset is available.",
      };
    }
    return {
      branch: "cascade",
      status: "blocked",
      receptorAvailable: false,
      rationale: receptor.registered
        ? "Explicit cascade route is blocked because the registered receptor file is not present in this source-only workspace."
        : "Explicit cascade route is blocked because the target has no registered receptor.",
    };
  }
  if (requestedRoute === "library") {
    return {
      branch: "library",
      status: "ready",
      receptorAvailable: receptor.available,
      rationale: "Library route was explicitly requested; docking is not part of this plan.",
    };
  }
  if (receptor.available) {
    return {
      branch: "cascade",
      status: "ready",
      receptorAvailable: true,
      rationale: "Auto route selected cascade because a registered receptor asset is available.",
    };
  }
  return {
    branch: "library",
    status: "ready",
    receptorAvailable: false,
    rationale: "Auto route selected the library branch because no usable receptor asset is available.",
  };
}

function plannedStages(route) {
  return [
    { id: "prepare", label: "准备输入", status: "planned", reason: "Attach and validate a molecule-set CSV." },
    { id: "score", label: "四级评分", status: "planned", reason: "Run strict L1/L2/L3/L4 scoring." },
    route.branch === "cascade"
      ? {
          id: "dock",
          label: "Docking cascade",
          status: route.status === "blocked" ? "blocked" : "planned",
          reason: route.rationale,
        }
      : { id: "dock", label: "Docking cascade", status: "skipped", reason: "Library branch does not run docking." },
    { id: "report", label: "证据报告", status: "planned", reason: "Write summary, provenance and manifest." },
  ];
}

function equivalentCli(targetId, runId, route) {
  const command = [
    "four-level-molecule \\",
    "  --input <molecule-set.csv> \\",
    `  --target ${JSON.stringify(targetId)} \\`,
    `  --mode ${route.branch} \\`,
    "  --asset-root '<asset-root>' \\",
    "  --strict-backends \\",
    `  --output runs/${runId}/scores.csv`,
  ];
  if (route.branch === "cascade") {
    command.splice(command.length - 1, 0, "  --cascade-top-n 300 \\");
  }
  return command.join("\n");
}

async function writePlanBundle({ runId, prompt, spec, route, stages, assetRequirements, createdAt }) {
  const runRoot = path.join(runsRoot, runId);
  await fs.mkdir(runRoot, { recursive: false });
  const payload = {
    schemaVersion: "open-molecule-lab.plan-bundle.v0.1",
    runId,
    createdAt,
    spec,
    route,
    stages,
    assetRequirements,
    executionBoundary: scientificBoundary,
  };
  const runJson = `${JSON.stringify(payload, null, 2)}\n`;
  const design = [
    "# Open Molecule Lab Plan",
    "",
    `- Run ID: ${runId}`,
    `- Target: ${spec.target.id}`,
    `- Requested route: ${spec.target.requestedRoute}`,
    `- Resolved branch: ${route.branch}`,
    `- Route status: ${route.status}`,
    `- Candidate count: ${spec.moleculeSet.expectedCandidateCount}`,
    `- Final selection count: ${spec.selection.finalSelectionCount}`,
    "",
    "## Boundary",
    "",
    scientificBoundary,
    "",
  ].join("\n");
  const files = {
    "run.json": runJson,
    "prompt.txt": `${prompt.trim()}\n`,
    "DESIGN.md": design,
  };
  for (const [name, content] of Object.entries(files)) {
    await fs.writeFile(path.join(runRoot, name), content, "utf8");
  }
  const manifest = Object.entries(files)
    .sort(([left], [right]) => left.localeCompare(right))
    .map(([name, content]) => `${sha256(content)}  ${name}`)
    .join("\n") + "\n";
  await fs.writeFile(path.join(runRoot, "MANIFEST.sha256"), manifest, "utf8");
  return {
    relativeRoot: path.posix.join("runs", runId),
    files: [...Object.keys(files), "MANIFEST.sha256"].sort(),
    manifestSha256: sha256(manifest),
  };
}

async function createPromptPlan(body) {
  if (!body || typeof body !== "object" || Array.isArray(body)) {
    const error = new Error("request body must be a JSON object");
    error.field = "body";
    throw error;
  }
  const prompt = String(body.prompt || "").trim();
  if (!prompt) {
    const error = new Error("prompt is required");
    error.field = "prompt";
    throw error;
  }
  if (prompt.length > 8000) {
    const error = new Error("prompt must not exceed 8000 characters");
    error.field = "prompt";
    throw error;
  }
  const candidatePool = normalizeInteger(body.candidatePool, "candidatePool", 1, 100000);
  const finalSelectionCount = normalizeInteger(body.finalSelectionCount, "finalSelectionCount", 1, 100);
  if (finalSelectionCount > candidatePool) {
    const error = new Error("finalSelectionCount must not exceed candidatePool");
    error.field = "finalSelectionCount";
    throw error;
  }
  const routeMode = String(body.routeMode || "auto");
  if (!new Set(["auto", "library", "cascade"]).has(routeMode)) {
    const error = new Error("routeMode must be auto, library, or cascade");
    error.field = "routeMode";
    throw error;
  }
  const targetId = inferTarget(prompt, body.target);
  const receptor = await receptorForTarget(targetId);
  const route = resolveRoute(routeMode, receptor);
  const createdAt = new Date().toISOString();
  const runId = `plan_${createdAt.replace(/[-:.TZ]/g, "").slice(0, 14)}_${randomBytes(4).toString("hex")}`;
  const spec = {
    schemaVersion: "open-molecule-lab.run-spec.v0.1",
    mode: "plan_only",
    project: {
      name: safeSlug(targetId, "molecule_project"),
      researchPrompt: prompt,
    },
    target: {
      id: targetId,
      text: targetId,
      requestedRoute: routeMode,
      resolvedBranch: route.branch,
    },
    moleculeSet: {
      attached: false,
      expectedCandidateCount: candidatePool,
    },
    selection: {
      finalSelectionCount,
    },
    execution: {
      strictBackends: true,
      worker: "local",
      seed: 42,
    },
  };
  const stages = plannedStages(route);
  const assetRequirements = modelStatusRows()
    .filter((row) => row.status !== "available")
    .map((row) => row.requirement);
  if (!assetRequirements.includes("validated molecule-set CSV")) {
    assetRequirements.unshift("validated molecule-set CSV");
  }
  await fs.mkdir(runsRoot, { recursive: true });
  const bundle = await writePlanBundle({ runId, prompt, spec, route, stages, assetRequirements, createdAt });
  return {
    ok: true,
    mode: "plan_only",
    runId,
    spec,
    route,
    stages,
    assetRequirements,
    executionBoundary: scientificBoundary,
    equivalentCli: equivalentCli(targetId, runId, route),
    bundle,
  };
}

async function verifyResumeEnvelope(run) {
  const runRoot = path.join(runsRoot, run.runId);
  const integrity = await artifactStore.verifyManifest(runRoot);
  if (!integrity.ok) {
    return { ok: false, code: "checkpoint_mismatch", resource: "MANIFEST.sha256" };
  }
  let spec;
  let specBytes;
  let inputBytes;
  try {
    specBytes = await fs.readFile(path.join(runRoot, "run-spec.json"));
    spec = JSON.parse(specBytes.toString("utf8"));
    inputBytes = await fs.readFile(path.join(runRoot, "inputs", "molecules.csv"));
  } catch {
    return { ok: false, code: "checkpoint_invalid", resource: "run-spec.json" };
  }
  if (sha256(specBytes) !== run.specSha256) {
    return { ok: false, code: "checkpoint_mismatch", resource: "run-spec.json" };
  }
  const moleculeSet = await artifactStore.getMoleculeSet(run.moleculeSetId);
  if (!moleculeSet) {
    return { ok: false, code: "checkpoint_invalid", resource: "molecule_set" };
  }
  const inputSha256 = sha256(inputBytes);
  if (
    inputSha256 !== moleculeSet.inputSha256
    || inputSha256 !== spec.moleculeSet?.inputSha256
    || moleculeSet.nRows !== spec.moleculeSet?.nCandidates
  ) {
    return { ok: false, code: "checkpoint_mismatch", resource: "inputs/molecules.csv" };
  }
  return { ok: true, spec, moleculeSet };
}

async function resumeRun(run) {
  const envelope = await verifyResumeEnvelope(run);
  if (!envelope.ok) {
    return artifactStore.updateRun(
      run.runId,
      {
        status: "blocked",
        finishedAt: new Date().toISOString(),
        error: {
          code: envelope.code,
          message: "Verified resume envelope mismatch",
          resource: envelope.resource,
        },
      },
      { event: "resume_blocked", status: "blocked", resource: envelope.resource },
    );
  }
  await artifactStore.appendEvent(run.runId, { event: "resume_requested", status: run.status });
  await artifactStore.refreshManifest(path.join(runsRoot, run.runId));
  const preflight = await runPreflight({
    sourceRoot,
    pythonPath,
    assetRoot: externalAssetRoot,
    routeBranch: run.route.branch,
    targetId: envelope.spec.target.id,
    expectedCandidateCount: envelope.spec.moleculeSet.nCandidates,
    actualCandidateCount: envelope.moleculeSet.nRows,
    sminaBin: process.env.SMINA_BIN || "",
    obabelBin: process.env.OBABEL_BIN || "",
  });
  if (!preflight.ok) {
    return artifactStore.updateRun(
      run.runId,
      {
        status: "blocked",
        finishedAt: new Date().toISOString(),
        resumePreflight: preflight,
        error: { code: "resume_preflight_failed", message: "Strict resume preflight failed" },
      },
      { event: "resume_blocked", status: "blocked", reason: "resume_preflight_failed" },
    );
  }
  let currentInputs;
  try {
    currentInputs = await buildStageInputs(run);
  } catch {
    return artifactStore.updateRun(
      run.runId,
      {
        status: "blocked",
        finishedAt: new Date().toISOString(),
        error: { code: "checkpoint_mismatch", message: "Runtime code or asset identity is unavailable" },
      },
      { event: "resume_blocked", status: "blocked", reason: "runtime_identity_unavailable" },
    );
  }
  const verification = await stageStore.verifyContiguous(run, currentInputs);
  if (!verification.ok) {
    return artifactStore.updateRun(
      run.runId,
      {
        status: "blocked",
        finishedAt: new Date().toISOString(),
        error: {
          code: verification.code,
          message: "Verified resume checkpoint mismatch",
          stage: verification.stage,
          resource: verification.resource,
        },
      },
      {
        event: "resume_blocked",
        status: "blocked",
        stage: verification.stage,
        resource: verification.resource,
      },
    );
  }
  if (!verification.nextStage) {
    return { ...run, resumeError: "run_already_complete" };
  }
  const queued = await artifactStore.updateRun(
    run.runId,
    {
      status: "queued",
      finishedAt: null,
      pid: null,
      currentStage: verification.nextStage,
      currentAttempt: null,
      error: null,
      resumePreflight: preflight,
      resumeCount: Number(run.resumeCount || 0) + 1,
    },
    { event: "resume_started", status: "queued", stage: verification.nextStage },
  );
  return workerManager.start(queued, { startStage: verification.nextStage });
}

async function handleApi(request, response, url) {
  if (request.method === "GET" && url.pathname === "/api/health") {
    return jsonResponse(response, 200, {
      ok: true,
      product: "open-molecule-lab",
      mode: "local_execution",
      sourceRoot: ".",
      cli: {
        available: existsSync(scoringEntry) && existsSync(benchmarkEntry),
        scoringEntry: "scoring/scoring.py",
        benchmarkEntry: "scientific_validation/four_level_cli_1kx10k/batch_cli.py",
      },
      assetManifestPresent: existsSync(
        externalAssetRoot ? path.join(externalAssetRoot, "ASSET_MANIFEST.json") : assetManifest,
      ),
    });
  }
  if (request.method === "GET" && url.pathname === "/api/model-status") {
    return jsonResponse(response, 200, { ok: true, models: modelStatusRows() });
  }
  if (request.method === "POST" && url.pathname === "/api/molecule-sets") {
    try {
      const result = await artifactStore.createMoleculeSet(await readJsonBody(request));
      return jsonResponse(response, result.created ? 201 : 200, { ok: true, ...result.metadata });
    } catch (error) {
      const known = error instanceof StoreError;
      return jsonResponse(
        response,
        known ? error.status : Number(error?.status) || 400,
        errorPayload(
          known ? error.code : error?.code || "invalid_molecule_set",
          error instanceof Error ? error.message : String(error),
          error?.field,
        ),
      );
    }
  }
  if (request.method === "GET" && url.pathname.startsWith("/api/molecule-sets/")) {
    const moleculeSetId = decodeURIComponent(url.pathname.slice("/api/molecule-sets/".length));
    const metadata = await artifactStore.getMoleculeSet(moleculeSetId);
    return metadata
      ? jsonResponse(response, 200, { ok: true, ...metadata })
      : jsonResponse(response, 404, errorPayload("molecule_set_not_found", "Molecule set not found"));
  }
  if (request.method === "POST" && url.pathname === "/api/runs") {
    try {
      const body = await readJsonBody(request);
      const plan = await artifactStore.getPlan(body.planRunId);
      if (!plan) throw new StoreError("Plan run not found", { code: "plan_not_found", field: "planRunId", status: 404 });
      const moleculeSet = await artifactStore.getMoleculeSet(body.moleculeSetId);
      if (!moleculeSet) throw new StoreError("Molecule set not found", { code: "molecule_set_not_found", field: "moleculeSetId", status: 404 });
      const preflight = await runPreflight({
        sourceRoot,
        pythonPath,
        assetRoot: externalAssetRoot,
        routeBranch: plan.route.branch,
        targetId: plan.spec.target.id,
        expectedCandidateCount: plan.spec.moleculeSet.expectedCandidateCount,
        actualCandidateCount: moleculeSet.nRows,
        sminaBin: process.env.SMINA_BIN || "",
        obabelBin: process.env.OBABEL_BIN || "",
      });
      let run = await artifactStore.createExecutionRun({
        planRunId: body.planRunId,
        plan,
        moleculeSet,
        preflight,
      });
      if (run.status === "queued") {
        try {
          await createPrepareAttempt(run);
          await workerManager.start(run);
        } catch {
          await artifactStore.updateRun(
            run.runId,
            {
              status: "failed",
              finishedAt: new Date().toISOString(),
              error: { code: "prepare_failed", message: "Failed to create verified prepare stage" },
            },
            { event: "worker_failed", status: "failed", stage: "prepare" },
          );
        }
        run = await artifactStore.getRun(run.runId);
      }
      return jsonResponse(response, 201, { ok: true, ...run });
    } catch (error) {
      const known = error instanceof StoreError;
      return jsonResponse(
        response,
        known ? error.status : Number(error?.status) || 400,
        errorPayload(
          known ? error.code : error?.code || "invalid_run",
          error instanceof Error ? error.message : String(error),
          error?.field,
        ),
      );
    }
  }
  if (request.method === "GET" && url.pathname.startsWith("/api/runs/") && url.pathname.endsWith("/results")) {
    const runId = decodeURIComponent(url.pathname.slice("/api/runs/".length, -"/results".length));
    const run = await artifactStore.getRun(runId);
    if (!run) return jsonResponse(response, 404, errorPayload("run_not_found", "Run not found"));
    if (run.status !== "complete") {
      return jsonResponse(response, 409, errorPayload("run_not_complete", `Run status is ${run.status}`));
    }
    try {
      const integrity = await artifactStore.verifyManifest(path.join(runsRoot, runId));
      if (!integrity.ok) {
        return jsonResponse(
          response,
          409,
          errorPayload("evidence_integrity_failed", "Run evidence does not match MANIFEST.sha256"),
        );
      }
      const summary = JSON.parse(await fs.readFile(path.join(runsRoot, runId, "results", "summary.json"), "utf8"));
      const rawOffset = Number(url.searchParams.get("offset") || 0);
      const rawLimit = Number(url.searchParams.get("limit") || 100);
      const view = url.searchParams.get("view") || "ranked";
      if (!Number.isInteger(rawOffset) || rawOffset < 0 || !Number.isInteger(rawLimit) || rawLimit < 1 || rawLimit > 200 || !new Set(["ranked", "failed", "all"]).has(view)) {
        return jsonResponse(response, 400, errorPayload("invalid_pagination", "offset must be >= 0, limit must be 1..200, view must be ranked/failed/all", "pagination"));
      }
      const pageProcess = spawnSync(
        pythonPath,
        [bridgePath, "page-results", "--input", path.join(runsRoot, runId, "results", "scores.csv"), "--offset", String(rawOffset), "--limit", String(rawLimit), "--view", view],
        { cwd: sourceRoot, encoding: "utf8", timeout: 30_000 },
      );
      const pageLine = String(pageProcess.stdout || "").trim().split("\n").filter(Boolean).at(-1);
      const page = pageLine ? JSON.parse(pageLine) : null;
      if (pageProcess.status !== 0 || !page?.ok) {
        return jsonResponse(response, 500, errorPayload("results_page_failed", page?.error || "result page failed"));
      }
      return jsonResponse(response, 200, { ...summary, ...page });
    } catch {
      return jsonResponse(response, 500, errorPayload("results_missing", "Completed run has no result summary"));
    }
  }
  if (request.method === "GET" && url.pathname.startsWith("/api/runs/") && url.pathname.endsWith("/stages")) {
    const runId = decodeURIComponent(url.pathname.slice("/api/runs/".length, -"/stages".length));
    const run = await artifactStore.getRun(runId);
    if (!run) return jsonResponse(response, 404, errorPayload("run_not_found", "Run not found"));
    return jsonResponse(response, 200, await stageStore.summarize(run));
  }
  if (request.method === "POST" && url.pathname.startsWith("/api/runs/") && url.pathname.endsWith("/resume")) {
    const runId = decodeURIComponent(url.pathname.slice("/api/runs/".length, -"/resume".length));
    const run = await artifactStore.getRun(runId);
    if (!run) return jsonResponse(response, 404, errorPayload("run_not_found", "Run not found"));
    if (!["failed", "cancelled"].includes(run.status)) {
      return jsonResponse(
        response,
        409,
        errorPayload("run_not_resumable", `Run status is ${run.status}`),
      );
    }
    const resumed = await resumeRun(run);
    if (resumed.resumeError === "run_already_complete") {
      return jsonResponse(response, 409, errorPayload("run_already_complete", "No stage remains"));
    }
    return jsonResponse(response, resumed.status === "blocked" ? 409 : 202, { ok: true, ...resumed });
  }
  if (request.method === "POST" && url.pathname.startsWith("/api/runs/") && url.pathname.endsWith("/cancel")) {
    const runId = decodeURIComponent(url.pathname.slice("/api/runs/".length, -"/cancel".length));
    const run = await artifactStore.getRun(runId);
    if (!run) return jsonResponse(response, 404, errorPayload("run_not_found", "Run not found"));
    if (!["queued", "running"].includes(run.status)) {
      return jsonResponse(response, 409, errorPayload("run_not_cancellable", `Run status is ${run.status}`));
    }
    const cancelled = await workerManager.cancel(runId);
    return jsonResponse(response, 202, { ok: true, ...cancelled });
  }
  if (request.method === "GET" && url.pathname.startsWith("/api/runs/")) {
    const runId = decodeURIComponent(url.pathname.slice("/api/runs/".length));
    const run = await artifactStore.getRun(runId);
    return run
      ? jsonResponse(response, 200, { ok: true, ...run })
      : jsonResponse(response, 404, errorPayload("run_not_found", "Run not found"));
  }
  if (request.method === "POST" && url.pathname === "/api/prompt-plan") {
    try {
      const body = await readJsonBody(request);
      return jsonResponse(response, 201, await createPromptPlan(body));
    } catch (error) {
      return jsonResponse(
        response,
        Number(error?.status) || 400,
        errorPayload(error?.code || "invalid_prompt_plan", error instanceof Error ? error.message : String(error), error?.field),
      );
    }
  }
  return jsonResponse(response, 404, errorPayload("not_found", "API endpoint not found"));
}

const mimeTypes = {
  ".css": "text/css; charset=utf-8",
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript; charset=utf-8",
  ".json": "application/json; charset=utf-8",
  ".svg": "image/svg+xml",
};

async function serveStatic(request, response, url) {
  const relative = url.pathname === "/" ? "index.html" : decodeURIComponent(url.pathname.slice(1));
  let filePath = path.resolve(distRoot, relative);
  if (!filePath.startsWith(`${distRoot}${path.sep}`) && filePath !== path.join(distRoot, "index.html")) {
    response.writeHead(403);
    return response.end("Forbidden");
  }
  try {
    const stat = await fs.stat(filePath);
    if (!stat.isFile()) throw new Error("not a file");
  } catch {
    filePath = path.join(distRoot, "index.html");
  }
  try {
    const content = await fs.readFile(filePath);
    response.writeHead(200, {
      "Content-Type": mimeTypes[path.extname(filePath)] || "application/octet-stream",
      "Content-Length": content.length,
    });
    response.end(content);
  } catch {
    response.writeHead(503, { "Content-Type": "text/plain; charset=utf-8" });
    response.end("Frontend is not built. Run npm run build first.\n");
  }
}

const server = createServer(async (request, response) => {
  try {
    const url = new URL(request.url || "/", `http://${request.headers.host || `${host}:${port}`}`);
    if (url.pathname.startsWith("/api/")) return await handleApi(request, response, url);
    return await serveStatic(request, response, url);
  } catch (error) {
    return jsonResponse(response, 500, errorPayload("internal_error", error instanceof Error ? error.message : String(error)));
  }
});

server.listen(port, host, () => {
  const address = server.address();
  const actualPort = typeof address === "object" && address ? address.port : port;
  process.stdout.write(`Open Molecule Lab listening on http://${host}:${actualPort}\n`);
});

let shuttingDown = false;
for (const signal of ["SIGINT", "SIGTERM"]) {
  process.on(signal, async () => {
    if (shuttingDown) return;
    shuttingDown = true;
    await workerManager.shutdown();
    server.close(() => process.exit(0));
  });
}
