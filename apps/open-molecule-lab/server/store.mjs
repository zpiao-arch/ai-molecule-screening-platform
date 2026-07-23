import { createHash, randomBytes } from "node:crypto";
import { spawnSync } from "node:child_process";
import fs from "node:fs/promises";
import path from "node:path";


export class StoreError extends Error {
  constructor(message, { code = "invalid_molecule_set", field, status = 400 } = {}) {
    super(message);
    this.code = code;
    this.field = field;
    this.status = status;
  }
}

export function sha256(value) {
  return createHash("sha256").update(value).digest("hex");
}

function manifestText(files) {
  return Object.entries(files)
    .sort(([left], [right]) => left.localeCompare(right))
    .map(([name, content]) => `${sha256(content)}  ${name}`)
    .join("\n") + "\n";
}

async function collectFiles(root, relative = "") {
  const directory = path.join(root, relative);
  const entries = await fs.readdir(directory, { withFileTypes: true });
  const files = {};
  for (const entry of entries) {
    const child = path.posix.join(relative.split(path.sep).join(path.posix.sep), entry.name);
    if (child === "MANIFEST.sha256") continue;
    if (entry.isDirectory()) {
      Object.assign(files, await collectFiles(root, child));
    } else if (entry.isFile()) {
      files[child] = await fs.readFile(path.join(root, child));
    }
  }
  return files;
}

async function refreshManifest(root) {
  await fs.writeFile(path.join(root, "MANIFEST.sha256"), manifestText(await collectFiles(root)), "utf8");
}

async function verifyManifest(root) {
  let manifest;
  try {
    manifest = await fs.readFile(path.join(root, "MANIFEST.sha256"), "utf8");
  } catch {
    return { ok: false, missing: ["MANIFEST.sha256"], mismatches: [], unexpected: [], malformed: [] };
  }
  const expected = new Map();
  const malformed = [];
  for (const line of manifest.split("\n").filter(Boolean)) {
    const match = line.match(/^([0-9a-f]{64})  (.+)$/);
    if (!match || path.isAbsolute(match?.[2] || "") || String(match?.[2] || "").split("/").includes("..")) {
      malformed.push(line);
      continue;
    }
    expected.set(match[2], match[1]);
  }
  const actualFiles = await collectFiles(root);
  const actualNames = new Set(Object.keys(actualFiles));
  const missing = [...expected.keys()].filter((name) => !actualNames.has(name)).sort();
  const unexpected = [...actualNames].filter((name) => !expected.has(name)).sort();
  const mismatches = [];
  for (const [name, digest] of expected.entries()) {
    if (!actualNames.has(name)) continue;
    if (sha256(actualFiles[name]) !== digest) mismatches.push(name);
  }
  mismatches.sort();
  return {
    ok: !missing.length && !unexpected.length && !mismatches.length && !malformed.length,
    missing,
    mismatches,
    unexpected,
    malformed,
  };
}

function safeLabel(value, fallback) {
  const normalized = String(value || "").normalize("NFKC").trim().slice(0, 160);
  return normalized || fallback;
}

function parseBridgeOutput(result) {
  const output = String(result.stdout || "").trim().split("\n").filter(Boolean).at(-1);
  if (!output) {
    throw new StoreError(
      "Python CSV validator did not start; configure OPEN_MOLECULE_PYTHON to a Python 3.11 runtime with the CLI dependencies",
      { code: "python_runtime_unavailable", status: 503 },
    );
  }
  try {
    return JSON.parse(output || "{}");
  } catch {
    throw new StoreError(`molecule validator returned invalid JSON: ${String(result.stderr || "").trim()}`, {
      code: "validator_failed",
      status: 500,
    });
  }
}

export function createArtifactStore({ dataRoot, runsRoot, sourceRoot, pythonPath, bridgePath }) {
  const moleculeSetsRoot = path.join(dataRoot, "molecule_sets");

  async function getMoleculeSet(moleculeSetId) {
    if (!/^molset_[0-9a-f]{16}$/.test(String(moleculeSetId || ""))) return null;
    try {
      const metadata = JSON.parse(
        await fs.readFile(path.join(moleculeSetsRoot, moleculeSetId, "metadata.json"), "utf8"),
      );
      return metadata;
    } catch {
      return null;
    }
  }

  async function getPlan(planRunId) {
    if (!/^plan_[0-9a-f_]+$/.test(String(planRunId || ""))) return null;
    try {
      return JSON.parse(await fs.readFile(path.join(runsRoot, planRunId, "run.json"), "utf8"));
    } catch {
      return null;
    }
  }

  async function getRun(runId) {
    if (!/^run_[0-9a-f_]+$/.test(String(runId || ""))) return null;
    try {
      return JSON.parse(await fs.readFile(path.join(runsRoot, runId, "run.json"), "utf8"));
    } catch {
      return null;
    }
  }

  async function listRuns() {
    const entries = await fs.readdir(runsRoot, { withFileTypes: true }).catch(() => []);
    const runs = [];
    for (const entry of entries) {
      if (!entry.isDirectory() || !entry.name.startsWith("run_")) continue;
      const run = await getRun(entry.name);
      if (run) runs.push(run);
    }
    return runs;
  }

  async function appendEvent(runId, event) {
    await fs.appendFile(
      path.join(runsRoot, runId, "events.jsonl"),
      `${JSON.stringify({ runId, timestampUtc: new Date().toISOString(), ...event })}\n`,
      "utf8",
    );
  }

  async function updateRun(runId, patch, event) {
    const current = await getRun(runId);
    if (!current) throw new StoreError("Run not found", { code: "run_not_found", status: 404 });
    const next = { ...current, ...patch };
    await fs.writeFile(path.join(runsRoot, runId, "run.json"), `${JSON.stringify(next, null, 2)}\n`, "utf8");
    if (event) await appendEvent(runId, event);
    await refreshManifest(path.join(runsRoot, runId));
    return next;
  }

  async function createMoleculeSet({ name, csvText, license }) {
    if (typeof csvText !== "string" || !csvText.trim()) {
      throw new StoreError("csvText is required", { field: "csvText" });
    }
    const input = Buffer.from(csvText, "utf8");
    const inputSha256 = sha256(input);
    const moleculeSetId = `molset_${inputSha256.slice(0, 16)}`;
    const existing = await getMoleculeSet(moleculeSetId);
    if (existing) return { created: false, metadata: existing };

    await fs.mkdir(moleculeSetsRoot, { recursive: true });
    const stagingRoot = path.join(moleculeSetsRoot, `.staging-${process.pid}-${randomBytes(4).toString("hex")}`);
    const destination = path.join(moleculeSetsRoot, moleculeSetId);
    await fs.mkdir(stagingRoot, { recursive: false });
    try {
      const inputPath = path.join(stagingRoot, "input.csv");
      await fs.writeFile(inputPath, input);
      const validationResult = spawnSync(
        pythonPath,
        [bridgePath, "validate-molecule-set", "--input", inputPath, "--max-rows", "100000"],
        { cwd: sourceRoot, encoding: "utf8", timeout: 30_000 },
      );
      const validation = parseBridgeOutput(validationResult);
      if (validationResult.status !== 0 || !validation.ok) {
        throw new StoreError(validation.error || "molecule set validation failed", {
          field: validation.field || "csvText",
        });
      }
      const metadata = {
        schemaVersion: "open-molecule-lab.molecule-set.v0.1",
        moleculeSetId,
        name: safeLabel(name, moleculeSetId),
        license: safeLabel(license, "unspecified"),
        inputSha256,
        nRows: validation.nRows,
        columns: validation.columns,
        createdAt: new Date().toISOString(),
      };
      const metadataText = `${JSON.stringify(metadata, null, 2)}\n`;
      await fs.writeFile(path.join(stagingRoot, "metadata.json"), metadataText, "utf8");
      await fs.writeFile(
        path.join(stagingRoot, "MANIFEST.sha256"),
        manifestText({ "input.csv": input, "metadata.json": metadataText }),
        "utf8",
      );
      try {
        await fs.rename(stagingRoot, destination);
      } catch (error) {
        if (error?.code !== "EEXIST" && error?.code !== "ENOTEMPTY") throw error;
        const concurrent = await getMoleculeSet(moleculeSetId);
        if (!concurrent) throw error;
        return { created: false, metadata: concurrent };
      }
      return { created: true, metadata };
    } finally {
      await fs.rm(stagingRoot, { recursive: true, force: true });
    }
  }

  async function createExecutionRun({ planRunId, plan, moleculeSet, preflight }) {
    const createdAt = new Date().toISOString();
    const runId = `run_${createdAt.replace(/[-:.TZ]/g, "").slice(0, 14)}_${randomBytes(4).toString("hex")}`;
    const runRoot = path.join(runsRoot, runId);
    const status = preflight.ok ? "queued" : "blocked";
    const spec = {
      ...plan.spec,
      mode: "execute",
      moleculeSet: {
        attached: true,
        id: moleculeSet.moleculeSetId,
        inputSha256: moleculeSet.inputSha256,
        nCandidates: moleculeSet.nRows,
      },
    };
    const specText = `${JSON.stringify(spec, null, 2)}\n`;
    const record = {
      schemaVersion: "open-molecule-lab.run.v0.1",
      runId,
      planRunId,
      moleculeSetId: moleculeSet.moleculeSetId,
      specSha256: sha256(specText),
      status,
      route: plan.route,
      preflight,
      createdAt,
      startedAt: null,
      finishedAt: status === "blocked" ? createdAt : null,
      error: status === "blocked" ? { code: "preflight_blocked", message: "Strict preflight failed" } : null,
    };
    const events = [
      { event: "run_created", runId, status: "queued", timestampUtc: createdAt },
      {
        event: preflight.ok ? "preflight_complete" : "preflight_blocked",
        runId,
        status,
        timestampUtc: createdAt,
      },
    ];
    await fs.mkdir(path.join(runRoot, "inputs"), { recursive: true });
    await fs.mkdir(path.join(runRoot, "results"), { recursive: true });
    await fs.mkdir(path.join(runRoot, "logs"), { recursive: true });
    await fs.copyFile(
      path.join(moleculeSetsRoot, moleculeSet.moleculeSetId, "input.csv"),
      path.join(runRoot, "inputs", "molecules.csv"),
    );
    await fs.writeFile(path.join(runRoot, "run-spec.json"), specText, "utf8");
    await fs.writeFile(path.join(runRoot, "preflight.json"), `${JSON.stringify(preflight, null, 2)}\n`, "utf8");
    await fs.writeFile(path.join(runRoot, "run.json"), `${JSON.stringify(record, null, 2)}\n`, "utf8");
    await fs.writeFile(
      path.join(runRoot, "events.jsonl"),
      `${events.map((event) => JSON.stringify(event)).join("\n")}\n`,
      "utf8",
    );
    await refreshManifest(runRoot);
    return record;
  }

  return {
    dataRoot,
    runsRoot,
    moleculeSetsRoot,
    createMoleculeSet,
    getMoleculeSet,
    getPlan,
    getRun,
    listRuns,
    updateRun,
    appendEvent,
    createExecutionRun,
    refreshManifest,
    verifyManifest,
  };
}
