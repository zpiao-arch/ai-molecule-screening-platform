import { createHash } from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";

import { computeStageFingerprint, stageSequence, STAGE_SCHEMA } from "./stage-contract.mjs";
import { atomicWriteJson } from "./store.mjs";


const STAGES = ["prepare", "score", "dock", "report"];
const TERMINAL_STATUSES = new Set(["complete", "failed", "cancelled", "blocked"]);
const FAILURE_STATUSES = new Set(["failed", "cancelled", "blocked"]);
const ATTEMPT_PATTERN = /^attempt-(\d{4})$/;
const SHA256_PATTERN = /^[0-9a-f]{64}$/;

function sha256(value) {
  return createHash("sha256").update(value).digest("hex");
}

function assertIdentifier(value, label) {
  if (!/^[A-Za-z0-9_-]+$/.test(String(value || ""))) {
    throw new TypeError(`invalid ${label}`);
  }
  return String(value);
}

function assertStage(stage) {
  if (!STAGES.includes(stage)) throw new TypeError(`unsupported stage: ${stage}`);
  return stage;
}

function attemptName(attempt) {
  if (!Number.isInteger(attempt) || attempt < 1 || attempt > 9999) {
    throw new TypeError("invalid stage attempt number");
  }
  return `attempt-${String(attempt).padStart(4, "0")}`;
}

function attemptRoot(runsRoot, runId, stage, attempt) {
  return path.join(runsRoot, runId, "stages", stage, attemptName(attempt));
}

function cloneJson(value) {
  return JSON.parse(JSON.stringify(value));
}

function containsAbsolutePath(value) {
  if (typeof value === "string") return path.isAbsolute(value);
  if (Array.isArray(value)) return value.some(containsAbsolutePath);
  if (value && typeof value === "object") return Object.values(value).some(containsAbsolutePath);
  return false;
}

function safeOutputPath(value) {
  const candidate = String(value || "").replaceAll("\\", "/");
  const normalized = path.posix.normalize(candidate);
  if (
    !candidate
    || candidate !== normalized
    || path.posix.isAbsolute(candidate)
    || normalized === ".."
    || normalized.startsWith("../")
    || !normalized.startsWith("outputs/")
  ) {
    throw new TypeError(`invalid stage output path: ${candidate}`);
  }
  return normalized;
}

function checkpointResource(output) {
  return output?.path || "checkpoint.json";
}

async function readCheckpoint(checkpointPath) {
  const bytes = await fs.readFile(checkpointPath);
  return { bytes, checkpoint: JSON.parse(bytes.toString("utf8")) };
}

function validCheckpoint(checkpoint, runId, stage, attempt) {
  return Boolean(
    checkpoint
      && checkpoint.schemaVersion === STAGE_SCHEMA
      && checkpoint.runId === runId
      && checkpoint.stage === stage
      && checkpoint.attempt === attempt
      && ["queued", "running", ...TERMINAL_STATUSES].includes(checkpoint.status)
      && SHA256_PATTERN.test(String(checkpoint.inputFingerprint || ""))
      && checkpoint.outputs
      && typeof checkpoint.outputs === "object"
      && !Array.isArray(checkpoint.outputs),
  );
}

export function createStageStore({
  runsRoot,
  refreshManifest = async () => {},
  appendEvent = async () => {},
  now = () => new Date().toISOString(),
}) {
  if (!path.isAbsolute(runsRoot)) throw new TypeError("runsRoot must be absolute");

  function runRoot(runId) {
    return path.join(runsRoot, assertIdentifier(runId, "run id"));
  }

  async function persist(runId, stage, attempt, checkpoint, eventName) {
    const root = attemptRoot(runsRoot, runId, stage, attempt);
    await atomicWriteJson(path.join(root, "checkpoint.json"), checkpoint);
    await appendEvent(runId, {
      event: eventName,
      stage,
      attempt,
      status: checkpoint.status,
    });
    await refreshManifest(runRoot(runId));
    return cloneJson(checkpoint);
  }

  async function numberedAttemptEntries(runId, stage) {
    const stageRoot = path.join(runRoot(runId), "stages", assertStage(stage));
    const entries = await fs.readdir(stageRoot, { withFileTypes: true }).catch((error) => {
      if (error?.code === "ENOENT") return [];
      throw error;
    });
    return entries
      .filter((entry) => entry.isDirectory() && ATTEMPT_PATTERN.test(entry.name))
      .map((entry) => ({ name: entry.name, attempt: Number(ATTEMPT_PATTERN.exec(entry.name)[1]) }))
      .sort((left, right) => left.attempt - right.attempt);
  }

  async function createAttempt(runId, stage, inputFingerprint, command, inputs = {}) {
    assertIdentifier(runId, "run id");
    assertStage(stage);
    if (!SHA256_PATTERN.test(String(inputFingerprint || ""))) {
      throw new TypeError("invalid stage input fingerprint");
    }
    if (!command || typeof command !== "object" || Array.isArray(command) || containsAbsolutePath(command)) {
      throw new TypeError("stage command must use logical paths");
    }
    if (!inputs || typeof inputs !== "object" || Array.isArray(inputs) || containsAbsolutePath(inputs)) {
      throw new TypeError("stage inputs must use logical resources");
    }

    const existing = await numberedAttemptEntries(runId, stage);
    const attempt = (existing.at(-1)?.attempt || 0) + 1;
    const root = attemptRoot(runsRoot, runId, stage, attempt);
    await fs.mkdir(path.dirname(root), { recursive: true });
    await fs.mkdir(root, { recursive: false });
    await fs.mkdir(path.join(root, "logs"));
    await fs.mkdir(path.join(root, "outputs"));
    const checkpoint = {
      schemaVersion: STAGE_SCHEMA,
      runId,
      stage,
      attempt,
      status: "queued",
      inputFingerprint,
      inputs: cloneJson(inputs),
      outputs: {},
      command: cloneJson(command),
      createdAt: now(),
      startedAt: null,
      finishedAt: null,
      error: null,
    };
    const persisted = await persist(runId, stage, attempt, checkpoint, "stage_attempt_queued");
    return { ...persisted, root };
  }

  async function loadMutable(runId, stage, attempt) {
    assertIdentifier(runId, "run id");
    assertStage(stage);
    const checkpointPath = path.join(attemptRoot(runsRoot, runId, stage, attempt), "checkpoint.json");
    const { checkpoint } = await readCheckpoint(checkpointPath);
    if (!validCheckpoint(checkpoint, runId, stage, attempt)) {
      throw new Error("invalid stage checkpoint");
    }
    if (TERMINAL_STATUSES.has(checkpoint.status)) {
      throw new Error(`terminal attempt cannot be modified: ${stage}/${attemptName(attempt)}`);
    }
    return checkpoint;
  }

  async function transition(runId, stage, attempt, status) {
    if (status !== "running") throw new TypeError(`unsupported stage transition: ${status}`);
    const checkpoint = await loadMutable(runId, stage, attempt);
    if (checkpoint.status !== "queued") {
      throw new Error(`invalid stage transition: ${checkpoint.status} -> ${status}`);
    }
    checkpoint.status = "running";
    checkpoint.startedAt = now();
    return persist(runId, stage, attempt, checkpoint, "stage_attempt_started");
  }

  async function complete(runId, stage, attempt, outputs) {
    const checkpoint = await loadMutable(runId, stage, attempt);
    if (checkpoint.status !== "running") {
      throw new Error(`invalid stage transition: ${checkpoint.status} -> complete`);
    }
    if (!outputs || typeof outputs !== "object" || Array.isArray(outputs) || !Object.keys(outputs).length) {
      throw new TypeError("complete stage requires declared outputs");
    }
    const root = attemptRoot(runsRoot, runId, stage, attempt);
    const hashedOutputs = {};
    for (const [logicalName, relativePath] of Object.entries(outputs).sort(([left], [right]) => left.localeCompare(right))) {
      if (!logicalName || logicalName.includes("/") || logicalName.includes("\\")) {
        throw new TypeError(`invalid stage output name: ${logicalName}`);
      }
      const safePath = safeOutputPath(relativePath);
      const bytes = await fs.readFile(path.join(root, safePath));
      hashedOutputs[logicalName] = { path: safePath, sha256: sha256(bytes) };
    }
    checkpoint.status = "complete";
    checkpoint.outputs = hashedOutputs;
    checkpoint.finishedAt = now();
    checkpoint.error = null;
    return persist(runId, stage, attempt, checkpoint, "stage_attempt_complete");
  }

  async function fail(runId, stage, attempt, error, { status = "failed" } = {}) {
    if (!FAILURE_STATUSES.has(status)) throw new TypeError(`unsupported terminal status: ${status}`);
    const checkpoint = await loadMutable(runId, stage, attempt);
    checkpoint.status = status;
    checkpoint.finishedAt = now();
    checkpoint.error = {
      code: String(error?.code || status),
      message: String(error?.message || "Stage attempt failed"),
    };
    const eventName = status === "cancelled"
      ? "stage_attempt_cancelled"
      : status === "blocked" ? "stage_attempt_blocked" : "stage_attempt_failed";
    return persist(runId, stage, attempt, checkpoint, eventName);
  }

  async function listAttempts(runId, stage) {
    const entries = await numberedAttemptEntries(runId, stage);
    const attempts = [];
    for (const entry of entries) {
      try {
        const { checkpoint } = await readCheckpoint(
          path.join(attemptRoot(runsRoot, runId, stage, entry.attempt), "checkpoint.json"),
        );
        if (!validCheckpoint(checkpoint, runId, stage, entry.attempt)) throw new Error("invalid");
        attempts.push(cloneJson(checkpoint));
      } catch {
        attempts.push({
          schemaVersion: STAGE_SCHEMA,
          runId,
          stage,
          attempt: entry.attempt,
          status: "blocked",
          startedAt: null,
          finishedAt: null,
          error: { code: "checkpoint_invalid", message: "Stage checkpoint is invalid" },
        });
      }
    }
    return attempts;
  }

  async function verifyAttempt(runId, stage, attempt, expectedFingerprint = null) {
    let record;
    try {
      record = await readCheckpoint(
        path.join(attemptRoot(runsRoot, runId, stage, attempt), "checkpoint.json"),
      );
    } catch {
      return { ok: false, code: "checkpoint_invalid", resource: "checkpoint.json" };
    }
    const { bytes, checkpoint } = record;
    if (!validCheckpoint(checkpoint, runId, stage, attempt) || checkpoint.status !== "complete") {
      return { ok: false, code: "checkpoint_invalid", resource: "checkpoint.json" };
    }
    if (expectedFingerprint && checkpoint.inputFingerprint !== expectedFingerprint) {
      return { ok: false, code: "checkpoint_mismatch", resource: "inputFingerprint" };
    }
    for (const [logicalName, output] of Object.entries(checkpoint.outputs)) {
      if (
        !output
        || typeof output !== "object"
        || !SHA256_PATTERN.test(String(output.sha256 || ""))
      ) {
        return { ok: false, code: "checkpoint_invalid", resource: `outputs/${logicalName}` };
      }
      let safePath;
      try {
        safePath = safeOutputPath(output.path);
      } catch {
        return { ok: false, code: "checkpoint_invalid", resource: `outputs/${logicalName}` };
      }
      try {
        const outputBytes = await fs.readFile(path.join(attemptRoot(runsRoot, runId, stage, attempt), safePath));
        if (sha256(outputBytes) !== output.sha256) {
          return { ok: false, code: "output_mismatch", resource: safePath };
        }
      } catch {
        return { ok: false, code: "output_mismatch", resource: safePath };
      }
    }
    return {
      ok: true,
      checkpoint: cloneJson(checkpoint),
      checkpointSha256: sha256(bytes),
    };
  }

  async function verifyContiguous(run, currentInputs) {
    const runId = assertIdentifier(run?.runId, "run id");
    const stages = stageSequence(run?.route?.branch);
    const reused = [];
    for (const stage of stages) {
      const attempts = await listAttempts(runId, stage);
      const invalid = attempts.find((attempt) => attempt.error?.code === "checkpoint_invalid");
      if (invalid) {
        return { ok: false, code: "checkpoint_invalid", stage, resource: "checkpoint.json" };
      }
      const checkpoint = attempts.filter((attempt) => attempt.status === "complete").at(-1);
      if (!checkpoint) return { ok: true, nextStage: stage, reused };
      if (!currentInputs?.[stage]) {
        return { ok: false, code: "checkpoint_invalid", stage, resource: "inputFingerprint" };
      }
      const verification = await verifyAttempt(
        runId,
        stage,
        checkpoint.attempt,
        computeStageFingerprint(currentInputs[stage]),
      );
      if (!verification.ok) {
        return {
          ok: false,
          code: verification.code === "output_mismatch" ? "checkpoint_mismatch" : verification.code,
          stage,
          resource: checkpointResource({ path: verification.resource }),
        };
      }
      reused.push({ stage, attempt: checkpoint.attempt });
    }
    return { ok: true, nextStage: null, reused };
  }

  async function summarize(run) {
    const runId = assertIdentifier(run?.runId, "run id");
    const branch = run?.route?.branch;
    stageSequence(branch);
    const rows = [];
    for (const stage of STAGES) {
      if (branch === "library" && stage === "dock") {
        rows.push({ stage, status: "skipped", attempts: 0 });
        continue;
      }
      const attempts = await listAttempts(runId, stage);
      const latest = attempts.at(-1);
      rows.push({
        stage,
        status: latest?.status || "waiting",
        attempts: attempts.length,
        startedAt: latest?.startedAt || null,
        finishedAt: latest?.finishedAt || null,
        errorCode: latest?.error?.code || null,
      });
    }
    return {
      ok: true,
      schemaVersion: "open-molecule-lab.stage-summary.v0.1",
      runId,
      resumable: ["failed", "cancelled"].includes(run?.status),
      stages: rows,
    };
  }

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
}
