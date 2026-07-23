import { spawn, spawnSync } from "node:child_process";
import fs from "node:fs/promises";
import path from "node:path";

import {
  canonicalJson,
  computeCodeIdentity,
  computeStageFingerprint,
  stageSequence,
} from "./stage-contract.mjs";
import { atomicWriteJson, sha256 } from "./store.mjs";


function tail(text, max = 4000) {
  return String(text || "").slice(-max);
}

function readBridgeJson(result) {
  const output = String(result.stdout || "").trim();
  try {
    return JSON.parse(output.split("\n").filter(Boolean).at(-1) || "{}");
  } catch {
    return { ok: false, error: `invalid bridge JSON: ${tail(output)}` };
  }
}

function redacted(value, replacements) {
  let text = String(value || "");
  for (const [absolutePath, token] of replacements) {
    if (absolutePath) text = text.split(absolutePath).join(token);
  }
  return text;
}

async function sanitizeLogs(logRoot, replacements) {
  for (const name of ["stdout.log", "stderr.log"]) {
    const logPath = path.join(logRoot, name);
    const text = await fs.readFile(logPath, "utf8").catch(() => "");
    await fs.writeFile(logPath, redacted(text, replacements), "utf8");
  }
}

async function atomicCopy(source, destination) {
  const temporary = `${destination}.tmp-${process.pid}-${Date.now()}`;
  try {
    await fs.copyFile(source, temporary);
    await fs.rename(temporary, destination);
  } finally {
    await fs.rm(temporary, { force: true });
  }
}

export function terminatePid(pid, signal = "SIGTERM") {
  const numericPid = Number(pid);
  if (!Number.isInteger(numericPid) || numericPid <= 0) return false;
  if (process.platform !== "win32") {
    try {
      process.kill(-numericPid, signal);
      return true;
    } catch (error) {
      if (error?.code !== "ESRCH") return false;
    }
  }
  try {
    process.kill(numericPid, signal);
    return true;
  } catch {
    return false;
  }
}

export function terminateProcessGroup(child, signal = "SIGTERM") {
  if (!child?.pid) return false;
  if (terminatePid(child.pid, signal)) return true;
  try {
    return child.kill(signal);
  } catch {
    return false;
  }
}

export function terminateOwnedPid(pid, { sourceRoot, pythonPath } = {}, signal = "SIGTERM") {
  const numericPid = Number(pid);
  if (!Number.isInteger(numericPid) || numericPid <= 0 || !sourceRoot) return false;
  const probe = spawnSync("ps", ["-p", String(numericPid), "-o", "command="], {
    encoding: "utf8",
    timeout: 2000,
  });
  if (probe.status !== 0) return false;
  const command = String(probe.stdout || "").trim();
  const expectedScript = path.join(path.resolve(sourceRoot), "scoring", "scoring.py");
  if (!command.includes(expectedScript)) return false;
  if (pythonPath && !command.includes(String(pythonPath))) return false;
  return terminatePid(numericPid, signal);
}

function wait(milliseconds) {
  return new Promise((resolve) => setTimeout(resolve, milliseconds));
}

function processGroupExists(pid) {
  const numericPid = Number(pid);
  if (!Number.isInteger(numericPid) || numericPid <= 0) return false;
  try {
    process.kill(process.platform === "win32" ? numericPid : -numericPid, 0);
    return true;
  } catch (error) {
    return error?.code === "EPERM";
  }
}

async function waitForSettled(entries, timeoutMs) {
  if (!entries.length) return true;
  return Promise.race([
    Promise.all(entries.map(([, entry]) => entry.donePromise)).then(() => true),
    wait(timeoutMs).then(() => false),
  ]);
}

function stagePolicy(stage, spec) {
  if (stage === "prepare") return { role: "validated_inputs" };
  if (stage === "score") return { role: "base_four_level_score", strictBackends: true };
  if (stage === "dock") {
    return {
      role: "docking_cascade",
      strictBackends: true,
      dockTopN: Number(spec.policy?.dockTopN || 300),
    };
  }
  return { role: "evidence_report" };
}

export function createStageInputBuilder({ store, stageStore, sourceRoot, assetRoot }) {
  return async function buildStageInputs(run) {
    const runRoot = path.join(store.runsRoot, run.runId);
    const spec = JSON.parse(await fs.readFile(path.join(runRoot, "run-spec.json"), "utf8"));
    const codeIdentity = await computeCodeIdentity(sourceRoot);
    const assetManifestBytes = await fs.readFile(path.join(assetRoot, "ASSET_MANIFEST.json"));
    const common = {
      schemaVersion: "open-molecule-lab.stage-input.v0.1",
      runSpecSha256: run.specSha256,
      moleculeSetSha256: spec.moleculeSet?.inputSha256,
      assetManifestSha256: sha256(assetManifestBytes),
      codeIdentity: codeIdentity.sha256,
      routeBranch: run.route.branch,
    };
    let upstreamCheckpointSha256 = null;
    let upstreamOutputSha256 = null;
    const inputs = {};
    for (const stage of stageSequence(run.route.branch)) {
      inputs[stage] = {
        ...common,
        stage,
        stagePolicy: stagePolicy(stage, spec),
        upstreamCheckpointSha256,
        upstreamOutputSha256,
      };
      const attempts = await stageStore.listAttempts(run.runId, stage);
      const complete = attempts.filter((attempt) => attempt.status === "complete").at(-1);
      if (!complete) break;
      const verification = await stageStore.verifyAttempt(
        run.runId,
        stage,
        complete.attempt,
        computeStageFingerprint(inputs[stage]),
      );
      if (!verification.ok) break;
      upstreamCheckpointSha256 = verification.checkpointSha256;
      upstreamOutputSha256 = sha256(Buffer.from(canonicalJson(verification.checkpoint.outputs), "utf8"));
    }
    return inputs;
  };
}

export function createWorkerManager({
  store,
  stageStore,
  buildStageInputs,
  sourceRoot,
  pythonPath,
  bridgePath,
  assetRoot,
  sminaBin,
  obabelBin,
}) {
  if (!stageStore || typeof buildStageInputs !== "function") {
    throw new TypeError("stageStore and buildStageInputs are required");
  }
  const processes = new Map();

  function replacementsFor(runRoot) {
    return [
      [assetRoot, "<ASSET_ROOT>"],
      [sourceRoot, "<SOURCE_ROOT>"],
      [runRoot, "<RUN_ROOT>"],
    ]
      .filter(([absolutePath]) => absolutePath)
      .sort(([left], [right]) => right.length - left.length);
  }

  function workerEnvironment() {
    return {
      ...process.env,
      FOUR_LEVEL_ASSET_ROOT: assetRoot,
      FOUR_LEVEL_ASSET_MANIFEST: path.join(assetRoot, "ASSET_MANIFEST.json"),
      FOUR_LEVEL_RECEPTOR_REGISTRY: path.join(sourceRoot, "scoring", "receptor_registry.json"),
      ...(sminaBin ? { SMINA_BIN: sminaBin } : {}),
      ...(obabelBin ? { OBABEL_BIN: obabelBin } : {}),
    };
  }

  function logicalCommand(stage, spec) {
    if (stage === "score") {
      return {
        executable: "OPEN_MOLECULE_PYTHON",
        args: [
          "scoring/scoring.py",
          "--input", "inputs/molecules.csv",
          "--output", "outputs/scores.csv",
          "--target", spec.target.id,
          "--mode", "library",
          "--strict-backends",
          "--asset-root", "<OPEN_MOLECULE_ASSET_ROOT>",
        ],
      };
    }
    return {
      executable: "OPEN_MOLECULE_PYTHON",
      args: [
        "scoring/scoring.py",
        "--input", "inputs/molecules.csv",
        "--base-scores", "<VERIFIED_SCORE_OUTPUT>",
        "--output", "outputs/scores.csv",
        "--target", spec.target.id,
        "--mode", "cascade",
        "--strict-backends",
        "--asset-root", "<OPEN_MOLECULE_ASSET_ROOT>",
        "--cascade-top-n", String(spec.policy?.dockTopN || 300),
      ],
    };
  }

  async function latestCompleteOutput(runId, stage, logicalName) {
    const attempts = await stageStore.listAttempts(runId, stage);
    const checkpoint = attempts.filter((attempt) => attempt.status === "complete").at(-1);
    if (!checkpoint) throw new Error(`missing complete ${stage} stage`);
    const verified = await stageStore.verifyAttempt(runId, stage, checkpoint.attempt);
    if (!verified.ok) throw new Error(`${stage} checkpoint verification failed: ${verified.code}`);
    const output = verified.checkpoint.outputs?.[logicalName];
    if (!output?.path) throw new Error(`${stage} checkpoint missing ${logicalName}`);
    return {
      checkpoint: verified.checkpoint,
      path: path.join(
        store.runsRoot,
        runId,
        "stages",
        stage,
        `attempt-${String(checkpoint.attempt).padStart(4, "0")}`,
        output.path,
      ),
    };
  }

  function summarizeResults(outputPath, inputPath) {
    const result = spawnSync(
      pythonPath,
      [bridgePath, "summarize-results", "--input", outputPath, "--expected-input", inputPath],
      { cwd: sourceRoot, encoding: "utf8", timeout: 30_000 },
    );
    return { process: result, summary: readBridgeJson(result) };
  }

  function spawnCommand(args, stdoutHandle, stderrHandle) {
    let child;
    try {
      child = spawn(pythonPath, args, {
        cwd: sourceRoot,
        env: workerEnvironment(),
        detached: process.platform !== "win32",
        stdio: ["ignore", stdoutHandle.fd, stderrHandle.fd],
      });
    } catch (error) {
      return {
        child: null,
        done: Promise.resolve({ child: null, code: null, signal: null, error }),
      };
    }
    const done = new Promise((resolve) => {
      let settled = false;
      const finish = (payload) => {
        if (settled) return;
        settled = true;
        resolve(payload);
      };
      child.once("error", (error) => finish({ child, code: null, signal: null, error }));
      child.once("exit", (code, signal) => finish({ child, code, signal, error: null }));
    });
    return { child, done };
  }

  async function failAttemptAndRun(run, entry, error, { status = "failed", signal = null } = {}) {
    const runRoot = path.join(store.runsRoot, run.runId);
    const replacements = replacementsFor(runRoot);
    const safeError = {
      code: String(error?.code || status),
      message: redacted(error?.message || "Stage attempt failed", replacements),
      ...(error?.exitCode !== undefined ? { exitCode: error.exitCode } : {}),
      ...(signal ? { signal } : {}),
      ...(error?.stderrTail ? { stderrTail: redacted(error.stderrTail, replacements) } : {}),
    };
    if (entry.stage && entry.attempt) {
      try {
        await stageStore.fail(run.runId, entry.stage, entry.attempt, safeError, { status });
      } catch (stageError) {
        if (!/terminal attempt/i.test(String(stageError?.message || stageError))) throw stageError;
      }
    }
    await store.updateRun(
      run.runId,
      {
        status,
        finishedAt: new Date().toISOString(),
        pid: null,
        currentStage: entry.stage,
        currentAttempt: entry.attempt,
        error: safeError,
      },
      { event: status === "cancelled" ? "worker_cancelled" : "worker_failed", status, stage: entry.stage },
    );
    entry.child = null;
    return false;
  }

  async function executeCommandStage(run, spec, stage, entry) {
    const runRoot = path.join(store.runsRoot, run.runId);
    const inputPath = path.join(runRoot, "inputs", "molecules.csv");
    const currentInputs = await buildStageInputs(await store.getRun(run.runId));
    const stageInput = currentInputs[stage];
    if (!stageInput) throw new Error(`stage input is unavailable for ${stage}`);
    const command = logicalCommand(stage, spec);
    const attempt = await stageStore.createAttempt(
      run.runId,
      stage,
      computeStageFingerprint(stageInput),
      command,
      stageInput,
    );
    entry.stage = stage;
    entry.attempt = attempt.attempt;
    await atomicWriteJson(path.join(attempt.root, "command.json"), command);
    await stageStore.transition(run.runId, stage, attempt.attempt, "running");
    if (entry.cancelRequested) {
      return failAttemptAndRun(run, entry, { code: "cancelled", message: "worker cancellation requested" }, {
        status: "cancelled",
      });
    }

    const outputPath = path.join(attempt.root, "outputs", "scores.csv");
    let baseScores = null;
    if (stage === "dock") baseScores = await latestCompleteOutput(run.runId, "score", "scores.csv");
    const actualArgs = [
      path.join(sourceRoot, "scoring", "scoring.py"),
      "--input", inputPath,
      ...(stage === "dock" ? ["--base-scores", baseScores.path] : []),
      "--output", outputPath,
      "--target", spec.target.id,
      "--mode", stage === "dock" ? "cascade" : "library",
      "--strict-backends",
      "--asset-root", assetRoot,
      ...(stage === "dock" ? ["--cascade-top-n", String(spec.policy?.dockTopN || 300)] : []),
    ];
    const stdoutHandle = await fs.open(path.join(attempt.root, "logs", "stdout.log"), "a");
    const stderrHandle = await fs.open(path.join(attempt.root, "logs", "stderr.log"), "a");
    let processResult;
    try {
      await store.updateRun(
        run.runId,
        {
          status: "running",
          currentStage: stage,
          currentAttempt: attempt.attempt,
          pid: null,
        },
        { event: "worker_stage_started", status: "running", stage, attempt: attempt.attempt },
      );
      const spawned = spawnCommand(actualArgs, stdoutHandle, stderrHandle);
      entry.child = spawned.child;
      if (spawned.child?.pid) {
        await store.updateRun(run.runId, { pid: spawned.child.pid });
      }
      processResult = await spawned.done;
    } finally {
      await stdoutHandle.close();
      await stderrHandle.close();
      await sanitizeLogs(path.join(attempt.root, "logs"), replacementsFor(runRoot));
    }
    entry.child = processResult?.child || null;

    if (entry.cancelRequested) {
      return failAttemptAndRun(run, entry, { code: "cancelled", message: "worker cancellation requested" }, {
        status: "cancelled",
        signal: processResult?.signal,
      });
    }
    if (processResult?.error) {
      return failAttemptAndRun(run, entry, {
        code: "worker_spawn_failed",
        message: String(processResult.error.message || processResult.error),
      });
    }
    if (processResult?.code !== 0) {
      const stderr = await fs.readFile(path.join(attempt.root, "logs", "stderr.log"), "utf8").catch(() => "");
      return failAttemptAndRun(run, entry, {
        code: "cli_failed",
        message: `${stage} CLI exited unsuccessfully`,
        exitCode: processResult?.code,
        stderrTail: tail(stderr),
      }, { signal: processResult?.signal });
    }
    const stat = await fs.stat(outputPath).catch(() => null);
    if (!stat?.isFile() || stat.size === 0) {
      return failAttemptAndRun(run, entry, {
        code: "missing_results",
        message: `${stage} CLI exited without a non-empty scores.csv`,
      });
    }
    const { process: bridgeProcess, summary } = summarizeResults(outputPath, inputPath);
    if (bridgeProcess.status !== 0 || !summary.ok) {
      return failAttemptAndRun(run, entry, {
        code: "result_identity_mismatch",
        message: summary.error || `${stage} result validation failed`,
      });
    }
    if (
      stage === "score"
      && (summary.rankingScoreField !== "final_score" || summary.columns?.includes("final_score_dock"))
    ) {
      return failAttemptAndRun(run, entry, {
        code: "base_score_mutated",
        message: "score stage emitted docking columns",
      });
    }
    if (
      stage === "dock"
      && (summary.rankingScoreField !== "final_score_dock" || Number(summary.structureDockingOk || 0) < 1)
    ) {
      return failAttemptAndRun(run, entry, {
        code: "docking_zero_success",
        message: "dock stage produced no successful structure docking rows",
      });
    }
    await stageStore.complete(run.runId, stage, attempt.attempt, { "scores.csv": "outputs/scores.csv" });
    await store.updateRun(
      run.runId,
      { pid: null, currentStage: null, currentAttempt: null },
      { event: "worker_stage_complete", status: "running", stage, attempt: attempt.attempt },
    );
    entry.child = null;
    entry.stage = null;
    entry.attempt = null;
    return true;
  }

  async function executeReportStage(run, spec, entry) {
    const runRoot = path.join(store.runsRoot, run.runId);
    const inputPath = path.join(runRoot, "inputs", "molecules.csv");
    const currentInputs = await buildStageInputs(await store.getRun(run.runId));
    const stageInput = currentInputs.report;
    if (!stageInput) throw new Error("stage input is unavailable for report");
    const command = { executable: "OPEN_MOLECULE_SERVER", args: ["report"] };
    const attempt = await stageStore.createAttempt(
      run.runId,
      "report",
      computeStageFingerprint(stageInput),
      command,
      stageInput,
    );
    entry.stage = "report";
    entry.attempt = attempt.attempt;
    await atomicWriteJson(path.join(attempt.root, "command.json"), command);
    await stageStore.transition(run.runId, "report", attempt.attempt, "running");
    await store.updateRun(
      run.runId,
      { currentStage: "report", currentAttempt: attempt.attempt, pid: null },
      { event: "worker_stage_started", status: "running", stage: "report", attempt: attempt.attempt },
    );
    try {
      if (entry.cancelRequested) throw Object.assign(new Error("worker cancellation requested"), { code: "cancelled" });
      const sourceStage = run.route.branch === "cascade" ? "dock" : "score";
      const verifiedResult = await latestCompleteOutput(run.runId, sourceStage, "scores.csv");
      const attemptScores = path.join(attempt.root, "outputs", "scores.csv");
      const attemptSummary = path.join(attempt.root, "outputs", "summary.json");
      await atomicCopy(verifiedResult.path, attemptScores);
      const { process: bridgeProcess, summary } = summarizeResults(attemptScores, inputPath);
      if (bridgeProcess.status !== 0 || !summary.ok) {
        throw Object.assign(new Error(summary.error || "report result validation failed"), {
          code: "result_identity_mismatch",
        });
      }
      if (
        run.route.branch === "cascade"
        && (summary.rankingScoreField !== "final_score_dock" || Number(summary.structureDockingOk || 0) < 1)
      ) {
        throw Object.assign(new Error("report has no successful docking rows"), { code: "docking_zero_success" });
      }
      const summaryText = `${JSON.stringify(summary, null, 2)}\n`;
      if (replacementsFor(runRoot).some(([absolutePath]) => summaryText.includes(absolutePath))) {
        throw Object.assign(new Error("report summary contains an absolute path"), { code: "evidence_incomplete" });
      }
      await atomicWriteJson(attemptSummary, summary);
      await atomicCopy(attemptScores, path.join(runRoot, "results", "scores.csv"));
      await atomicCopy(attemptSummary, path.join(runRoot, "results", "summary.json"));
      await stageStore.complete(run.runId, "report", attempt.attempt, {
        "scores.csv": "outputs/scores.csv",
        "summary.json": "outputs/summary.json",
      });
      await store.updateRun(
        run.runId,
        {
          status: "complete",
          finishedAt: new Date().toISOString(),
          pid: null,
          currentStage: null,
          currentAttempt: null,
          error: null,
          resultSummary: {
            nRows: summary.nRows,
            nRanked: summary.nRanked,
            nFailed: summary.nFailed,
            rankingScoreField: summary.rankingScoreField,
            structureDockingOk: summary.structureDockingOk,
          },
        },
        {
          event: "worker_complete",
          status: "complete",
          nRows: summary.nRows,
          nRanked: summary.nRanked,
          nFailed: summary.nFailed,
        },
      );
      const integrity = typeof store.verifyManifest === "function"
        ? await store.verifyManifest(runRoot)
        : { ok: true };
      if (!integrity.ok) {
        await store.updateRun(
          run.runId,
          {
            status: "failed",
            finishedAt: new Date().toISOString(),
            error: { code: "evidence_incomplete", message: "Run evidence manifest verification failed" },
          },
          { event: "worker_failed", status: "failed", stage: "report" },
        );
        return false;
      }
      entry.stage = null;
      entry.attempt = null;
      return true;
    } catch (error) {
      return failAttemptAndRun(
        run,
        entry,
        { code: error?.code || "evidence_incomplete", message: error?.message || String(error) },
        { status: error?.code === "cancelled" ? "cancelled" : "failed" },
      );
    }
  }

  async function executeSequence(run, startStage, entry) {
    const runRoot = path.join(store.runsRoot, run.runId);
    try {
      const spec = JSON.parse(await fs.readFile(path.join(runRoot, "run-spec.json"), "utf8"));
      const sequence = stageSequence(run.route.branch);
      const startIndex = sequence.indexOf(startStage);
      if (startIndex < 0 || startStage === "prepare") {
        throw Object.assign(new Error(`invalid worker start stage: ${startStage}`), { code: "checkpoint_invalid" });
      }
      for (const stage of sequence.slice(startIndex)) {
        const ok = stage === "report"
          ? await executeReportStage(run, spec, entry)
          : await executeCommandStage(run, spec, stage, entry);
        if (!ok) return;
      }
    } catch (error) {
      await failAttemptAndRun(run, entry, {
        code: error?.code || "worker_error",
        message: error?.message || String(error),
      });
    }
  }

  async function start(run, { startStage = null } = {}) {
    if (run.status !== "queued") return run;
    if (processes.has(run.runId)) return store.getRun(run.runId);
    const currentInputs = await buildStageInputs(run);
    const verification = await stageStore.verifyContiguous(run, currentInputs);
    if (!verification.ok) {
      return store.updateRun(
        run.runId,
        {
          status: "blocked",
          finishedAt: new Date().toISOString(),
          error: {
            code: verification.code,
            message: "Verified stage checkpoint mismatch",
            stage: verification.stage,
            resource: verification.resource,
          },
        },
        { event: "resume_blocked", status: "blocked", stage: verification.stage },
      );
    }
    const selectedStage = startStage || verification.nextStage;
    if (!selectedStage || (startStage && startStage !== verification.nextStage)) {
      return store.updateRun(
        run.runId,
        {
          status: "blocked",
          finishedAt: new Date().toISOString(),
          error: { code: "checkpoint_invalid", message: "Requested start stage is not the next verified stage" },
        },
        { event: "resume_blocked", status: "blocked", stage: selectedStage },
      );
    }

    let resolveDone;
    const donePromise = new Promise((resolve) => { resolveDone = resolve; });
    const entry = {
      child: null,
      stage: null,
      attempt: null,
      cancelRequested: false,
      donePromise,
    };
    processes.set(run.runId, entry);
    const running = await store.updateRun(
      run.runId,
      {
        status: "running",
        startedAt: run.startedAt || new Date().toISOString(),
        finishedAt: null,
        currentStage: selectedStage,
        currentAttempt: null,
        pid: null,
        error: null,
      },
      { event: "worker_started", status: "running", stage: selectedStage },
    );
    void executeSequence(running, selectedStage, entry).finally(() => {
      processes.delete(run.runId);
      resolveDone();
    });
    return running;
  }

  async function cancel(runId) {
    const entry = processes.get(runId);
    if (!entry) return store.getRun(runId);
    entry.cancelRequested = true;
    await store.appendEvent(runId, {
      event: "cancel_requested",
      status: "running",
      stage: entry.stage,
      attempt: entry.attempt,
    });
    if (entry.child) terminateProcessGroup(entry.child, "SIGTERM");
    return store.getRun(runId);
  }

  async function shutdown() {
    const entries = [...processes.entries()];
    for (const [, entry] of entries) {
      entry.cancelRequested = true;
      if (entry.child) terminateProcessGroup(entry.child, "SIGTERM");
    }
    if (await waitForSettled(entries, 3000)) return;
    for (const [, entry] of entries) {
      if (entry.child && entry.child.exitCode === null && entry.child.signalCode === null) {
        terminateProcessGroup(entry.child, "SIGKILL");
      }
    }
    await waitForSettled(entries, 3000);
  }

  async function terminateOrphan(pid, signal = "SIGTERM") {
    if (!terminateOwnedPid(pid, { sourceRoot, pythonPath }, signal)) return false;
    for (let attempt = 0; attempt < 20; attempt += 1) {
      if (!processGroupExists(pid)) return true;
      await wait(50);
    }
    terminatePid(pid, "SIGKILL");
    for (let attempt = 0; attempt < 20; attempt += 1) {
      if (!processGroupExists(pid)) return true;
      await wait(50);
    }
    return !processGroupExists(pid);
  }

  return {
    start,
    cancel,
    shutdown,
    terminateOrphan,
    hasProcess: (runId) => processes.has(runId),
  };
}
