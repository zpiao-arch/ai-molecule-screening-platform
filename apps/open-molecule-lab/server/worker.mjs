import { spawn, spawnSync } from "node:child_process";
import fs from "node:fs/promises";
import path from "node:path";


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

async function sanitizeLogs(runRoot, replacements) {
  for (const name of ["stdout.log", "stderr.log"]) {
    const logPath = path.join(runRoot, "logs", name);
    let text = await fs.readFile(logPath, "utf8").catch(() => "");
    for (const [absolutePath, token] of replacements) {
      if (absolutePath) text = text.split(absolutePath).join(token);
    }
    await fs.writeFile(logPath, text, "utf8");
  }
}

export function createWorkerManager({ store, sourceRoot, pythonPath, bridgePath, assetRoot, sminaBin, obabelBin }) {
  const processes = new Map();

  async function start(run) {
    if (run.status !== "queued") return run;
    const runRoot = path.join(store.runsRoot, run.runId);
    const spec = JSON.parse(await fs.readFile(path.join(runRoot, "run-spec.json"), "utf8"));
    const inputPath = path.join(runRoot, "inputs", "molecules.csv");
    const outputPath = path.join(runRoot, "results", "scores.csv");
    const logReplacements = [
      [assetRoot, "<ASSET_ROOT>"],
      [sourceRoot, "<SOURCE_ROOT>"],
      [runRoot, "<RUN_ROOT>"],
    ].sort(([left], [right]) => right.length - left.length);
    const logicalArgs = [
      "scoring/scoring.py",
      "--input",
      "inputs/molecules.csv",
      "--output",
      "results/scores.csv",
      "--target",
      spec.target.id,
      "--mode",
      run.route.branch,
      "--strict-backends",
      "--asset-root",
      "<OPEN_MOLECULE_ASSET_ROOT>",
      ...(run.route.branch === "cascade" ? ["--cascade-top-n", "300"] : []),
    ];
    await fs.writeFile(
      path.join(runRoot, "command.json"),
      `${JSON.stringify({ executable: "OPEN_MOLECULE_PYTHON", args: logicalArgs }, null, 2)}\n`,
      "utf8",
    );
    await store.refreshManifest(runRoot);

    const stdoutHandle = await fs.open(path.join(runRoot, "logs", "stdout.log"), "a");
    const stderrHandle = await fs.open(path.join(runRoot, "logs", "stderr.log"), "a");
    const environment = {
      ...process.env,
      FOUR_LEVEL_ASSET_ROOT: assetRoot,
      FOUR_LEVEL_ASSET_MANIFEST: path.join(assetRoot, "ASSET_MANIFEST.json"),
      FOUR_LEVEL_RECEPTOR_REGISTRY: path.join(sourceRoot, "scoring", "receptor_registry.json"),
      ...(sminaBin ? { SMINA_BIN: sminaBin } : {}),
      ...(obabelBin ? { OBABEL_BIN: obabelBin } : {}),
    };
    let child;
    try {
      child = spawn(
        pythonPath,
        [
          path.join(sourceRoot, "scoring", "scoring.py"),
          "--input",
          inputPath,
          "--output",
          outputPath,
          "--target",
          spec.target.id,
          "--mode",
          run.route.branch,
          "--strict-backends",
          "--asset-root",
          assetRoot,
          ...(run.route.branch === "cascade" ? ["--cascade-top-n", "300"] : []),
        ],
        { cwd: sourceRoot, env: environment, stdio: ["ignore", stdoutHandle.fd, stderrHandle.fd] },
      );
    } catch (error) {
      await stdoutHandle.close();
      await stderrHandle.close();
      return store.updateRun(
        run.runId,
        {
          status: "failed",
          finishedAt: new Date().toISOString(),
          error: { code: "worker_spawn_failed", message: String(error.message || error) },
        },
        { event: "worker_spawn_failed", status: "failed" },
      );
    }

    let settled = false;
    processes.set(run.runId, { child, cancelRequested: false });
    const running = await store.updateRun(
      run.runId,
      { status: "running", startedAt: new Date().toISOString(), pid: child.pid },
      { event: "worker_started", status: "running" },
    );

    child.once("error", async (error) => {
      if (settled) return;
      settled = true;
      await stdoutHandle.close();
      await stderrHandle.close();
      await sanitizeLogs(runRoot, logReplacements);
      processes.delete(run.runId);
      await store.updateRun(
        run.runId,
        {
          status: "failed",
          finishedAt: new Date().toISOString(),
          pid: null,
          error: { code: "worker_error", message: String(error.message || error) },
        },
        { event: "worker_error", status: "failed" },
      );
    });

    child.once("exit", async (code, signal) => {
      if (settled) return;
      settled = true;
      await stdoutHandle.close();
      await stderrHandle.close();
      await sanitizeLogs(runRoot, logReplacements);
      const entry = processes.get(run.runId);
      processes.delete(run.runId);
      const finishedAt = new Date().toISOString();
      if (entry?.cancelRequested) {
        await store.updateRun(
          run.runId,
          {
            status: "cancelled",
            finishedAt,
            pid: null,
            error: { code: "cancelled", message: "worker cancellation requested" },
          },
          { event: "worker_cancelled", status: "cancelled", signal: signal || null },
        );
        return;
      }
      if (code !== 0) {
        const stderr = await fs.readFile(path.join(runRoot, "logs", "stderr.log"), "utf8").catch(() => "");
        await store.updateRun(
          run.runId,
          {
            status: "failed",
            finishedAt,
            pid: null,
            error: { code: "cli_failed", exitCode: code, signal: signal || null, stderrTail: tail(stderr) },
          },
          { event: "worker_failed", status: "failed", exitCode: code, signal: signal || null },
        );
        return;
      }
      const stat = await fs.stat(outputPath).catch(() => null);
      if (!stat || stat.size === 0) {
        await store.updateRun(
          run.runId,
          {
            status: "failed",
            finishedAt,
            pid: null,
            error: { code: "missing_results", message: "CLI exited without a non-empty scores.csv" },
          },
          { event: "worker_failed", status: "failed", reason: "missing_results" },
        );
        return;
      }
      const summaryResult = spawnSync(
        pythonPath,
        [bridgePath, "summarize-results", "--input", outputPath],
        { cwd: sourceRoot, encoding: "utf8", timeout: 30_000 },
      );
      const summary = readBridgeJson(summaryResult);
      if (summaryResult.status !== 0 || !summary.ok) {
        await store.updateRun(
          run.runId,
          {
            status: "failed",
            finishedAt,
            pid: null,
            error: { code: "result_summary_failed", message: summary.error || "result summary failed" },
          },
          { event: "worker_failed", status: "failed", reason: "result_summary_failed" },
        );
        return;
      }
      await fs.writeFile(path.join(runRoot, "results", "summary.json"), `${JSON.stringify(summary, null, 2)}\n`, "utf8");
      await store.updateRun(
        run.runId,
        {
          status: "complete",
          finishedAt,
          pid: null,
          resultSummary: { nRows: summary.nRows, nRanked: summary.nRanked, nFailed: summary.nFailed },
        },
        {
          event: "worker_complete",
          status: "complete",
          nRows: summary.nRows,
          nRanked: summary.nRanked,
          nFailed: summary.nFailed,
        },
      );
    });
    return running;
  }

  async function cancel(runId) {
    const entry = processes.get(runId);
    if (!entry) return store.getRun(runId);
    entry.cancelRequested = true;
    await store.appendEvent(runId, { event: "cancel_requested", status: "running" });
    entry.child.kill("SIGTERM");
    return store.getRun(runId);
  }

  return { start, cancel, hasProcess: (runId) => processes.has(runId) };
}
