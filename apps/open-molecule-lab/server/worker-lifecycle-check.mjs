import { spawn } from "node:child_process";
import fs from "node:fs/promises";
import os from "node:os";
import path from "node:path";

import { createWorkerManager, terminateOwnedPid, terminateProcessGroup } from "./worker.mjs";


function assert(condition, message) {
  if (!condition) throw new Error(message);
}

function processExists(pid) {
  try {
    process.kill(pid, 0);
    return true;
  } catch (error) {
    return error?.code === "EPERM";
  }
}

async function waitForExit(pid) {
  for (let attempt = 0; attempt < 100; attempt += 1) {
    if (!processExists(pid)) return;
    await new Promise((resolve) => setTimeout(resolve, 25));
  }
  throw new Error(`process ${pid} survived group termination`);
}

const child = spawn(
  process.execPath,
  [
    "-e",
    [
      "const { spawn } = require('node:child_process');",
      "const grandchild = spawn(process.execPath, ['-e', 'setInterval(() => {}, 1000)'], { stdio: 'ignore' });",
      "process.stdout.write(String(grandchild.pid) + '\\n');",
      "setInterval(() => {}, 1000);",
    ].join(""),
  ],
  { detached: true, stdio: ["ignore", "pipe", "inherit"] },
);

const grandchildPid = Number(await new Promise((resolve, reject) => {
  child.stdout.once("data", (chunk) => resolve(String(chunk).trim()));
  child.once("error", reject);
}));

try {
  assert(Number.isInteger(child.pid), "worker PID missing");
  assert(Number.isInteger(grandchildPid), "grandchild PID missing");
  assert(terminateProcessGroup(child, "SIGTERM"), "group termination was not attempted");
  await Promise.all([waitForExit(child.pid), waitForExit(grandchildPid)]);

  const unrelated = spawn(process.execPath, ["-e", "setInterval(() => {}, 1000)"], {
    detached: true,
    stdio: "ignore",
  });
  try {
    assert(
      terminateOwnedPid(unrelated.pid, { sourceRoot: "/tmp/not-this-worker", pythonPath: "/tmp/not-python" }) === false,
      "orphan recovery accepted an unrelated reused PID",
    );
    assert(processExists(unrelated.pid), "orphan recovery terminated an unrelated process");
  } finally {
    terminateProcessGroup(unrelated, "SIGKILL");
  }

  const fixtureRoot = await fs.mkdtemp(path.join(os.tmpdir(), "open-molecule-worker-lifecycle-"));
  const sourceRoot = path.join(fixtureRoot, "source");
  const runId = "run_20260723000000_fixture";
  const runRoot = path.join(fixtureRoot, "runs", runId);
  await fs.mkdir(path.join(sourceRoot, "scoring"), { recursive: true });
  await fs.mkdir(path.join(runRoot, "inputs"), { recursive: true });
  await fs.mkdir(path.join(runRoot, "results"), { recursive: true });
  await fs.mkdir(path.join(runRoot, "logs"), { recursive: true });
  await fs.writeFile(
    path.join(sourceRoot, "scoring", "scoring.py"),
    "import subprocess, sys, time\nsubprocess.Popen([sys.executable, '-c', 'import time; time.sleep(60)'])\ntime.sleep(60)\n",
    "utf8",
  );
  await fs.writeFile(
    path.join(runRoot, "run-spec.json"),
    `${JSON.stringify({ target: { id: "fixture" } })}\n`,
    "utf8",
  );
  const run = { runId, status: "queued", route: { branch: "library" } };
  const store = {
    runsRoot: path.join(fixtureRoot, "runs"),
    run,
    async refreshManifest() {},
    async appendEvent() {},
    async getRun() { return this.run; },
    async updateRun(_runId, patch) {
      if (patch.status === "cancelled") await new Promise((resolve) => setTimeout(resolve, 800));
      this.run = { ...this.run, ...patch };
      return this.run;
    },
  };
  const manager = createWorkerManager({
    store,
    sourceRoot,
    pythonPath: process.env.OPEN_MOLECULE_PYTHON || "python3",
    bridgePath: path.join(sourceRoot, "missing-bridge.py"),
    assetRoot: fixtureRoot,
    sminaBin: "",
    obabelBin: "",
  });
  try {
    await manager.start(run);
    await manager.shutdown();
    assert(store.run.status === "cancelled", "shutdown returned before the cancelled state was persisted");
  } finally {
    await fs.rm(fixtureRoot, { recursive: true, force: true });
  }

  process.stdout.write(`${JSON.stringify({ ok: true, pid: child.pid, grandchildPid }, null, 2)}\n`);
} finally {
  if (processExists(child.pid)) process.kill(child.pid, "SIGKILL");
  if (processExists(grandchildPid)) process.kill(grandchildPid, "SIGKILL");
}
