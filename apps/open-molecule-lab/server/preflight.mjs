import { accessSync, constants, existsSync } from "node:fs";
import fs from "node:fs/promises";
import { spawnSync } from "node:child_process";
import path from "node:path";


const requiredModules = ["rdkit", "numpy", "pandas", "pyarrow", "sklearn", "torch", "unimol_tools"];

function check(id, label, required, passed, message, details) {
  return {
    id,
    label,
    required,
    status: passed ? "passed" : required ? "failed" : "skipped",
    message,
    ...(details ? { details } : {}),
  };
}

function parseLastJson(value) {
  try {
    return JSON.parse(String(value || "").trim());
  } catch {
    // Verification helpers may pretty-print JSON; try the last line next.
  }
  const line = String(value || "").trim().split("\n").filter(Boolean).at(-1);
  try {
    return JSON.parse(line || "{}");
  } catch {
    return null;
  }
}

function executable(value) {
  if (!value) return false;
  try {
    accessSync(value, constants.X_OK);
    return true;
  } catch {
    return false;
  }
}

async function cascadeAssets({ sourceRoot, assetRoot, targetId, sminaBin, obabelBin }) {
  if (!assetRoot) {
    return { receptor: false, receptorName: "", smina: false, obabel: false };
  }
  let receptorPath = "";
  try {
    const registry = JSON.parse(
      await fs.readFile(path.join(sourceRoot, "scoring", "receptor_registry.json"), "utf8"),
    );
    const entry = registry.entries?.[targetId];
    if (entry?.pdbqt) {
      receptorPath = path.isAbsolute(entry.pdbqt)
        ? entry.pdbqt
        : path.resolve(assetRoot, "scoring", entry.pdbqt);
    }
  } catch {
    receptorPath = "";
  }
  return {
    receptor: Boolean(receptorPath && existsSync(receptorPath)),
    receptorName: receptorPath ? path.basename(receptorPath) : "",
    smina: executable(sminaBin),
    obabel: executable(obabelBin),
  };
}

export async function runPreflight({
  sourceRoot,
  pythonPath,
  assetRoot,
  routeBranch,
  targetId,
  expectedCandidateCount,
  actualCandidateCount,
  sminaBin,
  obabelBin,
}) {
  const checks = [];
  const probeCode = [
    "import importlib.util,json,sys",
    `mods=${JSON.stringify(requiredModules)}`,
    "print(json.dumps({'version':list(sys.version_info[:3]),'modules':{m:bool(importlib.util.find_spec(m)) for m in mods}}))",
  ].join(";");
  const probe = spawnSync(pythonPath, ["-c", probeCode], {
    cwd: sourceRoot,
    encoding: "utf8",
    timeout: 15_000,
  });
  const python = probe.status === 0 ? parseLastJson(probe.stdout) : null;
  const versionOk = Boolean(python && python.version?.[0] === 3 && python.version?.[1] === 11);
  checks.push(check("python-version", "Python 3.11 runtime", true, versionOk, versionOk ? `Python ${python.version.join(".")}` : "Python 3.11 is required"));
  const missingModules = python
    ? requiredModules.filter((moduleName) => !python.modules?.[moduleName])
    : [...requiredModules];
  checks.push(check(
    "python-modules",
    "Four-level Python modules",
    true,
    missingModules.length === 0,
    missingModules.length === 0 ? "all required modules are importable" : `missing modules: ${missingModules.join(", ")}`,
    missingModules.length ? { missing: missingModules } : undefined,
  ));
  checks.push(check(
    "candidate-count",
    "MoleculeSet candidate count",
    true,
    expectedCandidateCount === actualCandidateCount,
    expectedCandidateCount === actualCandidateCount
      ? `${actualCandidateCount} sealed candidates match the plan`
      : `plan expects ${expectedCandidateCount}, molecule set contains ${actualCandidateCount}`,
  ));

  const manifestPath = assetRoot ? path.join(assetRoot, "ASSET_MANIFEST.json") : "";
  checks.push(check(
    "asset-root",
    "External asset root",
    true,
    Boolean(assetRoot && existsSync(assetRoot)),
    assetRoot && existsSync(assetRoot) ? "external asset root is configured" : "OPEN_MOLECULE_ASSET_ROOT is not configured",
  ));
  let assetVerification = null;
  if (assetRoot && existsSync(manifestPath)) {
    const verify = spawnSync(
      pythonPath,
      [
        path.join(sourceRoot, "scripts", "verify_assets.py"),
        "--asset-root",
        assetRoot,
        "--manifest",
        manifestPath,
      ],
      { cwd: sourceRoot, encoding: "utf8", timeout: 120_000 },
    );
    assetVerification = parseLastJson(verify.stdout);
  }
  checks.push(check(
    "asset-manifest",
    "External asset manifest",
    true,
    Boolean(assetVerification?.ok),
    assetVerification?.ok
      ? `${assetVerification.checked} manifested assets verified`
      : "asset manifest is missing, incomplete, or has SHA-256 mismatches",
    assetVerification && !assetVerification.ok
      ? { missing: assetVerification.missing || [], mismatches: assetVerification.mismatches || [] }
      : undefined,
  ));

  const cascade = await cascadeAssets({ sourceRoot, assetRoot, targetId, sminaBin, obabelBin });
  const needsDocking = routeBranch === "cascade";
  checks.push(check(
    "registered-receptor",
    "Registered receptor",
    needsDocking,
    cascade.receptor,
    needsDocking
      ? cascade.receptor
        ? `registered receptor ${cascade.receptorName} is available`
        : "cascade route requires a registered receptor under the asset root"
      : "library route does not require a receptor",
  ));
  checks.push(check(
    "docking-binaries",
    "smina and Open Babel",
    needsDocking,
    cascade.smina && cascade.obabel,
    needsDocking
      ? cascade.smina && cascade.obabel
        ? "smina and obabel are executable"
        : "cascade route requires executable SMINA_BIN and OBABEL_BIN"
      : "library route does not execute docking",
  ));

  return {
    schemaVersion: "open-molecule-lab.preflight.v0.1",
    ok: checks.every((item) => !item.required || item.status === "passed"),
    routeBranch,
    checkedAt: new Date().toISOString(),
    checks,
  };
}
