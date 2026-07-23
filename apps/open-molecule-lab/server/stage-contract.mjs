import { createHash } from "node:crypto";
import fs from "node:fs/promises";
import path from "node:path";


export const STAGE_SCHEMA = "open-molecule-lab.stage-attempt.v0.1";
export const CODE_IDENTITY_SCHEMA = "open-molecule-lab.code-identity.v0.1";

function sha256(value) {
  return createHash("sha256").update(value).digest("hex");
}

export function canonicalJson(value) {
  if (value === null || typeof value === "string" || typeof value === "boolean") {
    return JSON.stringify(value);
  }
  if (typeof value === "number") {
    if (!Number.isFinite(value)) throw new TypeError("canonical JSON requires finite numbers");
    return JSON.stringify(value);
  }
  if (Array.isArray(value)) {
    return `[${value.map((item) => canonicalJson(item)).join(",")}]`;
  }
  if (value && typeof value === "object" && Object.getPrototypeOf(value) === Object.prototype) {
    return `{${Object.keys(value)
      .sort()
      .map((key) => `${JSON.stringify(key)}:${canonicalJson(value[key])}`)
      .join(",")}}`;
  }
  throw new TypeError(`unsupported canonical JSON value: ${typeof value}`);
}

export function stageSequence(branch) {
  if (branch === "library") return ["prepare", "score", "report"];
  if (branch === "cascade") return ["prepare", "score", "dock", "report"];
  throw new TypeError(`unsupported route branch: ${branch}`);
}

export function computeStageFingerprint(input) {
  if (!input || typeof input !== "object" || Array.isArray(input)) {
    throw new TypeError("stage fingerprint input must be an object");
  }
  return sha256(Buffer.from(canonicalJson(input), "utf8"));
}

async function collectPythonFiles(root, relative = "") {
  const directory = path.join(root, relative);
  const entries = await fs.readdir(directory, { withFileTypes: true });
  const files = [];
  for (const entry of entries) {
    if (entry.name === "__pycache__" || entry.name.startsWith(".")) continue;
    const child = path.join(relative, entry.name);
    if (entry.isDirectory()) files.push(...await collectPythonFiles(root, child));
    else if (entry.isFile() && entry.name.endsWith(".py")) files.push(child);
  }
  return files;
}

export async function computeCodeIdentity(sourceRoot) {
  const root = path.resolve(sourceRoot);
  const scoringRoot = path.join(root, "scoring");
  const serverRoot = path.join(root, "apps", "open-molecule-lab", "server");
  const relativePaths = (await collectPythonFiles(scoringRoot)).map((relative) =>
    path.posix.join("scoring", relative.split(path.sep).join(path.posix.sep))
  );
  const serverEntries = await fs.readdir(serverRoot, { withFileTypes: true });
  for (const entry of serverEntries) {
    if (entry.isFile() && (entry.name.endsWith(".mjs") || entry.name.endsWith(".py"))) {
      relativePaths.push(path.posix.join("apps/open-molecule-lab/server", entry.name));
    }
  }
  relativePaths.push(
    "requirements.lock.txt",
    "requirements-runtime.txt",
    "apps/open-molecule-lab/package-lock.json",
  );

  const files = [];
  for (const relative of [...new Set(relativePaths)].sort()) {
    const bytes = await fs.readFile(path.join(root, relative));
    files.push({ path: relative, sha256: sha256(bytes) });
  }
  return {
    schemaVersion: CODE_IDENTITY_SCHEMA,
    sha256: sha256(Buffer.from(canonicalJson(files), "utf8")),
    files,
  };
}
