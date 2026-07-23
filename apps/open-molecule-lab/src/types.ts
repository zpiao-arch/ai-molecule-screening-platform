export type CapabilityStatus = "available" | "not_found" | "planned";

export type RouteMode = "auto" | "library" | "cascade";

export type RouteStatus = "ready" | "blocked";

export type StageStatus = "planned" | "skipped" | "blocked";

export interface CliHealth {
  available: boolean;
  scoringEntry: string;
  benchmarkEntry: string;
}

export interface LabHealth {
  ok: boolean;
  product: "open-molecule-lab";
  mode: "plan_only";
  sourceRoot: string;
  cli: CliHealth;
  assetManifestPresent: boolean;
}

export interface ModelStatusRow {
  id: string;
  label: string;
  layer: string;
  status: CapabilityStatus;
  role: string;
  requirement: string;
}

export interface RunSpec {
  schemaVersion: "open-molecule-lab.run-spec.v0.1";
  mode: "plan_only";
  project: {
    name: string;
    researchPrompt: string;
  };
  target: {
    id: string;
    text: string;
    requestedRoute: RouteMode;
    resolvedBranch: "library" | "cascade";
  };
  moleculeSet: {
    attached: false;
    expectedCandidateCount: number;
  };
  selection: {
    finalSelectionCount: number;
  };
  execution: {
    strictBackends: true;
    worker: "local";
    seed: number;
  };
}

export interface PlannedStage {
  id: string;
  label: string;
  status: StageStatus;
  reason: string;
}

export interface PromptPlan {
  ok: true;
  mode: "plan_only";
  runId: string;
  spec: RunSpec;
  route: {
    branch: "library" | "cascade";
    status: RouteStatus;
    receptorAvailable: boolean;
    rationale: string;
  };
  stages: PlannedStage[];
  assetRequirements: string[];
  executionBoundary: string;
  equivalentCli: string;
  bundle: {
    relativeRoot: string;
    files: string[];
    manifestSha256: string;
  };
}

export interface PromptError {
  ok: false;
  error: {
    code: string;
    message: string;
    field?: string;
  };
}
