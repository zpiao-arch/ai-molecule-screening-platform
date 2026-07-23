export type CapabilityStatus = "available" | "not_found" | "planned";

export type RouteMode = "auto" | "library" | "cascade";

export type RouteStatus = "ready" | "blocked";

export type StageStatus = "planned" | "skipped" | "blocked";

export type ExecutionStatus = "queued" | "running" | "complete" | "failed" | "blocked" | "cancelled";

export type ResultView = "ranked" | "failed" | "all";

export interface CliHealth {
  available: boolean;
  scoringEntry: string;
  benchmarkEntry: string;
}

export interface LabHealth {
  ok: boolean;
  product: "open-molecule-lab";
  mode: "local_execution";
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
  mode: "plan_only" | "execute";
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
  moleculeSet:
    | { attached: false; expectedCandidateCount: number }
    | { attached: true; id: string; inputSha256: string; nCandidates: number };
  selection: {
    finalSelectionCount: number;
  };
  execution: {
    strictBackends: true;
    worker: "local";
    seed: number;
  };
}

export interface PlanRunSpec extends Omit<RunSpec, "mode" | "moleculeSet"> {
  mode: "plan_only";
  moleculeSet: {
    attached: false;
    expectedCandidateCount: number;
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
  spec: PlanRunSpec;
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

export interface MoleculeSet {
  ok: true;
  schemaVersion: "open-molecule-lab.molecule-set.v0.1";
  moleculeSetId: string;
  name: string;
  license: string;
  inputSha256: string;
  nRows: number;
  columns: string[];
  createdAt: string;
}

export interface PreflightCheck {
  id: string;
  label: string;
  required: boolean;
  status: "passed" | "failed" | "skipped";
  message: string;
  details?: { missing?: string[]; mismatches?: string[] };
}

export interface ExecutionRun {
  ok: true;
  schemaVersion: "open-molecule-lab.run.v0.1";
  runId: string;
  planRunId: string;
  moleculeSetId: string;
  specSha256: string;
  status: ExecutionStatus;
  route: PromptPlan["route"];
  preflight: {
    schemaVersion: "open-molecule-lab.preflight.v0.1";
    ok: boolean;
    routeBranch: "library" | "cascade";
    checkedAt: string;
    checks: PreflightCheck[];
  };
  createdAt: string;
  startedAt: string | null;
  finishedAt: string | null;
  resultSummary?: { nRows: number; nRanked: number; nFailed: number };
  error?: { code: string; message: string; stderrTail?: string; exitCode?: number } | null;
}

export type ResultRow = Record<string, string | number | null> & {
  id: string;
  smiles: string;
  layer1_status?: string;
  layer2_status?: string;
  layer3_status?: string;
  layer4_status?: string;
  final_score?: number | null;
  final_score_dock?: number | null;
  gate_status?: string;
  gate_reason?: string | null;
};

export interface ResultSummary {
  ok: true;
  nRows: number;
  nRanked: number;
  nFailed: number;
  columns: string[];
  view: ResultView;
  offset: number;
  limit: number;
  total: number;
  rankingScoreField: "final_score" | "final_score_dock";
  rows: ResultRow[];
}
