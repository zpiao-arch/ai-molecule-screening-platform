import type { AttemptStatus, StageSummary, StageSummaryRow } from "./types";


const attemptStatus: AttemptStatus = "cancelled";
const stageRow: StageSummaryRow = {
  stage: "dock",
  status: attemptStatus,
  attempts: 2,
  startedAt: "2026-07-23T00:00:00.000Z",
  finishedAt: "2026-07-23T00:01:00.000Z",
  errorCode: "worker_interrupted",
};

export const stageSummaryContract: StageSummary = {
  ok: true,
  schemaVersion: "open-molecule-lab.stage-summary.v0.1",
  runId: "run_contract",
  resumable: true,
  stages: [stageRow],
};
