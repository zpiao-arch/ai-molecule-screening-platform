import {
  Activity,
  AlertTriangle,
  ArrowRight,
  ChevronDown,
  ChevronLeft,
  ChevronRight,
  ClipboardList,
  FileCheck2,
  FlaskConical,
  Loader2,
  Play,
  RotateCcw,
  ShieldCheck,
  Sparkles,
  Square,
  TerminalSquare,
  Upload,
} from "lucide-react";
import { useEffect, useMemo, useRef, useState } from "react";
import type { ChangeEvent, FormEvent, ReactNode } from "react";

import type {
  LabHealth,
  ExecutionRun,
  MoleculeSet,
  ModelStatusRow,
  PromptError,
  PromptPlan,
  ResultSummary,
  ResultView,
  RouteMode,
  StageSummary,
  StageSummaryRow,
} from "./types";

const scientificBoundary =
  "没有附加分子集并完成严格本机运行前，界面只显示可审计的 RunSpec；任何排序都不构成分子评分、对接证据、药效或安全性结论。";

const defaultPrompt =
  "为流感 A 型神经氨酸酶设计一次四级候选筛选。先检查本地资产，使用 1,000 个候选，最终保留 10 个分子供后续复核。";

const fallbackModels: ModelStatusRow[] = [
  {
    id: "l1-rdkit",
    label: "RDKit molecular quality",
    layer: "L1",
    status: "planned",
    role: "分子质量、描述符和结构合法性",
    requirement: "Python RDKit runtime",
  },
  {
    id: "l2-bindingdb",
    label: "BindingDB target binding",
    layer: "L2",
    status: "planned",
    role: "靶点感知的结合概率粗筛",
    requirement: "BindingDB L2 model asset",
  },
  {
    id: "l3-admet",
    label: "ADMET safety panel",
    layer: "L3",
    status: "planned",
    role: "ADMET 与毒性风险证据",
    requirement: "ADMET model assets",
  },
  {
    id: "l4-unimol",
    label: "UniMol reference similarity",
    layer: "L4",
    status: "planned",
    role: "三维表征与参考药物相似度",
    requirement: "UniMol weights and references",
  },
  {
    id: "dock-smina",
    label: "smina docking cascade",
    layer: "Dock",
    status: "planned",
    role: "有受体时对 L2 头部候选进行精排",
    requirement: "registered receptor, smina and obabel",
  },
];

function formatNumber(value: number) {
  return value.toLocaleString();
}

const stageLabels: Record<StageSummaryRow["stage"], string> = {
  prepare: "准备输入",
  score: "四级评分",
  dock: "Docking cascade",
  report: "证据报告",
};

function formatStageTime(value?: string | null) {
  if (!value) return "";
  const parsed = new Date(value);
  if (Number.isNaN(parsed.getTime())) return "";
  return parsed.toLocaleTimeString([], { hour: "2-digit", minute: "2-digit", second: "2-digit" });
}

function toneForStatus(status?: string) {
  const value = String(status || "").toLowerCase();
  if (value.includes("available") || value.includes("ready") || value.includes("passed") || value.includes("complete")) return "good";
  if (value.includes("blocked") || value.includes("not_found") || value.includes("error") || value.includes("failed")) return "bad";
  if (value.includes("skipped") || value.includes("review") || value.includes("cancelled")) return "warn";
  if (value.includes("running") || value.includes("queued") || value.includes("planned")) return "info";
  return "idle";
}

function App() {
  const [prompt, setPrompt] = useState(defaultPrompt);
  const [target, setTarget] = useState("CHEMBL2051");
  const [candidatePool, setCandidatePool] = useState(1000);
  const [finalSelectionCount, setFinalSelectionCount] = useState(10);
  const [routeMode, setRouteMode] = useState<RouteMode>("auto");
  const [health, setHealth] = useState<LabHealth | null>(null);
  const [models, setModels] = useState<ModelStatusRow[]>(fallbackModels);
  const [plan, setPlan] = useState<PromptPlan | null>(null);
  const [serviceError, setServiceError] = useState("");
  const [submitError, setSubmitError] = useState("");
  const [uploadError, setUploadError] = useState("");
  const [runError, setRunError] = useState("");
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [isUploading, setIsUploading] = useState(false);
  const [moleculeSet, setMoleculeSet] = useState<MoleculeSet | null>(null);
  const [run, setRun] = useState<ExecutionRun | null>(null);
  const [stageSummary, setStageSummary] = useState<StageSummary | null>(null);
  const [isResuming, setIsResuming] = useState(false);
  const [results, setResults] = useState<ResultSummary | null>(null);
  const [resultsLoading, setResultsLoading] = useState(false);
  const [auditOpen, setAuditOpen] = useState(false);
  const resultsRequestRef = useRef(0);
  const stageRequestRef = useRef(0);

  const runnerTone = health?.ok ? "good" : health ? "bad" : "info";
  const runnerLabel = health?.ok ? "local runner online" : health ? "server unavailable" : "connecting";
  const routeTone = toneForStatus(plan?.route.status);
  const runTone = toneForStatus(run?.status);
  const runActive = run?.status === "queued" || run?.status === "running";
  const canResume = Boolean(run && stageSummary?.runId === run.runId && stageSummary.resumable && !runActive);
  const countMatches = Boolean(plan && moleculeSet && plan.spec.moleculeSet.expectedCandidateCount === moleculeSet.nRows);
  const canRun = Boolean(
    health?.ok && plan && moleculeSet && countMatches && plan.route.status !== "blocked" && !runActive,
  );
  const equivalentCli = useMemo(
    () =>
      plan?.equivalentCli ||
      [
        "four-level-molecule \\",
        "  --input <molecule-set.csv> \\",
        `  --target ${target || "<target>"} \\`,
        "  --strict-backends \\",
        "  --output <run-dir>/scores.csv",
      ].join("\n"),
    [plan?.equivalentCli, target],
  );

  useEffect(() => {
    let cancelled = false;

    async function load() {
      try {
        const [healthResponse, modelResponse] = await Promise.all([
          fetch("/api/health", { cache: "no-store" }),
          fetch("/api/model-status", { cache: "no-store" }),
        ]);
        if (!healthResponse.ok) throw new Error(`health HTTP ${healthResponse.status}`);
        if (!modelResponse.ok) throw new Error(`model-status HTTP ${modelResponse.status}`);
        const healthPayload = (await healthResponse.json()) as LabHealth;
        const modelPayload = (await modelResponse.json()) as { models: ModelStatusRow[] };
        if (!cancelled) {
          setHealth(healthPayload);
          setModels(modelPayload.models?.length ? modelPayload.models : fallbackModels);
          setServiceError("");
        }
      } catch (error) {
        if (!cancelled) {
          setHealth(null);
          setModels(fallbackModels);
          setServiceError(error instanceof Error ? error.message : String(error));
        }
      }
    }

    void load();
    return () => {
      cancelled = true;
    };
  }, []);

  async function submitPrompt(event: FormEvent) {
    event.preventDefault();
    setSubmitError("");
    setIsSubmitting(true);
    try {
      const response = await fetch("/api/prompt-plan", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          prompt,
          target,
          candidatePool,
          finalSelectionCount,
          routeMode,
        }),
      });
      const payload = (await response.json()) as PromptPlan | PromptError;
      if (!response.ok || !payload.ok) {
        const message = payload.ok ? `prompt-plan HTTP ${response.status}` : payload.error.message;
        throw new Error(message);
      }
      resultsRequestRef.current += 1;
      stageRequestRef.current += 1;
      setPlan(payload);
      setRun(null);
      setStageSummary(null);
      setResults(null);
      setResultsLoading(false);
      setRunError("");
      setAuditOpen(false);
    } catch (error) {
      setSubmitError(error instanceof Error ? error.message : String(error));
    } finally {
      setIsSubmitting(false);
    }
  }

  async function uploadMoleculeSet(event: ChangeEvent<HTMLInputElement>) {
    const file = event.currentTarget.files?.[0];
    if (!file) return;
    resultsRequestRef.current += 1;
    stageRequestRef.current += 1;
    setResults(null);
    setResultsLoading(false);
    setUploadError("");
    setIsUploading(true);
    try {
      const csvText = await file.text();
      const response = await fetch("/api/molecule-sets", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ name: file.name, csvText, license: "user-supplied" }),
      });
      const payload = (await response.json()) as MoleculeSet | PromptError;
      if (!response.ok || !payload.ok) {
        throw new Error(payload.ok ? `molecule-sets HTTP ${response.status}` : payload.error.message);
      }
      setMoleculeSet(payload);
      setRun(null);
      setStageSummary(null);
    } catch (error) {
      setMoleculeSet(null);
      setUploadError(error instanceof Error ? error.message : String(error));
    } finally {
      setIsUploading(false);
      event.currentTarget.value = "";
    }
  }

  async function launchRun() {
    if (!plan || !moleculeSet || !canRun) return;
    resultsRequestRef.current += 1;
    stageRequestRef.current += 1;
    setRunError("");
    setResults(null);
    setStageSummary(null);
    setResultsLoading(false);
    try {
      const response = await fetch("/api/runs", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ planRunId: plan.runId, moleculeSetId: moleculeSet.moleculeSetId }),
      });
      const payload = (await response.json()) as ExecutionRun | PromptError;
      if (!response.ok || !payload.ok) {
        throw new Error(payload.ok ? `runs HTTP ${response.status}` : payload.error.message);
      }
      setRun(payload);
    } catch (error) {
      setRunError(error instanceof Error ? error.message : String(error));
    }
  }

  async function cancelRun() {
    if (!run || !runActive) return;
    try {
      const response = await fetch(`/api/runs/${run.runId}/cancel`, { method: "POST" });
      const payload = (await response.json()) as ExecutionRun | PromptError;
      if (!response.ok || !payload.ok) throw new Error(payload.ok ? `cancel HTTP ${response.status}` : payload.error.message);
      setRun(payload);
    } catch (error) {
      setRunError(error instanceof Error ? error.message : String(error));
    }
  }

  async function resumeRun() {
    if (!run || !canResume || isResuming) return;
    setRunError("");
    setIsResuming(true);
    resultsRequestRef.current += 1;
    setResults(null);
    try {
      const response = await fetch(`/api/runs/${run.runId}/resume`, { method: "POST" });
      const payload = (await response.json()) as ExecutionRun | PromptError;
      if (!payload.ok) throw new Error(payload.error.message);
      if (response.status !== 202 && response.status !== 409) {
        throw new Error(`resume HTTP ${response.status}`);
      }
      setRun(payload);
      const requestId = ++stageRequestRef.current;
      const stagesResponse = await fetch(`/api/runs/${payload.runId}/stages`, { cache: "no-store" });
      if (stagesResponse.ok) {
        const stagesPayload = (await stagesResponse.json()) as StageSummary;
        if (requestId === stageRequestRef.current && stagesPayload.runId === payload.runId) {
          setStageSummary(stagesPayload);
        }
      }
    } catch (error) {
      setRunError(error instanceof Error ? error.message : String(error));
    } finally {
      setIsResuming(false);
    }
  }

  async function loadResults(runId: string, view: ResultView = "ranked", offset = 0) {
    const requestId = ++resultsRequestRef.current;
    setResultsLoading(true);
    try {
      const params = new URLSearchParams({ view, offset: String(offset), limit: "50" });
      const response = await fetch(`/api/runs/${runId}/results?${params}`, { cache: "no-store" });
      if (!response.ok) throw new Error(`results HTTP ${response.status}`);
      const payload = (await response.json()) as ResultSummary;
      if (requestId === resultsRequestRef.current) setResults(payload);
    } catch (error) {
      if (requestId === resultsRequestRef.current) {
        setResults(null);
        setRunError(error instanceof Error ? error.message : String(error));
      }
    } finally {
      if (requestId === resultsRequestRef.current) setResultsLoading(false);
    }
  }

  useEffect(() => {
    if (!run) return undefined;
    let cancelled = false;
    const runId = run.runId;

    async function refreshRun() {
      const requestId = ++stageRequestRef.current;
      try {
        const [runResponse, stagesResponse] = await Promise.all([
          fetch(`/api/runs/${runId}`, { cache: "no-store" }),
          fetch(`/api/runs/${runId}/stages`, { cache: "no-store" }),
        ]);
        if (!runResponse.ok) throw new Error(`run HTTP ${runResponse.status}`);
        if (!stagesResponse.ok) throw new Error(`stages HTTP ${stagesResponse.status}`);
        const [runPayload, stagesPayload] = await Promise.all([
          runResponse.json() as Promise<ExecutionRun>,
          stagesResponse.json() as Promise<StageSummary>,
        ]);
        if (cancelled || requestId !== stageRequestRef.current || runPayload.runId !== runId) return;
        setRun(runPayload);
        setStageSummary(stagesPayload);
      } catch (error) {
        if (!cancelled && requestId === stageRequestRef.current) {
          setRunError(error instanceof Error ? error.message : String(error));
        }
      }
    }

    void refreshRun();
    const poll = runActive ? window.setInterval(() => void refreshRun(), 1000) : null;
    return () => {
      cancelled = true;
      if (poll !== null) window.clearInterval(poll);
    };
  }, [run?.runId, runActive]);

  useEffect(() => {
    if (!run || run.status !== "complete" || results) return;
    void loadResults(run.runId);
  }, [run, results]);

  return (
    <main className="product-shell">
      <header className="topbar">
        <a className="brand" href="/">
          <span className="brand-mark"><FlaskConical size={18} aria-hidden="true" /></span>
          <span>Open Molecule Lab</span>
        </a>
        <div className="topbar-actions">
          <span className="service-state"><span className={`dot ${runnerTone}`} />{runnerLabel}</span>
          <span className="service-state muted">local-first · strict CLI</span>
        </div>
      </header>

      <section className="intro-grid">
        <div className="intro-copy">
          <p className="eyebrow">prompt-to-runspec</p>
          <h1>描述科学问题，先生成可审计的计算计划。</h1>
          <p>{scientificBoundary}</p>
        </div>
        <div className="intro-metrics" aria-label="当前请求摘要">
          <Metric label="候选池" value={formatNumber(candidatePool)} />
          <Metric label="最终选择" value={formatNumber(finalSelectionCount)} />
          <Metric label="模式" value={routeMode} tone={routeMode === "cascade" ? "warn" : "info"} />
        </div>
      </section>

      <section className="workspace">
        <form className="request-console" aria-label="研究 prompt" onSubmit={submitPrompt}>
          <div className="panel-heading">
            <Sparkles size={18} aria-hidden="true" />
            <div>
              <h2>研究 prompt</h2>
              <p>描述目标、候选规模和需要保留的科学边界。</p>
            </div>
          </div>

          <label className="field">
            <span>研究请求</span>
            <textarea
              value={prompt}
              onChange={(event) => setPrompt(event.target.value)}
              placeholder="例如：为某个靶点设计四级筛选，先检查资产，再生成可复算的运行计划。"
            />
          </label>

          <label className="field">
            <span>靶点 ID 或描述</span>
            <input value={target} onChange={(event) => setTarget(event.target.value)} placeholder="CHEMBL2051 / EGFR / target text" />
          </label>

          <div className="field-grid">
            <label className="field">
              <span>候选池规模</span>
              <input min={1} max={100000} type="number" value={candidatePool} onChange={(event) => setCandidatePool(Number(event.target.value))} />
            </label>
            <label className="field">
              <span>最终保留数</span>
              <input min={1} max={100} type="number" value={finalSelectionCount} onChange={(event) => setFinalSelectionCount(Number(event.target.value))} />
            </label>
          </div>

          <label className="field">
            <span>路由策略</span>
            <select value={routeMode} onChange={(event) => setRouteMode(event.target.value as RouteMode)}>
              <option value="auto">Auto：按受体资产决定</option>
              <option value="library">Library：四级文库筛选</option>
              <option value="cascade">Cascade：要求真实受体和 docking</option>
            </select>
          </label>

          <label className="field upload-field">
            <span>候选分子 CSV</span>
            <span className="file-picker">
              <Upload size={16} aria-hidden="true" />
              <span>{isUploading ? "正在校验 CSV" : "选择 id,smiles 文件"}</span>
              <input type="file" accept=".csv,text/csv" onChange={uploadMoleculeSet} disabled={isUploading} />
            </span>
            {moleculeSet ? (
              <small className="upload-meta">
                <FileCheck2 size={14} aria-hidden="true" />
                {moleculeSet.name} · {formatNumber(moleculeSet.nRows)} 行 · {moleculeSet.inputSha256.slice(0, 12)}
              </small>
            ) : (
              <small className="upload-meta muted">上传后才会创建不可变 MoleculeSet。</small>
            )}
          </label>

          <button className="primary-action" type="submit" disabled={!health?.ok || isSubmitting}>
            {isSubmitting ? <Loader2 size={18} aria-hidden="true" /> : <Play size={18} aria-hidden="true" />}
            {isSubmitting ? "正在生成 RunSpec" : "生成执行计划"}
          </button>

          {(serviceError || submitError || uploadError || runError) && (
            <div className="inline-error">
              <AlertTriangle size={16} aria-hidden="true" />
              <span>{submitError || uploadError || runError || serviceError}</span>
            </div>
          )}
        </form>

        <aside className="run-panel" aria-label="执行计划">
          <div className="run-panel-head">
            <span className={`state-pill ${plan ? routeTone : "idle"}`}>{plan ? plan.route.status : "waiting"}</span>
            <span>{plan?.runId || "no plan"}</span>
          </div>
          <h2>{plan ? `${plan.spec.target.id} · ${plan.route.branch}` : "等待研究 prompt"}</h2>
          <p>{plan?.route.rationale || "提交后生成真实落盘的 RunSpec、阶段计划和 manifest。"}</p>

          <div className="stage-list">
            {run ? (
              stageSummary?.runId === run.runId ? stageSummary.stages.map((stage, index) => (
                <StageHistoryRow key={stage.stage} index={index} stage={stage} />
              )) : (
                <div className="stage-loading">
                  <Loader2 size={16} aria-hidden="true" />
                  <span>读取阶段记录</span>
                </div>
              )
            ) : (plan?.stages || [
              { id: "prepare", label: "准备输入", status: "planned", reason: "等待 prompt" },
              { id: "score", label: "四级评分", status: "planned", reason: "等待 prompt" },
              { id: "dock", label: "Docking cascade", status: "planned", reason: "等待 prompt" },
              { id: "report", label: "证据报告", status: "planned", reason: "等待 prompt" },
            ]).map((stage, index) => (
              <div className={`stage-line ${toneForStatus(stage.status)}`} key={stage.id} title={stage.reason}>
                <span className="stage-index">{index + 1}</span>
                <div className="stage-copy">
                  <strong>{stage.label}</strong>
                  <small>{stage.reason}</small>
                </div>
                <span className={`state-pill ${toneForStatus(stage.status)}`}>{stage.status}</span>
                <span className="stage-attempt-meta">planned</span>
              </div>
            ))}
          </div>

          <div className="cli-preview">
            <div className="cli-preview-title">
              <TerminalSquare size={15} aria-hidden="true" />
              <span>等价 CLI</span>
            </div>
            <pre>{equivalentCli}</pre>
          </div>
          <div className="run-actions">
            <button className="primary-action" type="button" onClick={() => void launchRun()} disabled={!canRun || runActive}>
              {runActive ? <Loader2 size={18} aria-hidden="true" /> : <Play size={18} aria-hidden="true" />}
              {runActive ? "四级 CLI 运行中" : run?.status === "complete" ? "已完成，可查看结果" : "运行四级 CLI"}
            </button>
            {runActive ? (
              <button className="icon-action on-dark" type="button" title="取消运行" aria-label="取消运行" onClick={() => void cancelRun()}>
                <Square size={16} aria-hidden="true" />
              </button>
            ) : null}
            {canResume ? (
              <button className="resume-action" type="button" onClick={() => void resumeRun()} disabled={isResuming}>
                {isResuming ? <Loader2 size={16} aria-hidden="true" /> : <RotateCcw size={16} aria-hidden="true" />}
                {isResuming ? "正在验证" : "从检查点恢复"}
              </button>
            ) : null}
          </div>
          {run ? (
            <div className="run-status-line">
              <span className={`state-pill ${runTone}`}>{run.status}</span>
              <span>{run.runId}</span>
            </div>
          ) : null}
        </aside>
      </section>

      <section className="result-strip" aria-label="计划摘要">
        <Metric label="RunSpec" value={plan ? "written" : "pending"} tone={plan ? "good" : ""} />
        <Metric label="Molecule set" value={moleculeSet ? `${formatNumber(moleculeSet.nRows)} rows` : "not attached"} tone={moleculeSet ? "good" : "warn"} />
        <Metric label="Branch" value={plan?.route.branch || "unresolved"} tone={routeTone} />
        <Metric label="Run" value={run?.status || "pending"} tone={run ? runTone : ""} />
        <Metric label="Evidence" value={run?.resultSummary ? `${run.resultSummary.nRanked}/${run.resultSummary.nRows} ranked` : plan ? `${plan.bundle.files.length} files` : "pending"} />
      </section>

      <section className="results-layout">
        <section className="candidate-section" aria-label="规范化研究意图">
          <div className="section-head">
            <div>
              <p className="eyebrow">normalized intent</p>
              <h2>{plan ? "RunSpec 已生成" : "尚未生成执行计划"}</h2>
            </div>
            <span className={`state-pill ${plan ? "good" : "idle"}`}>{plan?.mode || "plan_only"}</span>
          </div>
          <p className="section-note">{plan?.spec.project.researchPrompt || prompt}</p>

          {plan ? (
            <div className="plan-facts">
              <Fact label="Target" value={plan.spec.target.id} />
              <Fact label="Requested route" value={plan.spec.target.requestedRoute} />
              <Fact label="Resolved branch" value={plan.spec.target.resolvedBranch} />
              <Fact label="Candidate count" value={formatNumber(plan.spec.moleculeSet.expectedCandidateCount)} />
              <Fact label="Final selection" value={formatNumber(plan.spec.selection.finalSelectionCount)} />
              <Fact label="Bundle" value={plan.bundle.relativeRoot} />
            </div>
          ) : (
            <EmptyState title="等待 RunSpec" detail="提交 prompt 后，这里显示规范化参数和 bundle 位置。" />
          )}
          {run?.status === "blocked" ? (
            <div className="preflight-block">
              <div className="section-head compact">
                <h3>严格预检阻断</h3>
                <span className="state-pill bad">blocked</span>
              </div>
              <PreflightList checks={run.preflight.checks} />
            </div>
          ) : null}
          {results && run ? (
            <ResultSummaryView
              results={results}
              loading={resultsLoading}
              onPage={(view, offset) => void loadResults(run.runId, view, offset)}
            />
          ) : null}
        </section>

        <aside className="evidence-panel" aria-label="后端和证据要求">
          <PanelBlock icon={<ShieldCheck size={18} aria-hidden="true" />} title="四级后端">
            <div className="evidence-list">
              {models.map((model) => (
                <article className="evidence-row" key={model.id}>
                  <span className={`state-pill ${toneForStatus(model.status)}`}>{model.status}</span>
                  <strong>{model.layer} · {model.label}</strong>
                  <p>{model.role}</p>
                </article>
              ))}
            </div>
          </PanelBlock>

          <PanelBlock icon={<Activity size={18} aria-hidden="true" />} title="执行边界">
            <div className="evidence-list">
              <article className="evidence-row">
                <span className={`state-pill ${run ? runTone : "warn"}`}>{run?.status || "plan_only"}</span>
                <strong>{run ? "运行证据" : moleculeSet ? "等待执行" : "分子集尚未接入"}</strong>
                <p>{run?.status === "complete" ? "结果来自本机严格四级 CLI；失败行保留在结果摘要中。" : plan?.executionBoundary || scientificBoundary}</p>
              </article>
              {moleculeSet ? (
                <article className="evidence-row">
                  <span className="state-pill good">sealed</span>
                  <strong>{moleculeSet.moleculeSetId}</strong>
                  <p>{formatNumber(moleculeSet.nRows)} 行 · SHA-256 {moleculeSet.inputSha256.slice(0, 16)}…</p>
                </article>
              ) : null}
              {(run?.preflight.checks || []).map((check) => (
                <article className="evidence-row" key={check.id}>
                  <span className={`state-pill ${toneForStatus(check.status)}`}>{check.status}</span>
                  <strong>{check.label}</strong>
                  <p>{check.message}</p>
                </article>
              ))}
              {(plan?.assetRequirements || []).map((requirement) => (
                <article className="evidence-row" key={requirement}>
                  <span className="state-pill idle">required</span>
                  <strong>{requirement}</strong>
                </article>
              ))}
            </div>
          </PanelBlock>
        </aside>
      </section>

      <details className="audit-drawer" open={auditOpen} onToggle={(event) => setAuditOpen(event.currentTarget.open)}>
        <summary>
          <span><ClipboardList size={16} aria-hidden="true" /> 审计材料</span>
          <ChevronDown size={16} aria-hidden="true" />
        </summary>
        <div className="audit-grid">
          <AuditBlock title="研究 prompt"><p>{prompt}</p></AuditBlock>
          <AuditBlock title="RunSpec"><pre>{JSON.stringify(plan?.spec || {}, null, 2)}</pre></AuditBlock>
          <AuditBlock title="Route"><pre>{JSON.stringify(plan?.route || {}, null, 2)}</pre></AuditBlock>
          <AuditBlock title="MoleculeSet"><pre>{JSON.stringify(moleculeSet || {}, null, 2)}</pre></AuditBlock>
          <AuditBlock title="Run"><pre>{JSON.stringify(run || {}, null, 2)}</pre></AuditBlock>
          <AuditBlock title="Stages"><pre>{JSON.stringify(stageSummary || {}, null, 2)}</pre></AuditBlock>
          <AuditBlock title="Results"><pre>{JSON.stringify(results || {}, null, 2)}</pre></AuditBlock>
          <AuditBlock title="Bundle">
            {(plan?.bundle.files || []).map((file) => <p key={file}>{file}</p>)}
            {plan?.bundle.manifestSha256 ? <code>{plan.bundle.manifestSha256}</code> : null}
          </AuditBlock>
        </div>
      </details>

      <footer className="boundary">
        <AlertTriangle size={16} aria-hidden="true" />
        <span>{run?.status === "complete" ? "本次运行显示的是本机四级 CLI 产生的可追溯结果，不等同于生物学验证或临床结论。" : scientificBoundary}</span>
      </footer>
    </main>
  );
}

function StageHistoryRow({ index, stage }: { index: number; stage: StageSummaryRow }) {
  const started = formatStageTime(stage.startedAt);
  const finished = formatStageTime(stage.finishedAt);
  const timing = started && finished ? `${started} - ${finished}` : started || finished;
  return (
    <div className={`stage-line ${toneForStatus(stage.status)}`}>
      <span className="stage-index">{index + 1}</span>
      <div className="stage-copy">
        <strong>{stageLabels[stage.stage]}</strong>
        <small>{timing || (stage.status === "waiting" ? "尚未开始" : "已持久化阶段记录")}</small>
      </div>
      <span className={`state-pill ${toneForStatus(stage.status)}`}>{stage.status}</span>
      <span className="stage-attempt-meta">
        {stage.attempts ? `${stage.attempts} 次 attempt` : "0 attempt"}
        {stage.errorCode ? <code>{stage.errorCode}</code> : null}
      </span>
    </div>
  );
}

function Metric({ label, value, tone }: { label: string; value: string; tone?: string }) {
  return (
    <div className={`metric ${tone || ""}`}>
      <span>{label}</span>
      <strong>{value}</strong>
    </div>
  );
}

function Fact({ label, value }: { label: string; value: string }) {
  return (
    <div className="plan-fact">
      <span>{label}</span>
      <strong>{value}</strong>
    </div>
  );
}

function EmptyState({ title, detail }: { title: string; detail: string }) {
  return (
    <div className="empty-state">
      <FlaskConical size={22} aria-hidden="true" />
      <strong>{title}</strong>
      <p>{detail}</p>
    </div>
  );
}

function PreflightList({ checks }: { checks: ExecutionRun["preflight"]["checks"] }) {
  return (
    <div className="preflight-list">
      {checks.filter((check) => check.required && check.status !== "passed").map((check) => (
        <div className="preflight-item" key={check.id}>
          <span className={`state-pill ${toneForStatus(check.status)}`}>{check.status}</span>
          <span><strong>{check.label}</strong>{check.message}</span>
        </div>
      ))}
    </div>
  );
}

function ResultSummaryView({
  results,
  loading,
  onPage,
}: {
  results: ResultSummary;
  loading: boolean;
  onPage: (view: ResultView, offset: number) => void;
}) {
  const pageStart = results.total ? results.offset + 1 : 0;
  const pageEnd = Math.min(results.offset + results.rows.length, results.total);
  const canGoBack = results.offset > 0;
  const canGoForward = results.offset + results.limit < results.total;
  const showsFusedScore = results.rankingScoreField === "final_score_dock";
  return (
    <section className="result-summary" aria-label="四级结果摘要">
      <div className="section-head compact">
        <div>
          <p className="eyebrow">real cli result</p>
          <h3>结果摘要</h3>
        </div>
        <span className="state-pill good">{results.nRanked}/{results.nRows} ranked</span>
      </div>
      <div className="result-toolbar">
        <div className="result-tabs" aria-label="结果视图">
          <button
            className={results.view === "ranked" ? "active" : ""}
            type="button"
            disabled={loading}
            onClick={() => onPage("ranked", 0)}
          >
            Ranked {results.nRanked}
          </button>
          <button
            className={results.view === "failed" ? "active" : ""}
            type="button"
            disabled={loading || results.nFailed === 0}
            onClick={() => onPage("failed", 0)}
          >
            Failed {results.nFailed}
          </button>
        </div>
        <div className="result-pager">
          <button
            type="button"
            title="上一页"
            aria-label="上一页"
            disabled={loading || !canGoBack}
            onClick={() => onPage(results.view, Math.max(0, results.offset - results.limit))}
          >
            <ChevronLeft size={16} aria-hidden="true" />
          </button>
          <span>{pageStart}-{pageEnd} / {results.total}</span>
          <button
            type="button"
            title="下一页"
            aria-label="下一页"
            disabled={loading || !canGoForward}
            onClick={() => onPage(results.view, results.offset + results.limit)}
          >
            <ChevronRight size={16} aria-hidden="true" />
          </button>
        </div>
      </div>
      <div className="result-table-wrap">
        <table className="result-table">
          <thead>
            <tr>
              <th>ID</th>
              <th>L1</th>
              <th>L2</th>
              <th>L3</th>
              <th>L4</th>
              <th>Base</th>
              {showsFusedScore ? <th>Fused</th> : null}
              <th>Gate</th>
            </tr>
          </thead>
          <tbody>
            {results.rows.map((row) => (
              <tr key={row.id}>
                <td><strong>{row.id}</strong><small>{row.smiles}</small></td>
                <td>{row.layer1_status || "-"}</td>
                <td>{row.layer2_status || "-"}</td>
                <td>{row.layer3_status || "-"}</td>
                <td>{row.layer4_status || "-"}</td>
                <td>{typeof row.final_score === "number" ? row.final_score.toFixed(4) : "-"}</td>
                {showsFusedScore ? (
                  <td>{typeof row.final_score_dock === "number" ? row.final_score_dock.toFixed(4) : "-"}</td>
                ) : null}
                <td>{row.gate_status || "-"}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      {results.view === "failed" && !results.rows.length ? (
        <p className="result-warning">当前运行没有失败行。</p>
      ) : null}
    </section>
  );
}

function PanelBlock({ icon, title, children }: { icon: ReactNode; title: string; children: ReactNode }) {
  return (
    <section className="panel-block">
      <div className="panel-block-title">
        {icon}
        <h3>{title}</h3>
        <ArrowRight size={14} aria-hidden="true" />
      </div>
      {children}
    </section>
  );
}

function AuditBlock({ title, children }: { title: string; children: ReactNode }) {
  return (
    <section className="audit-block">
      <h3>{title}</h3>
      {children}
    </section>
  );
}

export default App;
