import {
  Activity,
  AlertTriangle,
  ArrowRight,
  ChevronDown,
  ClipboardList,
  FlaskConical,
  Loader2,
  Play,
  ShieldCheck,
  Sparkles,
  TerminalSquare,
} from "lucide-react";
import { useEffect, useMemo, useState } from "react";
import type { FormEvent, ReactNode } from "react";

import type {
  LabHealth,
  ModelStatusRow,
  PromptError,
  PromptPlan,
  RouteMode,
} from "./types";

const scientificBoundary =
  "当前原型只生成可审计的 RunSpec 和执行计划；没有附加分子集并运行严格后端前，不构成分子评分、对接证据、药效或安全性结论。";

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

function toneForStatus(status?: string) {
  const value = String(status || "").toLowerCase();
  if (value.includes("available") || value.includes("ready") || value.includes("planned")) return "good";
  if (value.includes("blocked") || value.includes("not_found") || value.includes("error")) return "bad";
  if (value.includes("skipped") || value.includes("review")) return "warn";
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
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [auditOpen, setAuditOpen] = useState(false);

  const runnerTone = health?.ok ? "good" : health ? "bad" : "info";
  const runnerLabel = health?.ok ? "plan server online" : health ? "server unavailable" : "connecting";
  const routeTone = toneForStatus(plan?.route.status);
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
      setPlan(payload);
      setAuditOpen(false);
    } catch (error) {
      setSubmitError(error instanceof Error ? error.message : String(error));
    } finally {
      setIsSubmitting(false);
    }
  }

  return (
    <main className="product-shell">
      <header className="topbar">
        <a className="brand" href="/">
          <span className="brand-mark"><FlaskConical size={18} aria-hidden="true" /></span>
          <span>Open Molecule Lab</span>
        </a>
        <div className="topbar-actions">
          <span className="service-state"><span className={`dot ${runnerTone}`} />{runnerLabel}</span>
          <span className="service-state muted">local-first · plan-only</span>
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

          <button className="primary-action" type="submit" disabled={!health?.ok || isSubmitting}>
            {isSubmitting ? <Loader2 size={18} aria-hidden="true" /> : <Play size={18} aria-hidden="true" />}
            {isSubmitting ? "正在生成 RunSpec" : "生成执行计划"}
          </button>

          {(serviceError || submitError) && (
            <div className="inline-error">
              <AlertTriangle size={16} aria-hidden="true" />
              <span>{submitError || serviceError}</span>
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
            {(plan?.stages || [
              { id: "prepare", label: "准备输入", status: "planned", reason: "等待 prompt" },
              { id: "cache-layers", label: "缓存 L1/L3/L4", status: "planned", reason: "等待 prompt" },
              { id: "score", label: "四级评分", status: "planned", reason: "等待 prompt" },
              { id: "dock", label: "Docking cascade", status: "planned", reason: "等待 prompt" },
              { id: "report", label: "证据报告", status: "planned", reason: "等待 prompt" },
            ]).map((stage, index) => (
              <div className={`stage-line ${toneForStatus(stage.status)}`} key={stage.id} title={stage.reason}>
                <span>{index + 1}</span>
                <strong>{stage.label}</strong>
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
        </aside>
      </section>

      <section className="result-strip" aria-label="计划摘要">
        <Metric label="RunSpec" value={plan ? "written" : "pending"} tone={plan ? "good" : ""} />
        <Metric label="Molecule set" value={plan?.spec.moleculeSet.attached ? "attached" : "not attached"} tone="warn" />
        <Metric label="Branch" value={plan?.route.branch || "unresolved"} tone={routeTone} />
        <Metric label="Evidence" value={plan ? `${plan.bundle.files.length} files` : "pending"} />
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
                <span className="state-pill warn">plan_only</span>
                <strong>分子集尚未接入</strong>
                <p>{plan?.executionBoundary || scientificBoundary}</p>
              </article>
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
          <AuditBlock title="Bundle">
            {(plan?.bundle.files || []).map((file) => <p key={file}>{file}</p>)}
            {plan?.bundle.manifestSha256 ? <code>{plan.bundle.manifestSha256}</code> : null}
          </AuditBlock>
        </div>
      </details>

      <footer className="boundary">
        <AlertTriangle size={16} aria-hidden="true" />
        <span>{scientificBoundary}</span>
      </footer>
    </main>
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
