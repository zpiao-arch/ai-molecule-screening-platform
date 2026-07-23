# Open Molecule Lab Verified Stage Resume Design

日期：2026-07-23  
状态：已批准  
上位规格：`docs/superpowers/specs/2026-07-22-open-molecule-lab-design.md`

## 1. 目标

把当前一次性执行的 Open Molecule Lab worker 改成可审计的阶段编排器，使失败、取消或服务中断后的运行可以从最后一个完整且哈希验证通过的阶段边界恢复。

本增量只承诺阶段边界恢复，不承诺恢复 Python 模型、RDKit、UniMol 或 smina 子进程的内存状态。运行比较、evidence bundle 下载和 10,000 候选后台批处理继续使用相同的 RunSpec、StageAttempt 和 manifest 契约，但不在本增量内实现。

## 2. 设计原则

1. 已完成阶段不可覆盖；重试创建新的 `StageAttempt`。
2. resume 前重新验证 RunSpec、分子输入、资产 manifest、代码身份、阶段输入和阶段输出。
3. 任何指纹或输出哈希不一致都会 `blocked`，不得静默重算并称为续跑。
4. `score` 与 `dock` 必须真正分离；cascade resume 不得为了方便重新计算 L1-L4。
5. UI 只显示已持久化的阶段状态，不根据时间或日志文本猜测进度。
6. library 路由没有 `dock` StageAttempt；它是 `skipped`，不是伪造的 `complete`。

## 3. 阶段模型

一次交互运行的阶段顺序固定为：

```text
prepare -> score -> dock (cascade only) -> report
```

### 3.1 prepare

输入：attached RunSpec、密封 MoleculeSet、strict preflight、资产 manifest、代码身份。  
输出：`run-spec.json`、`preflight.json`、密封输入引用和 prepare checkpoint。  
完成条件：所有 required preflight check 通过，RunSpec 和输入文件哈希与记录一致。

### 3.2 score

输入：prepare checkpoint 和密封 `id,smiles`。  
输出：只包含 L1-L4 基础结果的 `scores.csv`，不得包含结构对接列或 `final_score_dock`。  
完成条件：结果 ID/SMILES 与输入完全一致，无重复、缺失或额外行；正式状态列与 `final_score` 契约通过。

### 3.3 dock

仅 cascade 路由创建。  
输入：score checkpoint、基础 `scores.csv`、注册受体、对接盒、smina/obabel 身份和 dock policy。  
输出：结构对接列、`final_score_dock`、可选纯物理 reranked CSV。  
完成条件：基础列逐行保持不变，至少一个结构对接行为 `ok`，融合分和排序契约通过。零成功结果为 `failed`。

### 3.4 report

输入：library 的 score 输出或 cascade 的 dock 输出。  
输出：`results/scores.csv`、`results/summary.json`、失败摘要、复算命令、完整 run manifest。  
完成条件：最终结果身份契约通过，所有公开证据不含宿主绝对路径，manifest 覆盖本 run 的所有公开文件。

## 4. StageAttempt 契约

目录结构：

```text
runs/<run_id>/
  stages/
    prepare/attempt-0001/checkpoint.json
    score/attempt-0001/checkpoint.json
    dock/attempt-0001/checkpoint.json
    report/attempt-0001/checkpoint.json
  results/
  events.jsonl
  run.json
  MANIFEST.sha256
```

`checkpoint.json` 使用 `open-molecule-lab.stage-attempt.v0.1`：

```json
{
  "schemaVersion": "open-molecule-lab.stage-attempt.v0.1",
  "runId": "run_...",
  "stage": "score",
  "attempt": 1,
  "status": "complete",
  "inputFingerprint": "sha256",
  "inputs": {
    "runSpec": "sha256",
    "moleculeSet": "sha256",
    "assetManifest": "sha256",
    "codeIdentity": "sha256",
    "prepareCheckpoint": "sha256"
  },
  "outputs": {
    "scores.csv": "sha256"
  },
  "command": {
    "executable": "OPEN_MOLECULE_PYTHON",
    "args": ["logical", "relative", "arguments"]
  },
  "startedAt": "RFC3339",
  "finishedAt": "RFC3339",
  "error": null
}
```

状态限定为 `queued`、`running`、`complete`、`failed`、`cancelled`、`blocked`。只有 `complete` checkpoint 可以被后续阶段或 resume 复用。

## 5. 输入指纹

`inputFingerprint` 是以下规范化 JSON 的 SHA-256：

```text
schemaVersion
stage
runSpecSha256
moleculeSetSha256
assetManifestSha256
codeIdentity
routeBranch
stagePolicy
upstreamCheckpointSha256
upstreamOutputSha256
```

JSON 使用 UTF-8、排序键、无多余空白。路径不进入指纹；文件由逻辑资源 ID 和内容哈希表示。

`codeIdentity` 在所有环境中都使用运行时代码清单哈希。清单包含实际参与编排与科学计算的相对路径和 SHA-256：`scoring/**/*.py`、`apps/open-molecule-lab/server/*.{mjs,py}`、`requirements.lock.txt`、`requirements-runtime.txt` 和 Open Molecule Lab `package-lock.json`。构建缓存、文档、日志、运行目录和模型资产不进入代码清单；它们分别由 source release、run manifest 和 asset manifest 追踪。

若当前目录是 Git checkout，checkpoint 可附加 commit SHA 和 dirty 状态，但它们不替代文件清单哈希。任何必需运行时文件无法读取时 strict run `blocked`；不得使用时间戳、目录 mtime 或随机值代替代码身份。

## 6. 科学 CLI 边界

现有 `scoring.py` 增加受控的 `--base-scores` 输入：

- 未提供时执行 L1-L4 并生成基础结果；
- 仅 cascade 可提供 `--base-scores`；提供后跳过 `MoleculeScorer` 和所有 L1-L4 计算；
- 加载基础结果后验证 ID/SMILES、正式状态、基础列和输入哈希；
- docking 只能追加结构列和 `final_score_dock`，不能改写基础分、backend 或模型 ID。

工作台 score 阶段以 `--mode library` 运行基础评分。原 RunSpec 的 resolved branch 仍保留 cascade；score checkpoint 明确记录其角色是 `base_four_level_score`。dock 阶段再以 `--mode cascade --base-scores <verified-score-output>` 执行。

CLI 的默认单次运行行为保持兼容；没有 `--base-scores` 时继续支持当前 `auto/library/cascade` 用法。

## 7. Resume API 与状态机

新增：

```text
POST /api/runs/:id/resume
GET  /api/runs/:id/stages
```

resume 只接受 `failed`、`cancelled` 或 `worker_interrupted` 的运行。流程：

1. 验证 run manifest 和 `specSha256`；
2. 重新运行 strict preflight；
3. 从 prepare 开始寻找连续的 complete StageAttempt；
4. 重算每个 checkpoint 的 input fingerprint 和输出哈希；
5. 遇到第一个缺失、非 complete 或不匹配阶段时停止；
6. 指纹不匹配则将运行标记 `blocked/checkpoint_mismatch`；
7. 指纹匹配则为第一个未完成阶段创建递增 attempt 并排队。

resume 保持同一 `runId` 和不可变 RunSpec。它只追加 StageAttempt、事件和终态，不覆盖历史 attempt。参数、分子、资产或代码发生变化时必须创建新运行，而不是 resume。

## 8. 取消与重启恢复

取消信号作用于当前 StageAttempt 的进程组。退出后当前 checkpoint 写为 `cancelled`，已完成 checkpoint 保持不变。

服务启动时：

- 先验证并终止属于本项目的遗留进程组；
- 将遗留 `running` attempt 写为 `failed/worker_interrupted`；
- 不自动 resume；UI 提供显式恢复命令；
- resume 仍需完整执行第 7 节验证。

## 9. Artifact Store 与事件

Artifact Store 新增原子操作：

```text
createStageAttempt(runId, stage, fingerprint)
completeStageAttempt(runId, stage, attempt, outputs)
failStageAttempt(runId, stage, attempt, error)
verifyStageAttempt(runId, stage, attempt)
nextResumableStage(runId)
```

每次状态变化先原子写 checkpoint，再追加 JSONL 事件，最后刷新 manifest。事件至少包含 `stage_attempt_started`、`stage_attempt_complete`、`stage_attempt_failed`、`stage_attempt_cancelled`、`resume_requested`、`resume_blocked` 和 `resume_started`。

## 10. UI

Run Monitor 使用真实 StageAttempt 数据：

- 每个阶段显示 `status`、attempt 次数、开始/结束时间和失败类型；
- 当前阶段显示运行态，未开始阶段显示 waiting；
- failed/cancelled 运行在 checkpoint 验证通过后显示恢复按钮；
- checkpoint mismatch 显示 blocked 和具体不匹配资源；
- library 的 dock 显示 skipped；
- 不显示估算百分比或伪造资源消耗。

恢复按钮调用 resume API，成功后复用现有轮询；运行结果只在 report complete 后显示。

## 11. 错误处理

- checkpoint JSON 缺失或格式错误：`blocked/checkpoint_invalid`；
- 输入、资产、代码或上游输出哈希变化：`blocked/checkpoint_mismatch`；
- score 输出身份不一致：`failed/result_identity_mismatch`；
- dock 修改基础列：`failed/base_score_mutated`；
- dock 零成功：`failed/docking_zero_success`；
- report manifest 或脱敏失败：`failed/evidence_incomplete`；
- 用户取消：当前 attempt `cancelled`，run `cancelled`；
- 服务中断：当前 attempt `failed/worker_interrupted`，允许显式 resume。

错误消息不得包含本机绝对路径；完整本地路径只允许存在于进程内，不进入公开 run 文件。

## 12. 验证

### Source-only 契约

- 指纹对相同规范化输入稳定；任一输入哈希变化会改变指纹；
- complete attempt 的输出被修改后 resume 返回 409 并 blocked；
- failed/cancelled attempt 不被复用；
- library 跳过 dock；cascade 必须经过 dock；
- resume 创建递增 attempt，旧 attempt 字节不变；
- run manifest、事件顺序和路径脱敏成立。

### Python 契约

- `--base-scores` 不初始化 L1-L4 scorer；
- 基础结果缺行、重复、SMILES 变化或必要列缺失时失败；
- docking 只追加允许列，基础列逐值保持；
- 原单次 library/cascade CLI 结果保持兼容。

### 真实离线验证

1. 完成两分子 library 阶段链；
2. 完成两分子 cascade 阶段链并产生真实 smina 行；
3. 在 score complete 后中断 cascade；
4. resume 后验证 score attempt 文件和哈希不变，只新增 dock/report attempt；
5. 最终结果与无中断 cascade 在记录容差内逐列一致；
6. 修改 score 输出后 resume 必须 blocked，且不启动 smina。

## 13. 非目标

- 不恢复单个模型或 docking 子进程内部状态；
- 不自动重试科学阶段；
- 不允许换资产、代码或参数后继续原 run；
- 不在本增量实现运行比较、RO-Crate、GitHub 发布或 10,000 候选调度；
- 不改变四级分数或 docking fusion 公式。

## 14. 验收标准

1. library 和 cascade 均生成完整、连续的 StageAttempt 链。
2. cascade 的 score 输出可被 dock 阶段复用，L1-L4 不重复计算。
3. 中断后可从最后一个已验证阶段恢复，历史 attempt 不被覆盖。
4. 任一指纹或输出哈希变化都会 fail closed，且不会启动下游计算。
5. UI、run.json、checkpoint 和 events.jsonl 对阶段与恢复状态的表达一致。
6. 当前 Python、source-only、worker lifecycle、真实 library/cascade 和 compact snapshot 验证继续通过。
