# 候选分子筛查与验证工作台 CLI 前端对接说明

本文档面向前端/全栈开发者，用来把当前 CLI 整理成可操作的网页工作台。前端不需要理解全部药物研发细节，但必须严格区分：本系统当前输出的是 **computational screening / 计算筛选与验证规划**，不是药效、临床、安全性或真实结合力证明。

## 1. 项目定位

项目名称建议：

```text
候选分子筛查与验证工作台
```

一句话定位：

```text
面向药物发现早期阶段的靶点证据整理、候选分子筛查、真实化学库校验、验证计划与材料交付工作台。
```

核心流程：

```text
靶点选择
  -> 靶点证据包
  -> 候选分子输入/整理/过滤
  -> RDKit真实化学库校验
  -> 结果看板
  -> 质量门控和验证队列
  -> 交付包和前端页面规格
```

当前 CLI 已经实现 Stage 1-7，第八阶段是待开发的网页工作台。

## 2. 目录和入口

CLI 根目录：

```bash
./ai_mol_loop
```

主程序：

```bash
./ai_mol_loop/ai_mol_loop.py
```

推荐前端调用入口：

```bash
./ai_mol_loop/ai-mol-loop
```

示例项目目录：

```bash
./ai_mol_loop/demo_project
```

前端开发时建议把每个项目视为一个 workspace。所有阶段产物都写入项目目录。

## 3. 运行方式

查看总帮助：

```bash
ai_mol_loop/ai-mol-loop --help
```

查看某个命令帮助：

```bash
ai_mol_loop/ai-mol-loop stage4-real --help
ai_mol_loop/ai-mol-loop stage7-package --help
```

前端/后端调用命令时应使用：

```bash
cwd=.
```

命令成功标准：

- 进程 exit code 为 `0`
- 目标输出文件存在
- JSON 文件可解析
- CSV 表头符合对应阶段契约

命令失败标准：

- exit code 非 `0`
- stderr/stdout 中出现 `Missing config`、`No candidates found`、`Target not found`、`OpenAI API request failed` 等错误
- 目标输出文件未产生

## 4. 阶段总览

| Stage | CLI 命令 | 前端页面 | 主要作用 |
|---:|---|---|---|
| 0 | `init` | 项目登记 | 初始化项目目录和默认配置 |
| 1 | `target-select` / `brief-from-target` / `brief` | 靶点证据页 | 选择靶点并整理结构化靶点说明 |
| 2 | `evidence-refresh` / `evidence-stage2` | 靶点证据页 | 整理靶点证据源、PDB、控药、文献元数据 |
| 3 | `stage3-screen` | 候选漏斗 | 候选分子输入、整理、过滤、打分和反馈 |
| 4 | `stage4-real` | 真实库校验 | RDKit 描述符、控药相似性、SDF、对照面板、对接计划 |
| 5 | `stage5-dashboard` | 项目总览 | 生成前端聚合 JSON 和本地静态结果看板 |
| 6 | `stage6-validate` | 验证计划 | 质量门控、候选复核、验证队列、风险登记 |
| 7 | `stage7-package` | 交付材料 | 交付清单、执行摘要、复现手册、演示清单、Stage 8 规格 |
| 8 | 暂无实现 | 正式网页工作台 | 待开发的正式前端页面 |

## 5. 推荐前端页面

### 5.1 项目总览

对应数据：

- `stage5/dashboard_data.json`
- `stage6/round_N_validation_assets.json`
- `stage7/round_N_delivery_manifest.json`

展示内容：

- 项目名称、项目路径、当前轮次
- 目标靶点、疾病背景、PDB、参考配体
- 阶段状态标签
- 候选分子漏斗指标
- Stage 6 下一步动作
- Stage 7 交付状态

主要按钮：

- 运行 Stage 5：`stage5-dashboard`
- 运行 Stage 6：`stage6-validate`
- 运行 Stage 7：`stage7-package`
- 打开本地结果看板：`stage5/index.html`

### 5.2 靶点证据页

对应数据：

- `targets/target_selection.csv`
- `briefs/target_brief.json`
- `prompts/generator_prompt.md`
- `evidence/stage2_target_sources.csv`
- `evidence/stage2_closed_loop_assets.json`
- `reports/target_selection.md`
- `reports/stage2_evidence_report.md`

展示内容：

- 靶点候选列表
- 证据分
- 准备状态
- PDB、参考配体、控药
- 实验路径
- 生成模型可用输入

主要按钮：

- 选择靶点：`target-select`
- 由靶点目录整理说明：`brief-from-target`
- 手工创建靶点说明：`brief`
- 整理项目级证据包：`evidence-stage2`

### 5.3 候选漏斗

对应数据：

- `stage3/round_N_raw_candidates.csv`
- `filtered/round_N_filtered.csv`
- `candidates/round_N_candidates.csv`
- `scores/round_N_scores.csv`
- `ranked/round_N_ranked.csv`
- `feedback/round_N_feedback.json`
- `seeds/round_N_seeds.csv`
- `stage3/round_N_stage3_assets.json`
- `reports/stage3_round_N_report.md`

展示内容：

- 原始候选数量
- 通过过滤数量
- 排序后候选
- 进入下一轮 / 暂缓决策
- 反馈种子

主要按钮：

- 从 CSV/JSON/URL 导入候选
- 使用 OpenAI API 生成候选
- 运行筛选：`stage3-screen`
- 单独运行 `score` / `rank` / `feedback`

### 5.4 真实库校验

对应数据：

- `stage4/round_N_real_descriptors.csv`
- `stage4/round_N_similarity_to_controls.csv`
- `stage4/round_N_diverse_selection.csv`
- `stage4/round_N_benchmark_panel.csv`
- `stage4/round_N_benchmark_panel.sdf`
- `stage4/round_N_decoys.csv`
- `stage4/round_N_docking_inputs.csv`
- `stage4/round_N_docking_plan.json`
- `stage4/round_N_validation_metrics.json`
- `stage4/round_N_2d/`
- `stage4/round_N_stage4_assets.json`
- `reports/stage4_round_N_report.md`

展示内容：

- RDKit描述符表
- 分子 2D 图片
- 候选/控药/Decoy面板
- 对接后端状态
- 受体准备状态
- 预期外部分数CSV

主要按钮：

- 运行 Stage 4：`stage4-real`
- 尝试真实对接：`stage4-real --run-docking`
- 导入外部对接CSV后重跑评分：`score --external-scores`

### 5.5 验证计划

对应数据：

- `stage6/round_N_validation_assets.json`
- `stage6/round_N_quality_gates.csv`
- `stage6/round_N_hit_triage.csv`
- `stage6/round_N_assay_queue.csv`
- `stage6/round_N_risk_register.csv`
- `stage6/round_N_validation_runbook.md`
- `reports/stage6_round_N_validation_report.md`

展示内容：

- 通过/警示/失败质量门控
- 重点候选复核
- 计算对接队列
- Pose质量检查队列
- 湿实验计划行
- 风险登记

主要按钮：

- 运行 Stage 6：`stage6-validate`
- 打开验证报告
- 根据验证队列提醒用户下一步需要真实对接或pose检查

### 5.6 交付材料

对应数据：

- `stage7/round_N_delivery_manifest.json`
- `stage7/round_N_executive_summary.md`
- `stage7/round_N_reproducibility.md`
- `stage7/round_N_investor_demo_checklist.csv`
- `stage7/stage8_frontend_product_spec.md`
- `reports/stage7_round_N_delivery_report.md`

展示内容：

- 交付清单
- 缺失文件检查
- 项目执行摘要
- 复现命令
- 演示检查清单
- Stage 8 前端规格

主要按钮：

- 运行 Stage 7：`stage7-package`
- 打开执行摘要
- 打开复现手册
- 打开第八阶段规格

## 6. 常用命令详解

### 6.1 初始化项目

```bash
ai_mol_loop/ai-mol-loop init <project> [--force]
```

示例：

```bash
ai_mol_loop/ai-mol-loop init ai_mol_loop/demo_project
```

输出：

- `config.json`
- `seeds/round_0_seeds.csv`
- `README.md`

前端使用：

- 新建项目按钮调用此命令。
- 如果 `config.json` 已存在，不应默认传 `--force`，除非用户确认覆盖。

### 6.2 靶点选择

```bash
ai_mol_loop/ai-mol-loop target-select [project] --disease <query> --top 5 [--target <filter>] [--catalog <path>]
```

示例：

```bash
ai_mol_loop/ai-mol-loop target-select ai_mol_loop/demo_project --disease 甲流 --top 5
```

输出：

- `targets/target_selection.csv`
- `reports/target_selection.md`

CSV 关键字段：

- `rank`
- `target_id`
- `display_name`
- `score`
- `recommendation`
- `pdb_evidence_count`
- `pubmed_evidence_count`
- `evidence_score`
- `readiness`
- `known_drugs`
- `recommended_pdb`
- `recommendation_reason`

前端使用：

- 用表格展示靶点排序。
- 用户点击某个 target 后，调用 `brief-from-target --target <target_id>`。

### 6.3 由靶点生成 brief 和 prompt

```bash
ai_mol_loop/ai-mol-loop brief-from-target <project> \
  --disease <query> \
  [--target <target_id>] \
  [--free-text <requirement>] \
  [--max-heavy-atoms 55] \
  [--max-molecular-weight 550] \
  [--force]
```

示例：

```bash
ai_mol_loop/ai-mol-loop brief-from-target ai_mol_loop/demo_project \
  --disease 甲流 \
  --target influenza_a_h1n1_na \
  --free-text "Generate diverse virtual-screening candidates and keep known drugs as controls only." \
  --force
```

输出：

- `briefs/target_brief.json`
- `prompts/generator_prompt.md`
- 更新 `config.json`

前端使用：

- `target_brief.json` 用于展示靶点、口袋、设计约束。
- `generator_prompt.md` 可给 LLM 或生成模块展示/复制。

### 6.4 手工创建 brief

```bash
ai_mol_loop/ai-mol-loop brief <project> \
  --target-name <name> \
  --disease <disease> \
  --protein <protein> \
  --gene <gene> \
  --pdb-id <pdb> \
  --pocket-source <source> \
  --reference-ligand <ligand> \
  --key-residues "A; B; C" \
  --pocket <description> \
  --must-have "..." \
  --avoid "..." \
  --desired-properties "..." \
  --free-text "..." \
  --force
```

适用场景：

- 前端支持用户手工输入非流感靶点。
- 用户已有 PDB、口袋、参考配体、设计需求。

输出：

- `briefs/target_brief.json`
- `prompts/generator_prompt.md`
- 更新 `config.json`

### 6.5 刷新全局证据包

```bash
ai_mol_loop/ai-mol-loop evidence-refresh \
  --disease influenza \
  [--target <target_id>] \
  [--retmax 5] \
  [--timeout 20] \
  [--offline]
```

输出默认写入：

- `targets/influenza/evidence/evidence_summary.csv`
- `targets/influenza/evidence/<target_id>/evidence.json`
- `targets/influenza/evidence/<target_id>/pdb_entries.csv`
- `targets/influenza/evidence/<target_id>/pubmed_articles.csv`
- `targets/influenza/evidence/<target_id>/evidence_report.md`

前端使用：

- 可以做成“刷新证据库”按钮。
- 网络不可用时默认用 `--offline`，避免前端卡死。

### 6.6 生成项目级 Stage 2 证据资产

```bash
ai_mol_loop/ai-mol-loop evidence-stage2 <project> \
  --disease influenza \
  --top 5 \
  [--target <target_id>] \
  [--refresh] \
  [--offline]
```

输出：

- `evidence/stage2_target_sources.csv`
- `evidence/stage2_closed_loop_assets.json`
- `reports/stage2_evidence_report.md`

前端使用：

- 靶点证据页的核心数据源。
- 如果用户换了疾病/靶点，应重跑此命令。

### 6.7 候选分子输入、生成、过滤、评分和反馈

```bash
ai_mol_loop/ai-mol-loop stage3-screen <project> \
  --round N \
  [--n 30] \
  [--top 10] \
  [--source-csv <path>] \
  [--source-json <path>] \
  [--source-url <url>] \
  [--context-url <url>] \
  [--use-openai] \
  [--openai-model gpt-5.2] \
  [--api-key-file <path>] \
  [--api-key-env OPENAI_API_KEY] \
  [--prompt "..."] \
  [--external-scores <csv>] \
  [--no-score]
```

CSV 输入要求：

- 必须至少有 `smiles` 字段。
- 可选字段：`id`、`rationale`、`expected_interaction`、`design_family`、`risk_note`。

输出：

- `stage3/round_N_raw_candidates.csv`
- `filtered/round_N_filtered.csv`
- `candidates/round_N_candidates.csv`
- `scores/round_N_scores.csv`
- `ranked/round_N_ranked.csv`
- `feedback/round_N_feedback.json`
- `seeds/round_N_seeds.csv`
- `stage3/round_N_stage3_assets.json`
- `reports/stage3_round_N_report.md`

前端使用：

- 上传候选 CSV/JSON 后调用此命令。
- URL 输入候选时调用 `--source-url`。
- OpenAI 生成候选时，前端只传 key 来源，不要把 key 写入项目文件。

OpenAI key 安全要求：

- 推荐后端从环境变量读取 `OPENAI_API_KEY`。
- 或用户上传临时 key 文件，用 `--api-key-file`。
- 不要把 `sk-...` 保存进 JSON、CSV、日志、报告或前端 localStorage。

### 6.8 单独生成候选

```bash
ai_mol_loop/ai-mol-loop generate <project> --round N --n 50
```

或导入 CSV：

```bash
ai_mol_loop/ai-mol-loop generate <project> --round N --source-csv <path>
```

输出：

- `candidates/round_N_candidates.csv`

前端使用：

- 适合作为简单 demo，不如 `stage3-screen` 完整。

### 6.9 单独打分

```bash
ai_mol_loop/ai-mol-loop score <project> --round N [--external-scores <csv>] [--real-descriptors <csv>]
```

外部 docking CSV 推荐字段：

```csv
id,smiles,docking_score,pose_pass,backend,receptor,notes
mol_001,CCOC(=O)c1ccccc1,-7.8,true,vina,receptor.pdbqt,ok
```

输出：

- `scores/round_N_scores.csv`

关键字段：

- `id`
- `smiles`
- `validity_proxy`
- `qed_proxy`
- `sa_proxy`
- `lipinski_proxy`
- `docking_proxy`
- `pose_proxy`
- `novelty_proxy`
- `raw_docking_kcal_mol`
- `total_proxy`
- `score_source`

前端注意：

- `score_source` 包含 `proxy` 时，不得展示为真实实验结果。
- `score_source` 包含 `external_docking` 才能说明导入过外部 docking 分数。

### 6.10 单独排序

```bash
ai_mol_loop/ai-mol-loop rank <project> --round N --top 10
```

输出：

- `ranked/round_N_ranked.csv`
- `reports/round_N_summary.md`

关键字段：

- `rank`
- `id`
- `smiles`
- `total_proxy`
- `decision`

### 6.11 单独反馈

```bash
ai_mol_loop/ai-mol-loop feedback <project> --round N --top 10
```

输出：

- `seeds/round_N_seeds.csv`
- `feedback/round_N_feedback.json`

前端使用：

- 展示下一轮 seed 分子。
- 作为“开始下一轮”的输入来源。

### 6.12 Stage 4：真实化学库校验

```bash
ai_mol_loop/ai-mol-loop stage4-real <project> \
  --round N \
  --target <target_id> \
  --top 8 \
  --decoys 8 \
  [--input-csv <path>] \
  [--controls-csv <path>] \
  [--receptor-pdb <path>] \
  [--pdb-id <pdb_id>] \
  [--fetch-receptor] \
  [--docking-backend auto|vina|gnina] \
  [--run-docking] \
  [--rescore]
```

当前真实执行内容：

- RDKit SMILES 解析
- canonical SMILES
- 分子式、MW、LogP、TPSA、HBD/HBA、rotatable bonds、rings、formal charge、QED
- Lipinski violations
- Murcko scaffold
- Morgan fingerprint similarity
- 控药/decoy benchmark panel
- SDF 导出
- 2D PNG 分子图
- 对接计划

不会自动伪造内容：

- 不会伪造 Vina/GNINA 对接分数
- 不会伪造 PoseBusters 结果
- 不会伪造实验 assay 结果

输出：

- `stage4/round_N_real_descriptors.csv`
- `stage4/round_N_similarity_to_controls.csv`
- `stage4/round_N_diverse_selection.csv`
- `stage4/round_N_receptor_package.json`
- `stage4/round_N_benchmark_panel.csv`
- `stage4/round_N_benchmark_panel.sdf`
- `stage4/round_N_decoys.csv`
- `stage4/round_N_docking_inputs.csv`
- `stage4/round_N_docking_plan.json`
- `stage4/round_N_validation_metrics.json`
- `stage4/round_N_2d/`
- `stage4/round_N_ligands.sdf`
- `stage4/round_N_stage4_assets.json`
- `reports/stage4_round_N_report.md`

前端使用：

- 真实库校验页主要读取这些文件。
- 分子图库使用 `stage4/round_N_2d/*.png`。
- 对接状态读 `stage4/round_N_docking_plan.json` 和 `stage4/round_N_validation_metrics.json`。

### 6.13 Stage 5：前端聚合结果看板

```bash
ai_mol_loop/ai-mol-loop stage5-dashboard <project> \
  --round N \
  [--title "..."]
```

输出：

- `stage5/dashboard_data.json`
- `stage5/index.html`
- `stage5/styles.css`
- `stage5/app.js`
- `reports/stage5_dashboard_report.md`

`stage5/dashboard_data.json` 是前端最重要的数据入口。

顶层字段：

```json
{
  "schema_version": "0.1",
  "generated_at": "...",
  "stage": 5,
  "round": 4,
  "title": "...",
  "project": {},
  "target": {},
  "readiness": {},
  "metrics": {},
  "files": {},
  "tables": {},
  "molecule_images": [],
  "boundary": []
}
```

关键 `metrics`：

- `raw_candidates`
- `filtered_candidates`
- `candidate_rows`
- `scored_candidates`
- `ranked_candidates`
- `advanced`
- `valid_rdkit`
- `invalid_rdkit`
- `controls`
- `decoys`
- `benchmark_panel`
- `molecule_images`

关键 `readiness`：

- `target_evidence`
- `candidate_intake`
- `rdkit_validation`
- `control_panel`
- `decoy_panel`
- `docking`
- `feedback`

关键 `tables`：

- `stage2_targets`
- `ranked_top`
- `real_descriptors`
- `benchmark_panel`
- `control_similarity_top`
- `seeds`
- `validation_metrics`
- `docking_plan`
- `feedback`

前端使用：

- 项目总览可以优先只读这个 JSON。
- 如果要做更详细的子页面，再读各阶段 CSV/JSON。

### 6.14 Stage 6：质量门控与验证运营

```bash
ai_mol_loop/ai-mol-loop stage6-validate <project> \
  --round N \
  --top 8
```

输出：

- `stage6/round_N_validation_assets.json`
- `stage6/round_N_quality_gates.csv`
- `stage6/round_N_hit_triage.csv`
- `stage6/round_N_assay_queue.csv`
- `stage6/round_N_risk_register.csv`
- `stage6/round_N_validation_runbook.md`
- `reports/stage6_round_N_validation_report.md`

`validation_assets.json` 顶层字段：

```json
{
  "stage": 6,
  "round": 4,
  "overall_status": "computational_demo_ready_real_docking_missing",
  "docking_status": "skipped",
  "quality_gate_summary": {},
  "triage_count": 6,
  "queue_count": 13,
  "risk_count": 4,
  "files": {},
  "next_actions": [],
  "boundary": []
}
```

`quality_gates.csv` 关键字段：

- `gate_id`
- `gate_name`
- `status`
- `evidence`
- `required_next_step`

`status` 可取：

- `pass`
- `warn`
- `fail`

`hit_triage.csv` 关键字段：

- `rank`
- `id`
- `smiles`
- `total_score`
- `qed`
- `mw`
- `lipinski_violations`
- `rdkit_valid`
- `nearest_control`
- `control_similarity`
- `validation_tier`
- `priority`
- `next_action`
- `claim_allowed`

`assay_queue.csv` 关键字段：

- `queue_type`
- `priority`
- `id`
- `smiles`
- `target_id`
- `input_or_output`
- `acceptance_criterion`
- `owner_note`

前端使用：

- 验证计划页面读取这些文件。
- 如果 `real_docking` gate 是 `warn`，UI 要明确显示：真实对接尚未完成。

### 6.15 Stage 7：交付包和前端规格

```bash
ai_mol_loop/ai-mol-loop stage7-package <project> \
  --round N \
  [--title "..."]
```

输出：

- `stage7/round_N_delivery_manifest.json`
- `stage7/round_N_executive_summary.md`
- `stage7/round_N_reproducibility.md`
- `stage7/round_N_investor_demo_checklist.csv`
- `stage7/stage8_frontend_product_spec.md`
- `reports/stage7_round_N_delivery_report.md`

`delivery_manifest.json` 顶层字段：

```json
{
  "stage": 7,
  "round": 4,
  "project": "...",
  "delivery_status": "ready_for_demo_package",
  "stage6_status": "...",
  "metrics_snapshot": {},
  "readiness_snapshot": {},
  "deliverables": {},
  "claims_boundary": []
}
```

`deliverables` 每项格式：

```json
{
  "description": "...",
  "path": "...",
  "relative_path": "...",
  "exists": true
}
```

前端使用：

- 交付材料页面读取 manifest。
- 对 `exists=false` 的交付物显示红色缺失状态。
- Stage 8 前端开发者应直接阅读 `stage7/stage8_frontend_product_spec.md`。

### 6.16 诊断本地环境

```bash
ai_mol_loop/ai-mol-loop doctor [project]
```

输出：

- 本地仓库路径是否存在
- `vina`、`gnina`、`obabel`、`posebusters` 等可执行程序是否可用
- RDKit、numpy、pandas 等 Python 模块状态

前端使用：

- 可以做成环境检查页面。
- 如果对接后端不可用，Stage 4/6 仍可运行，但不能展示真实对接分数。

## 7. 推荐完整运行链路

从空项目开始：

```bash
PROJECT=./ai_mol_loop/demo_project

ai_mol_loop/ai-mol-loop init "$PROJECT"

ai_mol_loop/ai-mol-loop target-select "$PROJECT" --disease 甲流 --top 5

ai_mol_loop/ai-mol-loop brief-from-target "$PROJECT" \
  --disease 甲流 \
  --target influenza_a_h1n1_na \
  --force \
  --free-text "Generate diverse virtual-screening candidates and keep known drugs as controls only."

ai_mol_loop/ai-mol-loop evidence-stage2 "$PROJECT" --disease influenza --top 5

ai_mol_loop/ai-mol-loop stage3-screen "$PROJECT" --round 4 --n 30 --top 8

ai_mol_loop/ai-mol-loop stage4-real "$PROJECT" \
  --round 4 \
  --target influenza_a_h1n1_na \
  --top 8 \
  --decoys 8 \
  --rescore

ai_mol_loop/ai-mol-loop stage5-dashboard "$PROJECT" \
  --round 4 \
  --title "候选分子筛查与验证工作台"

ai_mol_loop/ai-mol-loop stage6-validate "$PROJECT" --round 4 --top 8

ai_mol_loop/ai-mol-loop stage7-package "$PROJECT" \
  --round 4 \
  --title "候选分子筛查与验证材料包"
```

前端应支持按阶段运行，也应支持“一键运行到 Stage 7”的处理流程。

## 8. 推荐本地 API 封装

前端不要直接在浏览器里执行 shell 命令。建议做一个本地后端，例如 Python FastAPI、Node.js Express、Tauri command 或 Electron main process。

推荐 API：

```text
GET  /api/projects
POST /api/projects
GET  /api/projects/:projectId/status
GET  /api/projects/:projectId/dashboard?round=N
POST /api/projects/:projectId/commands/init
POST /api/projects/:projectId/commands/target-select
POST /api/projects/:projectId/commands/brief-from-target
POST /api/projects/:projectId/commands/evidence-stage2
POST /api/projects/:projectId/commands/stage3-screen
POST /api/projects/:projectId/commands/stage4-real
POST /api/projects/:projectId/commands/stage5-dashboard
POST /api/projects/:projectId/commands/stage6-validate
POST /api/projects/:projectId/commands/stage7-package
GET  /api/projects/:projectId/files?path=...
```

命令请求示例：

```json
{
  "round": 4,
  "target": "influenza_a_h1n1_na",
  "top": 8,
  "decoys": 8,
  "rescore": true
}
```

命令响应示例：

```json
{
  "ok": true,
  "exitCode": 0,
  "command": "ai_mol_loop/ai-mol-loop stage4-real ...",
  "stdout": "...",
  "stderr": "",
  "outputs": [
    "stage4/round_4_real_descriptors.csv",
    "stage4/round_4_stage4_assets.json"
  ]
}
```

长任务处理：

- 后端创建 job id
- 前端轮询 `GET /api/jobs/:id`
- 后端持续保存 stdout/stderr
- 任务完成后重新读取产物 JSON/CSV

## 9. 前端状态模型

建议前端维护：

```ts
type ProjectState = {
  projectPath: string;
  currentRound: number;
  stage5?: ResultBoardData;
  stage6?: ValidationAssets;
  stage7?: DeliveryManifest;
  files: Record<string, FileStatus>;
  jobs: JobStatus[];
};
```

关键状态来源：

- Project status：检查 `config.json`
- Target status：检查 `briefs/target_brief.json` 和 `evidence/stage2_target_sources.csv`
- Candidate status：检查 `candidates/round_N_candidates.csv`
- Scoring status：检查 `ranked/round_N_ranked.csv`
- RDKit status：检查 `stage4/round_N_real_descriptors.csv`
- Result board status：检查 `stage5/dashboard_data.json`
- Validation status：检查 `stage6/round_N_validation_assets.json`
- Delivery status：检查 `stage7/round_N_delivery_manifest.json`

## 10. UI 文案边界

必须使用的安全文案：

```text
本工作台当前输出为计算筛选、候选优先级排序和验证规划，不代表真实药效、安全性、临床有效性或实验活性。
```

禁止文案：

- “发现了有效药物”
- “验证了药效”
- “证明可治疗”
- “真实结合力已确认”
- “临床有效”

推荐文案：

- “计算筛选候选”
- “RDKit 有效分子”
- “待 docking 验证”
- “待 pose QC”
- “已生成验证队列”
- “已生成交付包”
- “已知药物作为阳性/参考控药”

## 11. 前端优先级建议

第一版 MVP：

1. 项目总览：读取 `stage5/dashboard_data.json`
2. 验证计划：读取 Stage 6 文件
3. 交付材料：读取 Stage 7 manifest
4. 候选漏斗：展示排序候选和分子图
5. 靶点证据页：展示证据和靶点说明
6. 命令执行后端：封装 CLI

不要第一版就做复杂的三维结构浏览器。当前项目更需要把“靶点-候选-验证-交付”讲清楚。

## 12. 当前 demo 项目真实状态

当前 demo 项目已经生成到 Stage 7。

关键状态：

```text
Stage 5: 结果看板已生成
Stage 6: computational_demo_ready_real_docking_missing
Stage 7: ready_for_demo_package
对接: skipped
```

含义：

- 项目可以作为计算筛选演示材料展示。
- RDKit 层是真实本地化学库计算。
- 真实对接尚未运行，不能声称真实对接验证完成。
- 湿实验只在 Stage 6 中作为计划队列，不是已有实验结果。
