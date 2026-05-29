# AI 分子设计闭环 CLI 原型

这个目录提供一个轻量 CLI，用来把本地已经下载的 AI 分子生成项目和药物评估工具组织成一个闭环：

```text
靶点/口袋定义
  -> 分子生成
  -> 快速过滤
  -> docking / pose / proxy 评分
  -> 多目标排序
  -> top 分子反馈给下一轮生成
```

前三阶段可以在没有重型环境的情况下跑通工程链路；第四阶段已经接入本机可用的 RDKit，能生成真实化学描述符、Morgan 指纹相似性、Murcko scaffold、多样性选择和 docking-ready SDF；第五阶段把前四阶段资产聚合成可直接打开的本地产品仪表盘；第六阶段补质量门控、hit triage、验证队列和风险登记；第七阶段补交付包、复现手册、执行摘要和第八阶段前端产品规格。默认 `score` 仍是 `proxy`，只能用于流程演示和项目包装，不等同于真实药效或结合力证据；Vina/GNINA docking 只有在对应后端真实安装并运行后才会进入评分。

## 本地工具分工

生成侧：

- `REINVENT4`：主力 de novo design / scaffold hopping / reinforcement learning 优化器。
- `DrugEx`：多目标强化学习生成框架，可用 QSAR、约束和多目标优化。
- `GraphINVENT`：图生成模型，更适合展示“图神经网络生成分子”的路线。
- `MolBART`：SMILES/分子表示的 transformer 生成或变换思路。
- `guacamol`：分子生成基准测试和打分 benchmark。

验证侧：

- `DockStream`：把生成器和 docking 后端串起来，适合接 REINVENT。
- `AutoDock Vina`：开源 docking 基线，适合作为第一层结构打分。
- `GNINA`：CNN docking / rescoring，作为更有 AI 特征的结构评分工具。
- `PoseBusters`：检查 docked pose 是否物理合理，防止只看 docking 分数。
- `OpenFE`：对少量 top hits 做更高成本的自由能计算。
- `P2Rank`：找蛋白潜在口袋，辅助定义 docking box。

## 闭环如何成立

一个可展示、可逐步增强的闭环可以这样定义：

1. `target`：输入蛋白 PDB、已知配体或 P2Rank 预测口袋，得到 docking box。
2. `generator`：REINVENT4 / DrugEx / GraphINVENT 产生一批候选 SMILES。
3. `filter`：去重、基础合法性、类药性、合成可及性、简单性质过滤。
4. `real chemistry`：RDKit 计算真实描述符、指纹、QED、骨架，并导出 3D SDF。
5. `docking`：Vina / GNINA / DockStream 得到 docking affinity 或 CNN score。
6. `pose check`：PoseBusters 剔除不合理构象。
7. `rank`：用多目标权重合成总分，筛出 top molecules。
8. `feedback`：top molecules 作为下一轮 seed、奖励函数高分样本或微调数据。

真正“闭环”的关键不是一次性算出绝对正确的药效，而是每一轮都把验证结果反馈给生成器，让搜索空间向更符合目标的分子区域收缩。

## CLI 命令

查看闭环解释：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py explain
```

初始化项目：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py init ./ai_mol_loop/demo_project
```

第一阶段：根据疾病选择靶点：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py target-select ./ai_mol_loop/demo_project --disease 甲流 --top 5
```

甲流默认会把靶点按可用性和验证强度排序：

1. `influenza_a_h1n1_na`：A/H1N1 neuraminidase，第一推荐。
2. `influenza_a_pa_endonuclease`：PA cap-dependent endonuclease，第二推荐。
3. `influenza_a_m2`：M2 proton channel，历史对照，不建议作为主线。

乙流可以运行：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py target-select --disease 乙流 --top 5
```

从推荐靶点一键生成后续分子生成输入：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py brief-from-target ./ai_mol_loop/demo_project \
  --disease 甲流 \
  --force \
  --free-text "Use the selected influenza target as the first MVP target. Generate diverse virtual-screening candidates and keep known drugs as controls, not as claimed discoveries."
```

输出：

- `targets/target_selection.csv`
- `reports/target_selection.md`
- `briefs/target_brief.json`
- `prompts/generator_prompt.md`

第二阶段：刷新靶点证据库：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py evidence-refresh --disease influenza --retmax 3 --timeout 20
```

如果当前网络不可用，可以先生成离线证据包：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py evidence-refresh --disease influenza --offline
```

全局证据库输出：

- `targets/influenza/evidence/evidence_summary.csv`
- `targets/influenza/evidence/README.md`
- `targets/influenza/evidence/<target_id>/evidence.json`
- `targets/influenza/evidence/<target_id>/pdb_entries.csv`
- `targets/influenza/evidence/<target_id>/pubmed_articles.csv`
- `targets/influenza/evidence/<target_id>/evidence_report.md`

把第二阶段证据包整理进具体项目：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py evidence-stage2 ./ai_mol_loop/demo_project --disease influenza --top 5
```

如果希望一步完成“刷新证据包 + 写项目级第二阶段产物”，可以运行：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py evidence-stage2 ./ai_mol_loop/demo_project --disease influenza --top 5 --refresh --offline
```

项目级输出：

- `evidence/stage2_target_sources.csv`：靶点源矩阵，包含 evidence score、readiness、主 PDB、阳性对照、文献数和开放数据入口。
- `evidence/stage2_closed_loop_assets.json`：前端和后续 CLI 可直接读取的闭环资产包。
- `reports/stage2_evidence_report.md`：第二阶段靶点源报告。

证据库只保存元数据、来源链接和本地摘要，不批量保存受版权保护的全文。`evidence_score` 是“证据源可用性/闭环准备度”分数，不是活性、药效、毒理或临床有效性分数。

第三阶段：候选分子输入、生成、过滤、打分和反馈：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py stage3-screen ./ai_mol_loop/demo_project \
  --round 3 \
  --source-csv ./ai_mol_loop/demo_project/candidates/round_2_candidates.csv \
  --top 6
```

从 URL 输入候选分子：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py stage3-screen ./ai_mol_loop/demo_project \
  --round 4 \
  --source-url "https://example.org/candidates.csv" \
  --top 6
```

使用 OpenAI API 生成候选分子：

```bash
export OPENAI_API_KEY="<your-openai-api-key>"
python3 ./ai_mol_loop/ai_mol_loop.py stage3-screen ./ai_mol_loop/demo_project \
  --round 5 \
  --use-openai \
  --openai-model gpt-5.2 \
  --n 30 \
  --prompt "Generate diverse non-duplicate candidates for the selected influenza neuraminidase target. Keep controls as benchmarks only." \
  --top 8
```

也可以把 key 放在本地文件里：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py stage3-screen ./ai_mol_loop/demo_project \
  --round 5 \
  --use-openai \
  --api-key-file /path/to/openai_api_key.txt \
  --n 30 \
  --top 8
```

安全边界：

- 推荐使用 `OPENAI_API_KEY` 或 `--api-key-file`，不要把真实 `sk-...` 写进 README、JSON、CSV 或报告。
- CLI 只在本次进程中读取 key；输出文件只记录 key 来源，例如 `env:OPENAI_API_KEY`，不会保存 key 内容。
- OpenAI 生成只用于虚拟筛选候选，不生成合成步骤、剂量、临床声称或湿实验操作。

第三阶段输出：

- `stage3/round_N_raw_candidates.csv`：原始候选输入，来自 CSV/JSON/URL/OpenAI/proxy。
- `filtered/round_N_filtered.csv`：标准化、去重、RDKit/proxy 描述符、Lipinski 和风险标记。
- `candidates/round_N_candidates.csv`：通过过滤、进入 scoring 的候选。
- `scores/round_N_scores.csv`、`ranked/round_N_ranked.csv`、`feedback/round_N_feedback.json`：继续复用闭环评分和反馈。
- `stage3/round_N_stage3_assets.json`：前端可读的第三阶段资产包。
- `reports/stage3_round_N_report.md`：第三阶段报告。

第四阶段：接入本机真实化学库：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py stage4-real ./ai_mol_loop/demo_project \
  --round 4 \
  --target influenza_a_h1n1_na \
  --top 8 \
  --decoys 8 \
  --rescore
```

第四阶段当前真实执行的是 RDKit：

- 读取 `candidates/round_N_candidates.csv` 或 `--input-csv`。
- 用 RDKit 标准化 SMILES，计算 MW、Exact MW、LogP、TPSA、HBD/HBA、rotatable bonds、rings、formal charge、QED、Murcko scaffold。
- 用 Morgan fingerprint 计算候选分子与已知阳性/参考药的 Tanimoto 相似性。
- 基于 QED、控药相似性和候选间相似性生成多样化选择表。
- 为有效候选导出 `SDF`，可交给后续 OpenBabel/Meeko/Vina/GNINA。
- 生成受体结构包：记录 PDB ID、共晶/目录来源、binding-site 描述、受体准备状态和推荐 receptor preparation 步骤。
- 生成 benchmark panel：候选分子、阳性/参考控药和 decoy 放在同一个面板里，方便做 enrichment 和前端展示。
- 生成 `docking_inputs.csv` 和 `docking_plan.json`：明确哪些输入已经准备好、缺哪些工具、下一步命令是什么。
- 生成 2D 分子 PNG，供网页前端直接展示候选、控药和 decoy。
- 自动探测 `vina`、`gnina`、`obabel`、`meeko`、`posebusters` 等后端；未安装时只报告缺失，不伪造 docking 分数。
- `score` 会自动读取同轮 `stage4/round_N_real_descriptors.csv`，用真实 RDKit validity、QED、Lipinski 等字段替换对应 proxy 字段；docking/pose 分数仍只接受真实外部 CSV。
- 加 `--rescore` 后，第四阶段会一键执行 `score -> rank -> feedback`，把真实 RDKit 描述符接回闭环。
- 只有显式加 `--run-docking`，并且 receptor/ligand/backend 输入齐全时，CLI 才会尝试调用 Vina/GNINA；否则只写计划和状态。

第四阶段输出：

- `stage4/round_N_real_descriptors.csv`：RDKit 真实描述符和 scaffold。
- `stage4/round_N_similarity_to_controls.csv`：候选与 oseltamivir、zanamivir、peramivir、laninamivir 等控药的指纹相似性。
- `stage4/round_N_diverse_selection.csv`：适合进入下一步 docking 的多样化候选。
- `stage4/round_N_receptor_package.json`：靶点 PDB、binding site、受体准备状态和 receptor preparation 建议。
- `stage4/round_N_benchmark_panel.csv`：候选、阳性/参考控药和 decoy 的统一验证面板。
- `stage4/round_N_benchmark_panel.sdf`：候选、阳性/参考控药和 decoy 的统一结构面板。
- `stage4/round_N_decoys.csv`：本地 drug-like decoy 集。
- `stage4/round_N_docking_inputs.csv`：候选/控药进入 docking 前的输入状态表。
- `stage4/round_N_docking_plan.json`：Vina/GNINA/OpenBabel/Meeko/PoseBusters 的后端状态和建议命令。
- `stage4/round_N_validation_metrics.json`：panel 计数、控药相似性和 docking readiness。
- `stage4/round_N_2d/`：候选、控药和 decoy 的 2D PNG 图。
- `stage4/round_N_ligands.sdf`：RDKit 生成的 ligand SDF。
- `stage4/round_N_stage4_assets.json`：前端可读的第四阶段资产包，包含真实库/可选 docking 后端状态。
- `reports/stage4_round_N_report.md`：第四阶段报告。

如果已经有本地受体 PDB，可以加入：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py stage4-real /path/to/project \
  --round 4 \
  --target influenza_a_h1n1_na \
  --receptor-pdb /path/to/receptor.pdb \
  --docking-backend gnina \
  --run-docking
```

如果只是希望自动记录 PDB 来源但不下载受体，传 `--pdb-id 3TI6` 即可。若需要从 RCSB 拉取 PDB 文件，可显式加 `--fetch-receptor`。

第五阶段：产品化仪表盘：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py stage5-dashboard ./ai_mol_loop/demo_project \
  --round 4 \
  --title "甲流 NA AI 分子筛选闭环仪表盘"
```

第五阶段不会重新计算分子，也不会伪造 docking 分数。它会读取第一到第四阶段已有文件，生成一个静态网页产品层，用来展示：

- 靶点、证据源评分、PDB、阳性控药和靶点 readiness。
- 原始候选、过滤后候选、RDKit 有效候选、benchmark panel、advance 数量等漏斗指标。
- Top ranked candidates、真实 RDKit 描述符、候选/控药/decoy 面板。
- 2D 分子 PNG 图谱。
- Docking backend 状态、docking plan、预期外部评分 CSV。
- 反馈种子和下一轮闭环入口。
- 明确边界：computational screening only，不代表药效、安全性或临床有效性。

第五阶段输出：

- `stage5/dashboard_data.json`：前端可读的项目聚合数据包。
- `stage5/index.html`：可直接用浏览器打开的静态仪表盘。
- `stage5/styles.css`：本地样式文件，无外部依赖。
- `stage5/app.js`：本地渲染脚本，默认读取 HTML 内嵌数据，也保留 `dashboard_data.json` 入口。
- `reports/stage5_dashboard_report.md`：第五阶段交付报告。

打开方式：

```bash
open ./ai_mol_loop/demo_project/stage5/index.html
```

第六阶段：验证运营和质量门控：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py stage6-validate ./ai_mol_loop/demo_project \
  --round 4 \
  --top 8
```

第六阶段的目标不是证明药效，而是把“现在能不能继续讲这个项目”变成可检查的运营状态。它读取 Stage 5 聚合数据、ranked candidates、RDKit 描述符、benchmark panel 和 docking plan，然后产出：

- `stage6/round_N_validation_assets.json`：第六阶段总资产包，包含整体状态、质量门控摘要、下一步动作。
- `stage6/round_N_quality_gates.csv`：质量门控矩阵，覆盖 target evidence、candidate intake、RDKit validation、control panel、decoy panel、rank feedback、real docking、claim boundary。
- `stage6/round_N_hit_triage.csv`：top hit 分层表，标注 tier、优先级、RDKit 有效性、控药相似性和下一步动作。
- `stage6/round_N_assay_queue.csv`：验证队列，包含 computational docking、pose quality 和 wet-lab assay planning。wet-lab 行只作为规划，不声称已有实验结果。
- `stage6/round_N_risk_register.csv`：风险登记，例如 proxy score overclaim、docking missing、small decoy panel、wet-lab gap。
- `stage6/round_N_validation_runbook.md`：验证运营手册。
- `reports/stage6_round_N_validation_report.md`：第六阶段报告。

第六阶段会把缺失真实 docking 的状态标成 `warn`，而不是失败。含义是：产品 demo 可以继续，但不能做结构打分或真实结合力声明。

第七阶段：交付包、复现包和第八阶段规格：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py stage7-package ./ai_mol_loop/demo_project \
  --round 4 \
  --title "甲流 NA AI 分子筛选闭环交付包"
```

第七阶段把前面阶段的产物组织成可交付材料：

- `stage7/round_N_delivery_manifest.json`：交付清单，逐项检查 dashboard、Stage 6 资产、RDKit 描述符、benchmark panel、报告、复现手册、Stage 8 规格是否存在。
- `stage7/round_N_executive_summary.md`：执行摘要，用于比赛/汇报/团队 review 的项目口径。
- `stage7/round_N_reproducibility.md`：从靶点选择到 Stage 7 的复现命令。
- `stage7/round_N_investor_demo_checklist.csv`：演示检查清单，覆盖 demo story、technical depth、risk boundary、next milestone。
- `stage7/stage8_frontend_product_spec.md`：第八阶段完整前端产品说明。
- `reports/stage7_round_N_delivery_report.md`：第七阶段交付报告。

第八阶段前端页面说明已经写入：

```text
./ai_mol_loop/demo_project/stage7/stage8_frontend_product_spec.md
```

第八阶段建议做成真正的网页产品，而不是继续堆 CLI 报告。规格里定义了这些页面：

- `Project Command Center`：项目总控台，展示阶段进度、readiness、候选漏斗、下一步动作和关键文件入口。
- `Target Evidence Workspace`：靶点证据工作台，比较靶点、PDB、控药、证据分数和 assay path。
- `Candidate Funnel`：候选分子漏斗，从输入/生成到过滤、RDKit、排名、反馈。
- `Real Library Validation`：真实化学库验证页，展示 RDKit 描述符、控药相似性、benchmark panel、SDF/docking readiness。
- `Validation Operations`：第六阶段运营页，展示质量门控、hit triage、docking 队列、pose QC 队列、风险登记。
- `Delivery Room`：第七阶段交付页，展示执行摘要、manifest、复现命令、演示检查清单和导出入口。

生成结构化靶点需求和分子生成提示词：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py brief ./ai_mol_loop/demo_project \
  --target-name "SARS-CoV-2 Mpro" \
  --disease "COVID-19 antiviral discovery demo" \
  --protein "main protease" \
  --gene "nsp5" \
  --pdb-id "6LU7" \
  --pocket-source "co_crystal" \
  --reference-ligand "N3" \
  --key-residues "His41; Cys145; Gly143; His164; Glu166" \
  --pocket "catalytic cleft around the co-crystallized ligand" \
  --must-have "drug-like small molecule; hydrogen-bonding handle; synthetically accessible scaffold" \
  --avoid "reactive warheads unless explicitly justified; oversized macrocycles; PAINS-like motifs" \
  --desired-properties "MW < 500; balanced polarity; suitable for Vina/GNINA docking" \
  --free-text "Generate diverse virtual screening candidates for the Mpro catalytic pocket. Outputs are for computational ranking only." \
  --force
```

输出：

- `briefs/target_brief.json`：结构化靶点需求，给 CLI 和后续流程读取。
- `prompts/generator_prompt.md`：可喂给 LLM 或分子生成 agent 的合规提示词。

查看几类靶点需求提示词模板：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py prompt-examples ./ai_mol_loop/demo_project --write
```

跑两轮 demo 闭环：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py run-demo ./ai_mol_loop/demo_project --rounds 2 --n 24 --top 6
```

检查本地仓库和可执行程序：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py doctor ./ai_mol_loop/demo_project
```

## 产物文件

每轮都会产生这些文件：

- `candidates/round_N_candidates.csv`：候选分子。
- `scores/round_N_scores.csv`：proxy 或外部评分结果。
- `ranked/round_N_ranked.csv`：按多目标总分排序后的结果。
- `feedback/round_N_feedback.json`：反馈给下一轮的分子集合。
- `seeds/round_N_seeds.csv`：下一轮生成使用的 seed。
- `reports/round_N_summary.md`：中文项目汇报可直接引用的轮次摘要。
- `briefs/target_brief.json`：靶点、口袋、生成目标和合规边界。
- `prompts/generator_prompt.md`：面向分子生成模型的自然语言输入。
- `targets/target_selection.csv`：靶点排序表。
- `reports/target_selection.md`：靶点选择报告。

## 第一阶段靶点知识库

靶点知识库位于：

```text
./ai_mol_loop/targets/influenza
```

包含：

- `target_catalog.json`：流感靶点评分、结构、已知药物、验证路径。
- `known_drugs.csv`：已验证药物和它们在工作流里的角色。
- `pdb_structures.csv`：代表性 PDB 结构。
- `assay_plan.md`：计算验证和实验验证路径。
- `evidence_sources.json`：官方来源、开放数据库入口和 PubMed 检索式。
- `evidence/`：RCSB/PubMed 元数据证据包。

当前第一版推荐：

```text
疾病：甲流
靶点：Influenza A(H1N1) neuraminidase
推荐结构：3TI6
参考药物：oseltamivir, zanamivir, peramivir, laninamivir
验证路径：NA docking -> pose check -> neuraminidase inhibition assay
```

`M2` 保留在目录中，但标记为历史对照，因为 adamantane resistance 使它不适合作为当前主线 MVP 靶点。

第二阶段证据包已经接回 `target-select` 输出。运行靶点选择后，`target_selection.csv` 和 `target_selection.md` 会显示每个靶点已有多少 PDB 结构证据、PubMed 文献元数据、证据源评分和 readiness。

## 接入真实 docking 结果

如果已经用 Vina、GNINA、DockStream 或 PoseBusters 生成外部评分，可以整理成 CSV：

```csv
id,smiles,docking_score,pose_pass
r01_00001,c1ccccc1C(=O)O,-7.8,true
r01_00002,CCOC(=O)c1ccccc1,-6.4,false
```

然后运行：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py score /path/to/project --round 1 --external-scores /path/to/external_scores.csv
python3 ./ai_mol_loop/ai_mol_loop.py rank /path/to/project --round 1 --top 10
python3 ./ai_mol_loop/ai_mol_loop.py feedback /path/to/project --round 1 --top 10
```

`docking_score` 可以是 kcal/mol 形式的负数，也可以是 0 到 1 的归一化分数。`pose_pass` 可以填 `true/false` 或 0 到 1 的 pose 质量分。

如果已经跑过第四阶段，`score` 会默认读取同轮 RDKit 描述符：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py score /path/to/project --round 4
```

也可以显式指定描述符表：

```bash
python3 ./ai_mol_loop/ai_mol_loop.py score /path/to/project --round 4 --real-descriptors /path/to/round_4_real_descriptors.csv
```

## 项目包装口径

短期参赛或演示版本：

- 用 REINVENT4/DrugEx 作为生成器概念支撑。
- 用本 CLI 证明闭环数据流和自动筛选逻辑。
- 用 proxy 分数展示多目标筛选和反馈机制。
- 用少量真实 Vina/GNINA/PoseBusters 结果替换 top hits 的 proxy docking 字段，提高可信度。

更扎实版本：

- 选择一个明确靶点和公开 PDB 结构。
- 使用 P2Rank 或参考配体定义 pocket。
- 对每轮 top 100 做 Vina/GNINA docking。
- 对 top 20 做 PoseBusters。
- 对 top 3-5 做 OpenFE 或更高成本模拟。
- 把 top 分子作为 REINVENT4 scoring component 或 DrugEx reward 的正反馈。
