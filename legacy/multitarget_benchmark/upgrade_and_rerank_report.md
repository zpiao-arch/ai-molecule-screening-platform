# L2 靶点嵌入升级 与 对接精排 — 实验报告

> 承接 `auc_root_cause_report.md`（根因：L2 靶点条件化过弱，28/164 靶点 AUC<0.5）。
> 本轮按用户优先级执行：**① 靶点嵌入升级**（打根因）→ **② 柔性/刚性对接精排**（把信号转成真实 top-K 命中）。
> 环境：离线（无网），ESM-2 缓存为空不可用；解释器 `.venv_mlx_qwen35`。

---

## ① 靶点嵌入升级：序列嵌入 L2（BLOSUM62 + 1D CNN）

### 做法
- ESM-2 离线不可用 → 用**真实蛋白序列**替代 256 维 target_text 哈希：
  `chembl_id → uniprot → chembl_37.fa 序列 → BLOSUM62 进化残基编码 → 1D CNN → 128 维靶点向量`。
- 端到端训练 `分子520维(Morgan512+8desc) + 靶点128维 → MLP → P(结合)`，79,101 样本（BindingDB+ChEMBL），850 唯一靶点，序列覆盖 **100%**。
- 面板 164 靶点 **100% 落在训练 850 靶点内**（靶点是 in-distribution）。

### 关键发现（诚实的负面结果）
| 验证口径 | 旧哈希 L2 | 新序列 L2 |
|---|---|---|
| **随机 pair 留出 AUC**（分子有泄漏） | — | 0.815 |
| **分子不相交留出 AUC**（新分子，诚实） | — | **0.801** |
| **24k×164 面板检索 AUC 中位** | **0.674** | 0.553 |
| 面板反向靶点(AUC<0.5)修复 | — | 仅 3/28 |

- 序列嵌入在**药效学测定域（assay domain）**是显著更强的结合预测器（对新分子留出 AUC 0.80）。
- 但它在**药物库检索面板**上不升反降（0.55 < 0.67）。加正则（weight_decay + 分子不相交划分）后面板仍 0.55 —— **不是过拟合分子，而是训练/部署分布错配**：
  - 训练数据 = BindingDB/ChEMBL 测定配对；面板任务 = 从 2.4 万真实药物库里把已知活性排到前面，负样本是"未标注的其他药物"（噪声负例）。
  - **升级靶点表征无法弥合这层分布鸿沟。** 低面板 AUC 的根因比"靶点表征弱"更深。

### 工程落地
- `评分_work_package/评分/target_seq_embedding.py`：`TargetSeqEncoder`（128 维）。
- `评分_work_package/评分/l2_bindingdb.py`：新增 `Layer2BindingDBSeq`（`score()` / `predict_matrix()`）。
- `scoring.py`：`--l2` 支持 `bindingdb_seq / bindingdb / deeppurpose` 三选。**两种 L2 均保留可选**；面板检索场景建议 `bindingdb`，单分子亲和力预测场景建议 `bindingdb_seq`。
- 吞吐：序列版 `predict_matrix` 全向量化，**141k pair/s**（旧版 4,574 pair/s，×30）。

---

## ② 对接精排：CHEMBL2051 神经氨酸酶 / 受体 3TI6

### 做法
- 引擎：`smina 2020.12.10`（Vina 打分）+ `obabel 3.1.0`（SMILES→3D，pH 7.4 质子化），刚性受体 3TI6，保守活性口袋 box。
- 配体集：**21 个库内已知 NA 活性 + 100 个随机诱饵**；成功对接 **活性 16 + 诱饵 90**（15 个 gen3d 失败）。

### 结果：原始 affinity 反富集 → 配体效率校正后强富集
| 排序方式 | AUC | EF@10% | R@10 | R@20 | R@30 |
|---|---|---|---|---|---|
| L2（哈希） | 0.462 | 0.0 | 0.0 | 0.0 | 0.062 |
| 原始对接 affinity | **0.394** | 0.602 | 0.062 | 0.062 | 0.062 |
| L2 + 对接 | 0.407 | 0.0 | 0.0 | 0.0 | 0.188 |
| **对接 · 配体效率 LE** | **0.756** | 0.602 | 0.062 | 0.375 | **0.688** |

- **原始 Vina affinity 反富集**（AUC 0.39<0.5）：诱饵 affinity 中位 −7.90 优于活性 −7.55。
- **诊断**：affinity 与重原子数相关 **r=−0.70** —— Vina 偏好大而亲脂分子。诱饵重原子中位 28 / cLogP 3.4，NA 活性（唾液酸模拟物，小而强极性/两性离子）重原子中位 21 / cLogP 0.6。
- **校正**：配体效率 `LE = affinity/重原子数` 消除尺寸偏差后，对接 **AUC 0.39→0.756**，**recall@30 0.06→0.688**，全面超过 L2（0.46）。物理精排确实把弱 L2 信号转成了真实 top-K 命中。

### 结论 & 后续
1. **对接精排是有效的 top-K 引擎**，但**必须做尺寸/亲脂校正**（LE 或 SILE），否则 Vina 尺寸偏好会反富集，尤其对 NA 这类小极性活性口袋。
2. LE 会过度偏向极小片段（top 仍是 HAC 14-19 的小诱饵）——生产建议用 **SILE（size-independent LE）** 或 affinity 加轻度尺寸罚项取折中。
3. 后续增益：属性匹配诱饵（DUD-E 风格，消除混杂）、柔性侧链对接、极性感知重打分（如 gnina/CNN 打分，本机 gnina 未构建成功）。

---

## 产物清单
- 序列 L2：`评分/target_seq_embedding.py`、`评分/l2_bindingdb.py::Layer2BindingDBSeq`、`评分/models/bindingdb_l2_seq/l2_seq.pt`
- 训练：`scientific_validation/multitarget_benchmark/train_l2_seq.py`（分子不相交划分 + weight_decay）
- 面板对比：`hts_closed_loop_seq.py` → `hts_closed_loop_seq_results.json`
- 对接精排：`dock_rerank_na.py` → `dock_rerank_na_results.json`（含逐配体 affinity/LE/HAC）

---

## ③ 对接精排已接入 `scoring.py` CLI（用户要求同步）

此前 CLI 的对接腿 `from flexible_docking import FlexibleDocking` 指向一个**不存在的模块**，会被 `except` 静默吞掉 → 实际上对接从未真正运行过。本轮用已验证的 `dock_rerank` 模块替换它，并把 LE 校正做成一等公民。

### 新增模块
- `评分/dock_rerank.py`：`DockingReranker`（`prep_ligand` obabel `--gen3d -p 7.4` → `dock` smina → `ligand_efficiency = affinity/HAC`），二进制自动探测（`SMINA_BIN`/`OBABEL_BIN`/`DOCK_BIN_DIR` 或 `smina-local` conda env）。

### 新增 CLI 参数
| 参数 | 作用 |
|---|---|
| `--receptor <pdbqt>` | 受体路径；提供即启用 smina 对接腿（输出 affinity 列） |
| `--box-center` / `--box-size` | 对接盒中心/尺寸 x,y,z（默认 3TI6 NA 活性口袋） |
| `--flex` / `--mcce` | 柔性残基 / 受体 pH7.4 质子化 |
| `--dock-rerank` | 启用对接精排：对所有分子跑对接 + LE 校正重排，写出 `*_reranked.csv` 并打印 Top 榜单 |
| `--dock-mode le\|raw` | 重排信号：`le`=配体效率（默认，校正尺寸偏差）/ `raw`=原始 affinity |
| `--dock-fusion 0-1` | 把 LE 归一化信号按权重融合进 `final_score`（0=仅作独立列与重排，默认 0） |

### 输出
- 主 CSV（`-o`）追加列：`docking_affinity_kcal_mol`、`heavy_atoms`、`ligand_efficiency`、`dock_rerank_rank`、`dock_le_norm`。
- 重排 CSV（`<output>_reranked.csv`）：按 `dock_rerank_rank` 排的物理精排视图。

### 端到端冒烟测试（通过）
5 分子（含 zanamivir/oseltamivir 类 NA 抑制剂 + 诱饵）对 3TI6 跑通：
- 对接腿 5/5 成功，`ligand_efficiency` 与 `dock_rerank_rank` 正确生成；
- `--dock-mode le` 与 `--dock-fusion 0.5` 两条路径均无误，重排 CSV 正常写出。
- 注意：玩具集上纯 LE 会把极小片段（如 benzene, HAC=6）排到第一——这与 ② 验证一致：LE 在**群体富集层**（AUC 0.756）有效，并非逐分子保证；生产建议加分子量窗口或用 SILE。
