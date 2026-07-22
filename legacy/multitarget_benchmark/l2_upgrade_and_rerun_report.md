# L2 升级接入 + 24k药库 × 164靶点 重跑验证报告

> 日期: 2026-07-14
> 目标: 将失效的 L2 (DeepPurpose) 替换为真实可运行的靶点感知模型并接入 MoleculeScorer (权重 0.50),
> 并用 Boltz-2 多构象 + MCCE 喂柔性对接替换单构象刚性对接, 在同一组 (24k 候选药库 × 164 靶点面板) 上重跑验证 AUC 回到 0.5 以上。

---

## 1. 做了什么

### 1.1 替换了"那个 AUC 无效的模型"
- **原判**: `MoleculeScorer` 的 L2 槽位是 `Layer2DeepPurpose`, 因 DeepPurpose 未安装, 运行时返回 `docking_normalized=0.0` → 占 0.50 权重的 L2 实为空, 闭环"瞎选"。这正是对话中用户指的"之前那个 AUC 无效的模型"。
- **替换**: 新建 `评分_work_package/评分/l2_bindingdb.py` → `Layer2BindingDB`, 用本地 BindingDB 已对齐的 5 万条 `(药物SMILES, 靶点文本, 结合标签)` 训练靶点感知分类器, **保存权重** (上一轮只存了 AUC 数字、漏存权重, 本次已修正)。
  - 特征 (776 维, 自洽可复现): Morgan r2 fpSize=512 + 8 个 RDKit 描述符 + 256 维靶点文本哈希 (sklearn `FeatureHasher`)。
  - 模型: 逻辑回归 (留出 AUC 0.851) / MLP(256,128) (留出 AUC **0.935**); 推理用 MLP。
  - 训练/推理共用同一 `BindingDBFeature`, 保证特征空间一致。
- **接入**: `MoleculeScorer.__init__` 默认 `l2_method="bindingdb"`, 加载 `Layer2BindingDB`, **保留 0.50 权重** (line 698: `dock_norm*0.50`); `score_one/score_batch` 新增 `target_text` 透传; `Layer2DeepPurpose` 保留为 legacy。冒烟测试: CA 抑制剂 P(结合)=0.98, 乙醇=0.001 → 靶点感知生效。

### 1.2 重跑协议 (同一组面板, 真实模型, 不偷看 actives)
- 候选库: 24,268 个本地药 (≥ 用户要求的 1 万)。
- 面板: 164 个 ChEMBL 靶点 (≥10 正 + ≥10 负)。
- 靶点文本对齐: 面板 `name` → BindingDB `target_text` (全串 + 尾部规范名两遍匹配), **116/164 (71%)** 映射到模型训练时见过的靶点。
- 负类: 面板已知阴性 + 500 个共享随机库药 (排除面板已知药, 避免泄漏)。
- 指标: 逐靶点算 ROC-AUC (活性 vs 阴性+decoy), 用 (a) 真实 L2 分数, (b) 集成 `final_score = 0.20·L1 + 0.50·L2 + 0.20·L3`, 及 (c) 无 L2 的 L1 基线作对照。

### 1.3 柔性对接 M2 管线 (Boltz-2 + MCCE + smina --flex)
- 新建 `scientific_validation/multitarget_benchmark/flexible_docking.py` → `FlexibleDocking` 编排类:
  - `boltz_multiconf()` 钩子 (多构象受体), `mcce_protonate()` 钩子 (质子化态), `dock()` 核心用 **smina `--flexres`** 柔性对接。
- **环境现实 (必须如实说明)**:
  - ✅ **smina `--flex`**: 可用 (smina-local 环境, Open Babel 3.1.0)。
  - ⚠️ **Boltz-2**: 当前环境 `import boltz` 深层依赖缺 `einops` 直接 `ModuleNotFoundError`, 且真实蛋白结构预测需 GPU, CPU 不可行 → 钩子保留, 环境就绪后接入。
  - ⚠️ **MCCE (质子化)**: 仓库内 `mcce_lite/` 是 ML DPO 训练/评估框架, **非质子化工具**; 真实 MCCE 二进制未提供 → 钩子保留, 提供后接入。目前用 `obabel --addh --pH 7.4` 作代理加氢。
- 因此**本环境可立即运行的是 smina `--flex` 柔性对接腿** (核心升级), 已在 NA(3TI6) 上演示 (刚性 vs 柔性对照, 见第 3 节)。Boltz/MCCE 上游为可选、缺失时优雅降级。

---

## 2. 重跑结果: AUC 是否回到 0.5 以上?

**是 — 中位数与多数靶点击穿 0.5, 证实了"替换无效模型"有效。**

| 判别器 | 靶点数 | 中位数 AUC | 均值 AUC | >0.5 占比 | >0.7 占比 |
|---|---|---|---|---|---|
| 无 L2 基线 (L1/QED, "瞎选") | 116 | **0.493** | 0.495 | 49.1% | 10.3% |
| 真实 L2 (集成权重 0.50) | 116 | **0.596** | 0.582 | **66.4%** | 27.6% |
| 集成 final_score (L1+L2+L3) | 116 | **0.608** | 0.593 | **73.3%** | 26.7% |
| 未映射靶点 L2 (用 name 作文本) | 48 | 0.548 | 0.545 | 62.5% | 22.9% |

- **对照上轮**: 通用(QED) AUC=0.469 → 本次无 L2 基线 0.493 (一致, 贴随机); 升级后真实 L2 **中位数 0.596、集成 0.608**, 明确越过 0.5。
- **诚实校准**: 上轮"靶点感知 0.99"是**逐靶点最大 Tanimoto-to-actives (留一法, 偷看 actives)** 的乐观上界; 本次是**真实训练模型在外部 ChEMBL 靶点上的泛化**, BindingDB 训练的 L2 跨分布到 ChEMBL 靶点达 ~0.60 中位数是合理且可信的, 而非虚高。
- **集成 final 略优于纯 L2** (0.608 vs 0.596, 73.3% vs 66.4% 过 0.5): L1/L3 (类药性/安全性) 对"真实药"额外提供少量信号。

**结论**: L2 升级达成 "AUC 回到 0.5 以上" 的目标 (中位数 + 多数靶点), 但仍非完美 —— 约 1/3 映射靶点的 AUC 仍在 0.5 以下, 说明单靠 BindingDB 训练的 L2 在部分靶点 (尤其化学型难分、或靶点命名/分布偏移大者) 仍有提升空间。

---

## 3. 柔性对接演示 (NA / 3TI6, smina --flex)

> 运行: `flexible_docking.demo_na_flexible()`, 结果见 `na_flexible_docking_demo.json`。

- 受体: `3TI6_protein_only_obabel.pdbqt` (刚性); 柔性残基: H1N1 NA 活性腔 14 个催化残基 (`A:118,119,151,152,179,227,247,248,292,293,368,370,398,401`), 经 smina `--flexres` 自动拆分为刚性+柔性。
- 配体: oseltamivir carboxylate (NA 活性), caffeine (非 NA 对照) — 各跑 刚性(对照) 与 柔性 两种对接, 记录最佳亲和力 (kcal/mol)。(演示用精简参数: 4 个核心活性腔残基作柔性 + exhaustiveness 4, 以快速验证管线)
- 意义: 证明柔性对接腿已接入并可执行 (取代单构象刚性); 真实 Boltz-2/MCCE 上游待环境就绪后补全, 即可形成完整 "多构象受体 + 质子化态 → 柔性对接" 物理接地判别器。

| 配体 | 刚性亲和力 (kcal/mol) | 柔性亲和力 (kcal/mol) |
|---|---|---|
| oseltamivir carboxylate (NA 活性) | **-6.0** | -5.0 |
| caffeine (非 NA 对照) | -5.9 | -4.7 |

**关键解读**: 刚性对接下 NA 活性药 (-6.0) 与非活性对照 caffeine (-5.9) **几乎打平** —— 这与之前 "刚性单构象 smina 在 NA 上反富集 (AUC 0.33)" 的结论一致: 单构象对接无法区分活性/非活性。柔性化后侧链获得自由度, 最优 pose 略变 (得分略升, 即结合略弱)。这进一步说明:**对接物理腿本身不足以判别靶点特异性, 必须依赖已接入的 L2 靶点感知分; 而 Boltz-2 多构象 + MCCE 质子化态 (捕捉诱导契合/质子化依赖结合) 是让物理腿真正具备判别力的下一步升级方向。**

---

## 4. 交付物清单

| 文件 | 说明 |
|---|---|
| `评分_work_package/评分/l2_bindingdb.py` | 新建 L2 靶点感知模块 (特征+模型+推理) |
| `评分_work_package/评分/models/bindingdb_l2/l2_model.joblib` | **保存的 L2 权重** (MLP+LR) |
| `评分_work_package/评分/scoring.py` | 已接入 `Layer2BindingDB` (默认 L2, 权重 0.50) |
| `scientific_validation/multitarget_benchmark/train_l2_bindingdb.py` | L2 训练脚本 |
| `scientific_validation/multitarget_benchmark/rerun_panel_real_l2.py` | 全面板重跑脚本 |
| `scientific_validation/multitarget_benchmark/real_l2_rerun_results.json` | 重跑逐靶点 AUC + 汇总 |
| `scientific_validation/multitarget_benchmark/flexible_docking.py` | Boltz-2+MCCE+smina --flex 编排模块 + NA 演示 |
| `scientific_validation/multitarget_benchmark/na_flexible_docking_demo.json` | NA 柔性对接演示结果 |

---

## 5. 下一步建议 (针对仍低于 0.5 的 ~1/3 靶点)

1. **扩大 L2 训练面**: 融合 ChEMBL 51G 对齐产物 + OpenTargets 药-靶点证据, 提升对冷门靶点 (GPCR/病毒靶) 的覆盖与泛化。
2. **补全 Boltz-2 + MCCE 上游**: 安装 `einops` 并验证 `boltz.predict`; 提供真实 MCCE 二进制, 形成 "多构象 + 质子化态 → 柔性对接" 物理接地判别器, 与 L2 并联。
3. **靶点命名对齐**: 为 48 个未映射靶点补 ChEMBL↔BindingDB 同义词表 (UniProt 桥接), 把覆盖率从 71% 提到接近 100%。
4. **集成判别器**: 将柔性对接物理分与 L2 靶点感知分在 MoleculeScorer 中融合 (加权或学习融合), 在 NA 等已有结构/数据的靶点上先验证 lift。
