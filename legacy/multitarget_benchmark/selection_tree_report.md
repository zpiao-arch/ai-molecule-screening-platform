# 药物分子打分：多模型 vs 选择树 — 调查结论与 CLI 集成

> 调查问题（用户）：这种"坍缩 bug"是否需要**同时训练不同模型**来修复？
> 多版本模型都有，是否可能需要一个**选择树**，按规则合适运行并**整合进 CLI**？

## 1. 调查：这个 bug 需要多模型来修吗？—— 不需要"在版本间挑"

我们手上有 3 个 L2 版本，在 **471 个标注靶 / 10万药库** 同口径下对比：

| L2 版本 | 中位 AUC | 反向率 | 判读 |
|---|---|---|---|
| **哈希基线** (bindingdb, 7/14) | **0.6889** | 19.1% | ✓ 通用最佳 |
| CE 重训 (检索负采样 K=2) | 0.5705 | 38.2% | ✗ 聚合退化 |
| BPR 重训 (配对排序损失) | 0.5065 | ~49% | ✗ 分数尺度失校准, ≈随机 |

**关键诊断：**
- **CE 重训**把"好靶"也坑了：`CHEMBL4069` 0.999→0.001。根因 = 二元 CE + 极端检索负例 → 模型学成"几乎都是非活性→输出≈0"。
- **BPR 重训**训练 loss 收敛到 0.022，但评测≈随机。诊断显示其**分数尺度失校准**：NA 活性 −3.0 vs 库 −7.0（好），而 `CHEMBL4069` 活性 −31.5 vs 库 −8.2（活性被压到库底）。BPR 只约束"配对内 活性>阴性"，对绝对尺度无约束 → 多数靶活性被推到全库最低。
- **结论**：单模型"重训修复"在聚合域**两次都失败**。原始哈希模型（在域测定负例上训练）反而是最好的通用模型。

### "按版本选择"的天花板（决定性证据）
假设有一个**完美路由器**逐靶选 哈希/CE 中最优的那个：
- Oracle 中位 AUC = **0.7343**，仅比哈希基线高 **+0.045**。
- 用廉价特征 `n_pos`（推理期可得）做路由规则：最好也只能恢复到 oracle 的 **10%**，且 `n_pos≥20` 的阈值反而让中位降到 0.663（比哈希更差）。
- → **在 L2 版本间做选择树是无价值的**：天花板极低（≈0.734）且无法用可获取的信号恢复。

## 2. 真正有价值的选择树：按"可用证据"分支（非按模型版本）

真正拉开差距的不是"哪个 L2 版本"，而是**该靶点"有没有可用的 3D 受体"**：

| 分支 | 路径 | 实测 AUC |
|---|---|---|
| 无受体 → 文库分支 | L2 + ADMET + UniMol | 聚合中位 **0.689** |
| 有受体 → 级联分支 | L2 粗筛 + smina 对接精排(LE) + 融合 | NA 平衡集融合 **0.9242** |

这是基于**可用证据类型**的规则分支，干净、可解释、已被数据证明。判别树见 `selection_tree_diagram.svg`。

## 3. CLI 集成（已完成并验证）

新增模块 `评分_work_package/评分/pipeline_router.py` + 受体注册表 `receptor_registry.json`，并接入 `scoring.py`：

- **`--mode auto|cascade|library`**（默认 `auto`）：`auto` 下经 `route()` 查受体注册表自动路由。
- 命中受体 → 自动注入 `--receptor` / 对接盒 / `--dock-rerank` / 默认融合权重 **0.5**（复现 NA 0.9242）。
- **`L2_MODEL_PATH` opt-in**（沿用你之前问的）：任一分支可改用以检索重训模型（NA 上 0.485→0.935），默认不启用（因聚合退化）。
- **扩展方式**：在 `receptor_registry.json` 的 `entries` 添加靶点 → `pdbqt` 路径 + 对接盒中心/尺寸，该靶即自动走级联。

### 用法
```bash
# 自动: NA 自动走级联对接精排, EGFR 自动走文库分支
python scoring.py -i mols.csv -o out.csv -t CHEMBL2051 --mode auto
python scoring.py -i mols.csv -o out.csv -t EGFR       --mode auto

# 强制分支
python scoring.py -i mols.csv -o out.csv -t X --mode cascade   # 需注册表有受体
python scoring.py -i mols.csv -o out.csv -t X --mode library

# opt-in: 让 NA 用检索重训 L2 (预筛更强)
L2_MODEL_PATH=评分_work_package/评分/models/bindingdb_l2_retrieval/l2_model.joblib \
  python scoring.py -i mols.csv -o out.csv -t CHEMBL2051 --mode auto
```

### 验证结果
| 测试 | 路由 | 结果 |
|---|---|---|
| EGFR / `--mode auto` | **library** | L2+ADMET 输出完整, 无对接 |
| CHEMBL2051 / `--mode auto` | **cascade** | 自动命中受体 → 4/4 成功对接 → 融合重排; zanamivir 被 L2 正确识别(0.68) |

## 4. 边界与下一步
- **受体限制**：当前注册表仅 NA(3TI6) 有受体；其它靶仍走文库分支（0.689）。要扩更多靶的 0.8+ 闭环，需补充对应 pdbqt + 对接盒。
- **BPR 尺度失校准**：值得后续修（全库 softmax / 难负采样 / 加 sigmoid 校准），但属锦上添花——即便修好，单模型仍难稳超哈希（见 §1 oracle）。
- **环境注记**：UniMol 在本环境缺 `unimol_tools` 模块（L4 跳过，不影响 L1/L2/L3）；ADMET 模型因 sklearn 1.7.2↔1.6.1 版本差有 `InconsistentVersionWarning`（无害）。
