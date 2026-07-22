# Four-Level Scoring Core

本目录是正式小批量 CLI 的评分核心。当前发布并验证的路径为：

| 层 | 权重 | 正式后端 | 输出 |
|---|---:|---|---|
| L1 | 20% | RDKit 分子描述符、QED、规则化 SA/Lipinski | `layer1_score` |
| L2 | 50% | BindingDB 靶点文本哈希 + sklearn MLP | `docking_normalized` |
| L3 | 20% | 4 个 ADMET sklearn 模型 + RDKit/SMARTS 规则 | `admet_score`、`bbb_prob` |
| L4 | 10% | UniMol 表征 + FDA 参考嵌入 | `unimol_score` |

综合分公式固定为：

```text
final_score = 0.20*L1 + 0.50*L2 + 0.20*L3 + 0.10*L4
```

## 运行

从仓库根目录运行。先按根目录 `README.md` 合并离线资产并完成 runtime doctor：

```bash
PYTHONPATH=scoring python scoring/scoring.py \
  --input examples/hiv_candidates.csv \
  --target HIV-1_protease \
  --strict-backends \
  --output outputs/hiv_scores.csv
```

正式 CLI 的 `--l2` 只暴露随发布包验证的 `bindingdb`。源码中保留的 `Layer2BindingDBSeq` 和 `Layer2DeepPurpose` 是历史/实验实现，需要额外模型、序列数据或第三方环境，不属于开箱即用承诺。

## 输入契约

CSV 必须同时包含非空的 `id` 和 `smiles` 列；`id` 必须唯一。缺表头、空值或重复 ID 会在加载阶段报错，不会进入评分或对接。

```csv
id,smiles
mol_001,CC(=O)Oc1ccccc1C(=O)O
mol_002,c1ccccc1N
```

## 严格模式

`--strict-backends` 会在 L1/L2/L3/L4 任一后端加载失败或单分子状态非 `ok` 时退出。非严格模式可能保留失败行，但每层均输出：

- `layerN_status`
- `layerN_backend`
- `layerN_model_asset_id`

因此不能仅凭 `final_score` 判断某行是否完成真实四层计算。

正式 `bindingdb` 权重缺失时，非严格模式不会切换到未打包的 DeepPurpose；L2 会明确输出 `BindingDB-L2-unavailable` 和 `failed:backend_unavailable:*` 状态，L2 贡献为 0。需要完整四层结果时应使用 `--strict-backends` 并先通过 runtime doctor。

## 结构对接

注册表当前只有 CHEMBL2051 的有效受体。`--mode auto` 命中有效注册资产时进入级联；`--mode cascade` 强制级联，受体、二进制或零成功对接都会使命令失败。

手工受体必须同时提供对接盒，显式参数优先于注册表：

```bash
PYTHONPATH=scoring python scoring/scoring.py \
  --input examples/flu_candidates.csv \
  --target CHEMBL2051 \
  --mode cascade \
  --receptor /path/to/receptor.pdbqt \
  --box-center=-28.914,14.334,20.794 \
  --box-size=23.585,20.45,24.18 \
  --strict-backends \
  --output outputs/flu_cascade.csv
```

历史材料中的 `0.935` 来自人工平衡/特定展示集，不是本 CLI 在真实 1000×10000 池上的通用能力。当前冻结运行的统计边界见根目录 `VALIDATION.md`。

## 主要输出字段

| 字段 | 含义 |
|---|---|
| `layer1_score` | L1 分子质量分，0-1 |
| `docking_normalized` | L2 BindingDB 结合概率，0-1 |
| `admet_score` | L3 综合安全性分，0-1 |
| `bbb_prob` | BBBP 模型的穿透概率，不是 logBB |
| `unimol_score` | L4 FDA 参考相似度分，0-1 |
| `final_score` | 四层固定权重综合分，0-1 |
| `gate_status` | 规则门槛 `PASS`/`FAIL` |
| `structure_docking_status` | 结构对接状态（级联时） |
| `docking_affinity_kcal_mol` | smina affinity（级联成功时） |
| `ligand_efficiency` | affinity / heavy atoms |

没有显式确认标签时，CLI 末尾的虚筛 benchmark 固定显示 `not_evaluated`，不会用模型自己的 top 分子充当阳性。

## 资产与限制

模型和受体由配套离线资产包提供；smina/obabel 需自行安装并通过 `DOCK_BIN_DIR` 或各自的 `SMINA_BIN`/`OBABEL_BIN` 配置。全部分数是计算预测，不替代冷启动、时间外、骨架外或湿实验验证。
