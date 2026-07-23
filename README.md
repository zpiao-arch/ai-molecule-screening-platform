# Four-Level Molecule CLI

可审计的四级药物分子筛选 CLI：L1 分子质量、L2 靶点结合、L3 ADMET、L4 UniMol 药物相似度，并提供质量门槛、受体级联和可复算的 `1000 × 10000` 批量验证工具。

> 发布状态：当前 `LICENSE` 只允许 source-available 使用，不授予开源再分发权。合规边界见 [COMPLIANCE.md](COMPLIANCE.md)。

## 发布边界

本仓库是完整 source-available 源码包，包含主 CLI、批量验证模块、测试、案例、训练/预计算脚本、历史分析脚本归档以及冻结 run 的紧凑审计快照。模型权重、完整 ChEMBL/BindingDB 数据、受体文件和 10M 候选级 score 分区不放入 GitHub 源码包。对应的本机离线资产包为：

`../four-level-molecule-cli-offline-assets`

运行 CLI 时可通过 `--asset-root /path/to/offline-assets` 显式指向离线资产根；不要把权重、受体或二进制复制进源码包。源码包自带的 `scoring/data/` 只是小型案例输入，不能替代完整 benchmark 数据。smina/obabel 是平台相关二进制，不随资产包分发，需自行安装并设置 `DOCK_BIN_DIR`，也可分别设置 `SMINA_BIN` 和 `OBABEL_BIN`。legacy 行级数据不随 GitHub 源码归档分发，见 [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md)。

源码包同时包含 prompt-first 的 [Open Molecule Lab](apps/open-molecule-lab/)。它可以在本机上传 `id,smiles` CSV、执行严格预检、调用现有四级 CLI，并输出带状态和校验 manifest 的运行证据包。缺少外部资产时运行会 `blocked`，不会生成假分数。

完整文件盘点见 [FULL_SOURCE_INVENTORY.md](FULL_SOURCE_INVENTORY.md)。`legacy/multitarget_benchmark/` 是历史实验和诊断脚本归档，不是当前 1000×10000 CLI 的必要运行路径。

## 已验证范围

- `1000` 个靶点 × `10000` 候选，共 `10000000` 条四级 score。
- L1/L2/L3/L4 failure 均为 `0`，fit-pair contamination 为 `0`。
- 正式 CLI 与批量路径的 20 分子等价性复核最大差为 `0`，容差 `1e-4`。
- CHEMBL2051 的 top-300 对接包含逐分子状态、affinity、ligand efficiency 和命令参数。
- 冻结运行基线测试：`43 passed`；当前 source-only 默认测试：`97 passed, 10 deselected`；Open Molecule Lab source-only contract、worker 进程组生命周期、前端 build，以及两分子的真实 library + smina cascade 离线运行均已单独验证。测试标记和运行命令见下文。

结果的统计边界见 [VALIDATION.md](VALIDATION.md)，compact snapshot 的证据在 [validation/frozen_run](validation/frozen_run) 中。运行 `python -m scientific_validation.four_level_cli_1kx10k.verify_snapshot --snapshot-dir validation/frozen_run` 验证源码包快照；完整 1 GB 本机 run 仍使用 `verify_run --strict`。pair-heldout 结果不是冷靶点、冷分子、时间外推或湿实验验证；当前受体注册表只有 CHEMBL2051。

Open Molecule Lab 的两分子真实本机检查：

```bash
OPEN_MOLECULE_PYTHON=/path/to/python3.11 \
OPEN_MOLECULE_ASSET_ROOT=/path/to/four-level-molecule-cli-offline-assets \
SMINA_BIN=/path/to/smina \
OBABEL_BIN=/path/to/obabel \
npm --prefix apps/open-molecule-lab run real-run-check
```

## 环境

推荐 Python 3.11：

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install .
python -m pip install pytest
```

离线资产就位后先检查资产：

```bash
python scripts/verify_assets.py \
  --asset-root ../four-level-molecule-cli-offline-assets \
  --manifest ../four-level-molecule-cli-offline-assets/ASSET_MANIFEST.json
```

安装后也可直接使用 `four-level-molecule`, `four-level-benchmark`, `four-level-doctor` 和 `four-level-verify-snapshot` 四个命令；未安装时保留 `PYTHONPATH=scoring python scoring/scoring.py` 的脚本用法。

严格运行时预检：

```bash
python -m scientific_validation.four_level_cli_1kx10k.runtime_doctor \
  --scoring-dir scoring --strict
```

## 单分子/小批量评分

```bash
PYTHONPATH=scoring python scoring/scoring.py \
  --input examples/hiv_candidates.csv \
  --target HIV-1_protease \
  --asset-root ../four-level-molecule-cli-offline-assets \
  --strict-backends \
  --output outputs/hiv_scores.csv
```

## 1000×10000 批量闭环

数据资产就位后，在仓库根目录执行：

```bash
export DOCK_BIN_DIR=/path/to/smina-local/bin
python -m scientific_validation.four_level_cli_1kx10k.batch_cli \
  --run-dir runs/$(date -u +%Y%m%dT%H%M%SZ) \
  prepare --targets 1000 --pool-size 10000 --seed 42
```

随后对同一 `--run-dir` 依次执行：

```bash
python -m scientific_validation.four_level_cli_1kx10k.batch_cli \
  --run-dir runs/<RUN_ID> cache-layers --strict --resume

python -m scientific_validation.four_level_cli_1kx10k.batch_cli \
  --run-dir runs/<RUN_ID> score --resume

python -m scientific_validation.four_level_cli_1kx10k.batch_cli \
  --run-dir runs/<RUN_ID> dock --target CHEMBL2051 --top-n 300 --resume
```

最后运行：

```bash
python -m scientific_validation.four_level_cli_1kx10k.batch_cli \
  --run-dir runs/<RUN_ID> report

python -m scientific_validation.four_level_cli_1kx10k.verify_run \
  --run-dir runs/<RUN_ID> --strict
```

`report` 会重新计算审计报告、补齐 `DESIGN.md`/失败清单/报告 checkpoint，并在完整运行目录生成完整 `MANIFEST.sha256`。源码包内的 compact snapshot 使用单独的 `snapshot.json` 和 `verify_snapshot`；阶段输入和代码哈希会写入 checkpoint；哈希不匹配时续跑会硬失败。

单批 CLI 没有确认阳性/阴性标签时，末尾的虚拟筛选 benchmark 会显示 `not_evaluated`；不会把模型自己排序的 top-5 当作阳性。

## 测试

源码包默认测试不读取外部模型、数据或平台二进制：

```bash
python -m pytest -q
```

合并离线资产后，只运行资产集成测试：

```bash
python -m pytest -q -o addopts='' -m "integration_assets and not integration_docking"
```

配置 `DOCK_BIN_DIR` 后运行完整测试（含 smina/obabel）：

```bash
DOCK_BIN_DIR=/path/to/smina-local/bin python -m pytest -q -o addopts=''
```

构建 UTF-8 安全的 GitHub 源码归档：

```bash
python scripts/build_source_release.py \
  --source-dir . \
  --output ../packages/four-level-molecule-cli-source-YYYYMMDD.zip \
  --asset-dir ../four-level-molecule-cli-offline-assets \
  --asset-output ../packages/four-level-molecule-cli-offline-assets-YYYYMMDD.tar.gz
```

## 许可证

源码包当前没有授予开源再分发许可，见 [LICENSE](LICENSE) 和 [COMPLIANCE.md](COMPLIANCE.md)。UniMol、RDKit、PyTorch、ChEMBL、BindingDB、ADMET 权重和 smina/obabel 各自受其上游许可约束；只有在权利人明确选择许可证后，才能将仓库标为 open source。
