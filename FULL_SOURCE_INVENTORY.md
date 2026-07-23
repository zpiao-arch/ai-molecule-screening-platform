# Full Source Inventory

本发布包按“可运行主路径 / 历史分析归档 / 验证证据 / 外部资产”分层，避免把本机大文件误当作源码。

## Included

- `scoring/`: 当前正式四层 CLI、L2 BindingDB、L3 ADMET、L4 UniMol、路由、对接精排、靶点解析、序列嵌入和脚本。
- `apps/open-molecule-lab/`: prompt-first 本地工作台、MoleculeSet 校验、严格预检、后台 CLI worker、状态/结果 API 和前端结果浏览器；不包含运行产物或 `node_modules`。
- `scoring/data/` 和 `examples/`: 小型案例输入。输出目录由 CLI 按需创建；旧字段/失败状态的历史 CSV 不作为正式样例分发。
- `scientific_validation/four_level_cli_1kx10k/`: 1000×10000 批量 CLI、数据契约、缓存、指标、报告、manifest 和 verifier。
- `tests/`: 四级 CLI 契约、数据集、指标和运行验证测试。
- `pytest.ini`：默认 source-only 测试过滤和 `integration_assets` / `integration_docking` 标记定义。
- `scripts/build_source_release.py`：确定性、UTF-8 安全的 GitHub 源码归档构建器；`scripts/verify_assets.py` 校验外部资产哈希。
- `scoring/asset_paths.py`：显式解析外部模型、UniMol 缓存和受体根，避免把离线资产复制进源码包。
- `COMPLIANCE.md`：当前 source-available 许可、第三方资产和分发边界。
- `legacy/multitarget_benchmark/`: 与本项目历史结果相关的训练、诊断、级联、批量和报告源码；机器可读的 `.csv`/`.json`/`.parquet` 行级研究产物不进入 GitHub 源码归档，边界见 `THIRD_PARTY_DATA.md`。
- `validation/frozen_run/`: compact `snapshot.json`、summary、报告、环境/失败/等价性证据、target metrics、阶段 checkpoint、对接前后指标和只覆盖已包含文件的 manifest；不包含 10M score 分区、1000 个 pool 文件或原始训练表。
- `docs/`: 已确认的四级 CLI 设计规格。

## Excluded

- `scoring/models/` 和其他训练权重：由离线资产包提供。
- ChEMBL/BindingDB 原始或完整对齐数据、100k library 完整特征、UniMol 全模型目录中的非必需权重。
- `scoring/receptors/` 和 smina/obabel 二进制。
- `validation/frozen_run/scores/`、`pool_manifest/`、layer cache shards、候选级大产物，以及 `legacy/**/*.csv|json|parquet`。
- Python bytecode、pytest 缓存、临时日志和本机绝对路径。
- Open Molecule Lab 的 `data/`、`runs/`、worker 日志和结果分区；公开 bundle 只保留相对路径和脱敏日志。

历史归档中的 `data_lake/...`、`评分_work_package/...` 等字符串是原实验脚本的输入/输出约定，不是当前发布包要求的目录；需要重跑历史脚本时应通过环境变量或本地路径重新配置。

源码包内所有历史归档中的本机路径已替换为 `<validated-workspace>`；这些历史脚本若要重新执行，应先按当前环境改写输入路径。
