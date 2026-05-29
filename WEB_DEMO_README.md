# Web Demo 使用说明

## 一键启动

在项目根目录运行：

```bash
./start_web.sh
```

脚本会启动 `webapp/server.py`，默认地址是 `http://localhost:8765/`。如果 8765 端口已有服务，脚本会先检查 `http://localhost:8765/api/health`，避免重复启动。

## 健康检查

浏览器或终端访问：

```text
http://localhost:8765/api/health
http://localhost:8765/api/system/health
http://localhost:8765/api/product/bootstrap-health
```

`/api/health` 返回 `status: ok` 表示后端、前端入口和项目目录可用。`/api/system/health` 返回 `overall_status`、运行环境、前端资源、计算工具和交付链路检查结果。
`/api/product/bootstrap-health` 会把启动脚本、Demo 预检、Stage 8 下一步和交付能力合并成产品级自检包。

也可以在项目根目录运行：

```bash
python3 health_check.py
```

该命令会读取 `/api/system/health`，输出各检查项、`recommended_commands` 和合规边界。

## 产品指挥台 Demo

打开网页后进入「产品指挥台」。Stage 8 页面提供：

- Demo 模式：优先加载 `flu_na_real_demo`。
- 一键生成 Demo 包：自动串联 Stage 5 Dashboard、Stage 6 验证运营、Stage 7 交付包。
- 完整交付导出：点击「完整交付导出」会生成 Stage 8 全量 ZIP、manifest、系统健康、Target Pack 校验和科学就绪度材料。
- 后台任务：点击「后台导出任务」会提交 `/api/jobs`，前端轮询 `/api/jobs/{job_id}` 并回填导出结果。
- 下载材料：交付清单、执行摘要、质量门、候选分子排序表等可直接打开或下载。
- 阶段入口：点击阶段卡可跳转到对应 Stage 页面继续查看细节。

## Product Ops API

常用接口：

```text
GET  /api/system/health
GET  /api/product/bootstrap-health
GET  /api/generator-adapters
GET  /api/jobs
POST /api/jobs
GET  /api/jobs/{job_id}
GET  /api/projects/{name}/target-pack/validate?round=1
GET  /api/projects/{name}/stage4/operator-guide?round=1
GET  /api/projects/{name}/stage8/action-plan?round=1
GET  /api/projects/{name}/scientific-readiness?round=1
POST /api/projects/{name}/stage8/full-export
```

前端应优先用「环境检查」页判断可演示状态，用「产品指挥台」页触发完整交付导出。

## 标准交付构建

需要打包单个项目时运行：

```bash
python3 scripts/build_product_delivery.py --project flu_na_real_demo --round 1
```

脚本会先运行测试，然后生成 Stage 8 full export、bootstrap health、Stage 8 action plan、generator adapters 清单和标准交付 zip。演示赶时间时可加 `--skip-tests`，正式交付建议保留测试。

## 靶点与生成器扩展

`/api/target-catalog` 已从甲流目录扩展为多疾病演示目录，保留 influenza NA/PA/M2，同时加入 EGFR、BACE1、HIV-1 protease 等公开结构和控药更充分的靶点。

`/api/generator-adapters` 汇总当前可接入的候选生成来源：内置 proxy SMILES、OpenAI prompt 入口、REINVENT4、DrugEx、DiffSBDD、ColabFold。外部仓库默认作为 Stage 3 CSV/结构包适配器，不直接声称模型训练或真实药效。

## 边界

本系统只做计算筛选、真实库校验、验证规划和交付材料打包。不声称药效、活性、毒性、安全性、剂量、临床收益或真实疗效。
