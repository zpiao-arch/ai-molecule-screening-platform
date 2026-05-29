# AI Molecule Screening Platform

候选分子筛查与验证工作台：一个面向早期药物发现演示和工程验证的 CLI + FastAPI Web 平台。项目把靶点证据整理、候选分子生成/导入、化学规则过滤、RDKit 描述符校验、docking 计划/可选真实 docking、对照校准、Dashboard 和交付报告组织成一个可复现的计算筛选流程。

> This repository supports computational screening and validation planning only. It does not prove biological activity, potency, toxicity, safety, dosing, clinical usefulness, or therapeutic efficacy.

## What Is Included

- `ai_mol_loop/`: core CLI workflow for Stage 1-8 assets.
- `webapp/`: FastAPI backend and single-page frontend.
- `webapp/static/vendor/3Dmol-min.js`: local 3Dmol.js viewer bundle for receptor/pose display.
- `ai_mol_loop/targets/influenza/`: small public seed evidence for influenza targets, known drugs, PDB IDs, and offline evidence summaries.
- `scripts/download_public_data.py`: selected public-data downloader for small local caches.
- `scripts/download_ai_molecule_repos_api.py`: optional source mirror helper for REINVENT4 / DockStream / GraphINVENT / DrugEx style projects.
- `tests/`, `webapp/tests/`, `ai_mol_loop/tests/`: regression tests and API/frontend contract tests.

## What Is Not Included

- Runtime projects and generated results: `webapp/projects/`, `outputs/`, `deliverables/`.
- Large third-party repositories and binaries: `ai_drug_eval_tools/`, `ai_molecule_design_repos/`, Vina/GNINA builds.
- Large structure/data caches: PDB/PDBQT/SDF docking outputs, full ChEMBL/BindingDB/DrugBank dumps.
- API keys, `.env` files, commercial software configs, or restricted databases.

## Quick Start

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt

python -m py_compile webapp/server.py ai_mol_loop/ai_mol_loop.py
./start_web.sh
```

Open [http://localhost:8765/](http://localhost:8765/).

For optional real-chemistry and docking features, install `requirements-optional.txt` and external tools such as AutoDock Vina, OpenBabel, Meeko, and PoseBusters. Without them, the app still runs in proxy/planning mode and will not fabricate docking scores.

## Minimal CLI Demo

```bash
PROJECT=demo_runs/flu_na_demo

python ai_mol_loop/ai_mol_loop.py init "$PROJECT"
python ai_mol_loop/ai_mol_loop.py target-select "$PROJECT" --disease 甲流 --top 5
python ai_mol_loop/ai_mol_loop.py brief-from-target "$PROJECT" --disease 甲流 --force
python ai_mol_loop/ai_mol_loop.py evidence-stage2 "$PROJECT" --disease influenza --top 5 --offline
python ai_mol_loop/ai_mol_loop.py run-demo "$PROJECT" --rounds 2 --n 24 --top 6
python ai_mol_loop/ai_mol_loop.py stage5-dashboard "$PROJECT" --round 2
python ai_mol_loop/ai_mol_loop.py stage6-validate "$PROJECT" --round 2
python ai_mol_loop/ai_mol_loop.py stage7-package "$PROJECT" --round 2
```

See [REPRODUCIBLE_DEMO.md](REPRODUCIBLE_DEMO.md) for the longer influenza NA path and Stage 4 notes.

## Public Data

The repository stores only small seed metadata. Selected public structures can be cached locally:

```bash
python scripts/download_public_data.py --out data/public_cache --pdb 3TI6 --pdb 6FS6
```

See [DATA_SOURCES.md](DATA_SOURCES.md) for source stability and local-storage policy.

## GitHub Release Policy

This repo is intended to stay small and reviewable. Generated assets should be produced by commands, not committed. If a demo result needs to be shared, use a GitHub Release asset or external storage instead of committing large files to Git history.

## License

Project code is released under Apache-2.0. Bundled third-party assets are documented in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
