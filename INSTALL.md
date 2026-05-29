# Installation

## Core Environment

Python 3.9+ is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run a basic import/compile check:

```bash
python -m py_compile webapp/server.py ai_mol_loop/ai_mol_loop.py scripts/download_public_data.py
python -m unittest tests/test_public_data_downloader.py -v
```

## Start Web App

```bash
./start_web.sh
```

Default URL: [http://localhost:8765/](http://localhost:8765/).

Health checks:

```bash
curl http://localhost:8765/api/health
python health_check.py
```

## Optional Chemistry/Docking Environment

The core app can run without these tools. Install them only when you want real descriptors, ligand preparation, docking, pose QC, or benchmark workflows.

```bash
python -m pip install -r requirements-optional.txt
```

Common external executables:

- `vina`: AutoDock Vina docking backend.
- `gnina`: optional CNN docking/rescoring backend.
- `obabel`: OpenBabel ligand/receptor conversion.
- `mk_prepare_ligand.py`: Meeko ligand preparation.
- `bust`: PoseBusters pose quality check.

The app reports missing tools in the environment page and Stage 4 operator guide. Missing tools should produce `planning_only` or `skipped` states, not fake docking scores.

## OpenAI API

OpenAI candidate generation is optional. Use environment variables or a local key file. Do not commit keys.

```bash
export OPENAI_API_KEY="<your-openai-api-key>"
python ai_mol_loop/ai_mol_loop.py stage3-screen demo_runs/flu_na_demo --round 2 --use-openai --n 24 --top 6
```

The code is designed to use API keys in memory for the request and avoid writing them into outputs.
