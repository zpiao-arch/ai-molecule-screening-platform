# Reproducible Demo

This demo uses influenza A H1N1 neuraminidase (`influenza_a_h1n1_na`) and the public RCSB structure `3TI6` as the primary target story. Oseltamivir is treated as a known positive control, not as a newly discovered molecule.

## 1. Install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Optional real-chemistry environment:

```bash
python -m pip install -r requirements-optional.txt
```

## 2. Cache Selected Public Structure Data

```bash
python scripts/download_public_data.py --out data/public_cache --pdb 3TI6
```

This creates a small local cache and manifest. It is not required for the proxy-only CLI path, but it documents a real structure source.

## 3. Run CLI Closed-Loop Demo

```bash
PROJECT=demo_runs/flu_na_demo

python ai_mol_loop/ai_mol_loop.py init "$PROJECT"
python ai_mol_loop/ai_mol_loop.py target-select "$PROJECT" --disease 甲流 --top 5
python ai_mol_loop/ai_mol_loop.py brief-from-target "$PROJECT" --disease 甲流 --force
python ai_mol_loop/ai_mol_loop.py evidence-stage2 "$PROJECT" --disease influenza --top 5 --offline
python ai_mol_loop/ai_mol_loop.py run-demo "$PROJECT" --rounds 2 --n 24 --top 6
```

## 4. Stage 4 Real-Library Assets

If RDKit is installed:

```bash
python ai_mol_loop/ai_mol_loop.py stage4-real "$PROJECT" \
  --round 2 \
  --target influenza_a_h1n1_na \
  --top 8 \
  --decoys 8 \
  --rescore
```

This produces descriptors, similarity-to-control tables, decoys, SDF assets, receptor package metadata, and a docking plan. It does not claim real docking unless `--run-docking` is used with a valid receptor, prepared ligands, and a real Vina/GNINA executable.

## 5. Dashboard and Delivery

```bash
python ai_mol_loop/ai_mol_loop.py stage5-dashboard "$PROJECT" --round 2
python ai_mol_loop/ai_mol_loop.py stage6-validate "$PROJECT" --round 2
python ai_mol_loop/ai_mol_loop.py stage7-package "$PROJECT" --round 2
```

## 6. Web Demo

```bash
./start_web.sh
```

Open [http://localhost:8765/](http://localhost:8765/). Use the project registration page to create or select a project, then run stages from the frontend.

## Expected Interpretation

Expected outputs are computational-screening artifacts: target evidence matrices, candidate tables, descriptor tables, docking plans, dashboard summaries, and delivery reports. They are not wet-lab or clinical validation.
