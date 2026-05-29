# Open Source Manifest

## Repository Scope

This repository contains source code, tests, small public seed data, and reproducibility instructions for the AI Molecule Screening Platform.

## Core Code

- `ai_mol_loop/ai_mol_loop.py`: CLI and Stage 1-8 workflow logic.
- `webapp/server.py`: FastAPI backend, project APIs, Stage APIs, artifact endpoints, and health checks.
- `webapp/static/index.html`: web UI.
- `scripts/download_public_data.py`: selected public data cache builder.
- `scripts/build_product_delivery.py`: local delivery-package helper.
- `scripts/frontend_closed_loop_e2e.py`: optional Playwright E2E driver.

## Seed Data

- `ai_mol_loop/targets/influenza/target_catalog.json`
- `ai_mol_loop/targets/influenza/known_drugs.csv`
- `ai_mol_loop/targets/influenza/pdb_structures.csv`
- `ai_mol_loop/targets/influenza/evidence/`

These are lightweight metadata snapshots, not full public databases.

## Excluded from Git

- generated projects and exports;
- large third-party repositories;
- external docking binaries;
- downloaded PDB/PDBQT/SDF caches;
- secrets and `.env` files;
- commercial/restricted datasets.

## Recommended GitHub Repository

- Owner: `zpiao-arch`
- Repository: `ai-molecule-screening-platform`
- Visibility: public after final secret scan.

## Release Package

Use the generated zip under `open_source_release/` for upload if GitHub CLI or SSH authentication is unavailable.
