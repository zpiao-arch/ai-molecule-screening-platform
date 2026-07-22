# Four-Level CLI Release Hardening V2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the eight confirmed release defects and rebuild a trustworthy source-available GitHub package without changing scientific scores.

**Architecture:** Centralize asset trust and CLI scalar validation in small reusable helpers, keep the complete-run verifier unchanged, and add a distinct compact-snapshot verifier. Preserve script execution while making `scoring` installable as a package. Public packaging excludes unreviewed machine-readable legacy data but retains historical source and reports.

**Tech Stack:** Python 3.11, argparse, hashlib/JSON manifests, pytest, setuptools/PEP 621, GitHub Actions, RDKit, pandas/pyarrow, joblib, PyTorch/Uni-Mol, smina/Open Babel.

---

### Task 1: Enforce trusted asset loading and align runtime doctor

**Files:**
- Create: `scoring/asset_integrity.py`
- Modify: `scoring/l2_bindingdb.py`
- Modify: `scoring/scoring.py`
- Modify: `scoring/scripts/unimol_scorer.py`
- Modify: `scoring/scripts/unimol_embedding.py`
- Modify: `scientific_validation/four_level_cli_1kx10k/runtime_doctor.py`
- Test: `tests/test_four_level_cli_contract.py`

- [ ] **Step 1: Write failing tests** for a valid asset, mismatched digest, uncovered asset, explicit unsafe opt-in, doctor manifest failure, and `SMINA_BIN`/`OBABEL_BIN` resolution.

```python
def test_asset_verification_rejects_mismatch_before_joblib_load(tmp_path, monkeypatch):
    model = tmp_path / "scoring/models/bindingdb_l2/model.joblib"
    model.parent.mkdir(parents=True)
    model.write_bytes(b"not trusted")
    manifest = tmp_path / "ASSET_MANIFEST.json"
    manifest.write_text(json.dumps({"files": {"scoring/models/bindingdb_l2/model.joblib": "0" * 64}}))
    monkeypatch.setenv("FOUR_LEVEL_ASSET_MANIFEST", str(manifest))
    with pytest.raises(AssetIntegrityError, match="sha256 mismatch"):
        verify_asset(model)
```

- [ ] **Step 2: Run targeted tests and confirm expected failures.**

Run: `python -m pytest -q tests/test_four_level_cli_contract.py -k 'asset_verification or doctor_honors_individual'`

- [ ] **Step 3: Implement `AssetIntegrityError`, manifest discovery, cached SHA verification, and unsafe opt-in reporting.** The public API is:

```python
class AssetIntegrityError(RuntimeError):
    pass

def verify_asset(path: str | Path, *, manifest_path: str | Path | None = None) -> dict[str, object]:
    """Return path/relative/sha/manifest evidence or raise before deserialization."""

def verify_manifest(root: str | Path, manifest_path: str | Path) -> dict[str, object]:
    """Verify every manifest file and report missing, mismatches, and checked count."""
```

- [ ] **Step 4: Call `verify_asset` immediately before every joblib/pickle/NumPy/PyTorch load used by production scoring.** Keep training-script output writes out of this gate.

- [ ] **Step 5: Make runtime doctor use the shared manifest verifier and production-equivalent binary lookup.** `overall` must require clean manifest verification unless unsafe opt-in is explicitly active, in which case doctor reports `unsafe_override` and strict mode fails.

- [ ] **Step 6: Run targeted and full source tests.**

### Task 2: Harden CLI argument and CSV contracts

**Files:**
- Modify: `scoring/scoring.py`
- Modify: `scientific_validation/four_level_cli_1kx10k/batch_cli.py`
- Test: `tests/test_four_level_cli_contract.py`

- [ ] **Step 1: Write failing tests** proving BOM CSV acceptance; parser rejection of zero/negative counts; forced-library rejection of explicit `--dock-mode` and `--cascade-top-n`; auto-library rejection of `--dock-fusion`; and explicit receptor plus valid box selecting cascade in auto mode.

```python
def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return parsed

def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return parsed
```

- [ ] **Step 2: Run the new tests and confirm each fails for the old behavior.**

- [ ] **Step 3: Use `utf-8-sig`; wire positive/non-negative validators into both parsers; change docking-only defaults to `None`; validate effective route conflicts after routing; normalize omitted values only after validation.**

- [ ] **Step 4: Run targeted parser/CLI tests and the complete default suite.**

### Task 3: Record failed stage checkpoints

**Files:**
- Modify: `scientific_validation/four_level_cli_1kx10k/provenance.py`
- Modify: `scientific_validation/four_level_cli_1kx10k/batch_cli.py`
- Modify: `scientific_validation/four_level_cli_1kx10k/verify_run.py`
- Test: `tests/test_four_level_run_verifier.py`

- [ ] **Step 1: Write failing tests** for writing/validating a failed checkpoint and for `batch_cli.main()` converting an existing `in_progress` checkpoint to `failed` when a stage raises.

```python
def fail_stage_checkpoint(run_dir: str | Path, stage: str, error: BaseException) -> Path | None:
    """Preserve existing inputs/outputs and atomically record a bounded failure."""
```

- [ ] **Step 2: Confirm failures because `failed` is currently rejected and main leaves `in_progress`.**

- [ ] **Step 3: Accept `failed`, add error metadata and JSONL event, wrap stage dispatch, and keep `require_complete=True` behavior unchanged.** If no checkpoint exists, do not fabricate inputs.

- [ ] **Step 4: Verify retry with matching input can replace `failed` with `in_progress`, then run provenance/verifier tests.**

### Task 4: Make receptor protonation atomic and stale-proof

**Files:**
- Modify: `scoring/dock_rerank.py`
- Modify: `scoring/scoring.py`
- Test: `tests/test_four_level_cli_contract.py`

- [ ] **Step 1: Write a failing test** where Open Babel returns 1 while an adjacent stale `.protonated.pdb` exists; assert `DockingReranker` is never constructed and the stale file is untouched.

```python
@contextmanager
def protonated_receptor(receptor: str | Path, obabel_bin: str, *, timeout: int = 300):
    with tempfile.TemporaryDirectory(prefix="four-level-receptor-") as directory:
        output = Path(directory) / "receptor.pdb"
        result = subprocess.run(
            [obabel_bin, str(Path(receptor).resolve()), "-O", str(output), "-p", "7.4"],
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        if result.returncode != 0 or not output.is_file() or output.stat().st_size == 0:
            detail = (result.stderr or result.stdout or "no diagnostic output").strip()[:500]
            raise RuntimeError(f"Open Babel receptor protonation failed: {detail}")
        yield str(output)
```

- [ ] **Step 2: Confirm red, implement the context manager, integrate it around `dock_all`, and update `--mcce` help text.**

- [ ] **Step 3: Add success/timeout cleanup tests and run docking contract tests.**

### Task 5: Build a verifiable compact snapshot

**Files:**
- Create: `scientific_validation/four_level_cli_1kx10k/verify_snapshot.py`
- Create: `validation/frozen_run/snapshot.json`
- Replace: `validation/frozen_run/MANIFEST.sha256`
- Modify: `README.md`
- Modify: `VALIDATION.md`
- Modify: `FULL_SOURCE_INVENTORY.md`
- Test: `tests/test_four_level_run_verifier.py`

- [ ] **Step 1: Write failing snapshot tests** for missing/tampered included files, cross-document count disagreement, nonzero layer failures, parity failure, malformed complete checkpoints, and valid compact evidence.

```python
def verify_snapshot(snapshot_dir: str | Path) -> dict[str, object]:
    """Verify included bytes and summary consistency without claiming row-level replay."""
```

- [ ] **Step 2: Confirm red, implement manifest parsing and summary/checkpoint/docking consistency checks.**

- [ ] **Step 3: Create `snapshot.json` using the full run's verified counts and full-manifest SHA-256; regenerate compact `MANIFEST.sha256` from included files except itself.**

- [ ] **Step 4: Verify compact snapshot succeeds and the untouched complete 1 GB run still succeeds with `verify_run --strict`.**

### Task 6: Add packaging, CI, and public-data exclusions

**Files:**
- Create: `pyproject.toml`
- Create: `scoring/__init__.py`
- Create: `scientific_validation/__init__.py`
- Create: `.github/workflows/ci.yml`
- Create: `THIRD_PARTY_DATA.md`
- Modify: local imports in `scoring/scoring.py`, `scoring/pipeline_router.py`, and `scoring/scripts/*.py`
- Modify: `scripts/build_source_release.py`
- Modify: `COMPLIANCE.md`
- Modify: `README.md`
- Modify: `FULL_SOURCE_INVENTORY.md`
- Test: `tests/test_four_level_cli_contract.py`

- [ ] **Step 1: Write failing tests** that the source builder excludes `legacy/**/*.csv|json|parquet` but retains `.py`, `.sh`, `.md`, and `.svg`; and that package-style imports expose `scoring.scoring.main`.

- [ ] **Step 2: Confirm red, add dual relative/script imports and PEP 621 entry points.**

```toml
[project.scripts]
four-level-molecule = "scoring.scoring:main"
four-level-benchmark = "scientific_validation.four_level_cli_1kx10k.batch_cli:main"
four-level-doctor = "scientific_validation.four_level_cli_1kx10k.runtime_doctor:main"
four-level-verify-snapshot = "scientific_validation.four_level_cli_1kx10k.verify_snapshot:main"
```

- [ ] **Step 3: Add CI source tests, wheel build/install, entry-point smoke, deterministic archive build, and archive-policy checks.**

- [ ] **Step 4: Document source-available status and excluded upstream-derived row-level data without asserting unverified third-party license rights.**

- [ ] **Step 5: Build wheel/sdist, install wheel into an isolated target, and exercise all four `--help` entry points.**

### Task 7: Final scientific and release verification

**Files:**
- Modify only if a verification failure proves a defect.
- Rebuild: `packages/four-level-molecule-cli-github-release-hardening-v2-20260722.zip`
- Rebuild: `packages/four-level-molecule-cli-offline-assets-release-hardening-v2-20260722.tar.gz`

- [ ] **Step 1: Run `python -m pytest -q`, asset integration tests, docking integration tests, and `pip check`.**

- [ ] **Step 2: Verify clean/tampered asset behavior and strict runtime doctor with both binary environment styles.**

- [ ] **Step 3: Run compact snapshot verification and complete 1 GB strict run verification.**

- [ ] **Step 4: Recompute the 20-molecule parity gate and one complete 10,000-row target comparison; require maximum numeric difference zero and zero status mismatches.**

- [ ] **Step 5: Rebuild archives twice and require byte-identical outputs, valid extraction, no host paths, no prohibited legacy data, no model/data assets in the GitHub ZIP, and clean asset manifest verification in the offline archive.**

- [ ] **Step 6: Publish archive sizes and SHA-256 digests; leave the original verified archives untouched.**
