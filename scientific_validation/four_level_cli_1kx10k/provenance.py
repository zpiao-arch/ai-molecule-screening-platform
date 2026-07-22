from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def canonical_json(payload: Any) -> str:
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def payload_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _checkpoint_path(run_dir: str | Path, stage: str) -> Path:
    return Path(run_dir) / "checkpoints" / f"{stage}.json"


def write_stage_checkpoint(
    run_dir: str | Path,
    stage: str,
    *,
    inputs: dict[str, Any],
    status: str,
    outputs: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> Path:
    if status not in {"in_progress", "complete", "failed"}:
        raise ValueError("checkpoint status must be in_progress, complete, or failed")
    destination = _checkpoint_path(run_dir, stage)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": stage,
        "status": status,
        "input_fingerprint": payload_sha256(inputs),
        "inputs": inputs,
        "outputs": outputs or {},
        "metadata": metadata or {},
        "updated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    temporary = destination.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    os.replace(temporary, destination)
    event = {
        "checkpoint": str(destination.relative_to(Path(run_dir))),
        "event": "stage_checkpoint",
        "input_fingerprint": payload["input_fingerprint"],
        "stage": stage,
        "status": status,
        "timestamp_utc": payload["updated_at_utc"],
    }
    log_path = Path(run_dir) / "run.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(canonical_json(event) + "\n")
    return destination


def fail_stage_checkpoint(
    run_dir: str | Path,
    stage: str,
    error: BaseException,
    *,
    max_message_length: int = 1000,
) -> Path | None:
    """Atomically turn an existing in-progress checkpoint into a failed checkpoint."""
    destination = _checkpoint_path(run_dir, stage)
    if not destination.is_file():
        return None
    payload = json.loads(destination.read_text(encoding="utf-8"))
    timestamp = datetime.now(timezone.utc).isoformat()
    message = str(error).strip() or type(error).__name__
    payload["status"] = "failed"
    payload["updated_at_utc"] = timestamp
    metadata = dict(payload.get("metadata") or {})
    metadata.update(
        {
            "error_type": type(error).__name__,
            "error_message": message[:max_message_length],
            "failed_at_utc": timestamp,
        }
    )
    payload["metadata"] = metadata
    temporary = destination.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True), encoding="utf-8")
    os.replace(temporary, destination)
    event = {
        "checkpoint": str(destination.relative_to(Path(run_dir))),
        "event": "stage_checkpoint",
        "input_fingerprint": payload["input_fingerprint"],
        "stage": stage,
        "status": "failed",
        "timestamp_utc": timestamp,
    }
    log_path = Path(run_dir) / "run.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(canonical_json(event) + "\n")
    return destination


def validate_stage_checkpoint(
    run_dir: str | Path,
    stage: str,
    *,
    inputs: dict[str, Any],
    require_complete: bool = False,
) -> dict[str, Any]:
    path = _checkpoint_path(run_dir, stage)
    if not path.is_file():
        return {"ok": False, "reason": "missing", "path": str(path)}
    payload = json.loads(path.read_text(encoding="utf-8"))
    actual = payload_sha256(inputs)
    expected = payload.get("input_fingerprint")
    reasons = []
    if actual != expected:
        reasons.append("input fingerprint mismatch")
    if require_complete and payload.get("status") != "complete":
        reasons.append("checkpoint is not complete")
    return {
        "ok": not reasons,
        "reason": "; ".join(reasons) if reasons else "ok",
        "path": str(path),
        "checkpoint": payload,
        "actual_input_fingerprint": actual,
    }


def require_stage_checkpoint(
    run_dir: str | Path,
    stage: str,
    *,
    inputs: dict[str, Any],
    require_complete: bool = False,
) -> dict[str, Any]:
    result = validate_stage_checkpoint(
        run_dir,
        stage,
        inputs=inputs,
        require_complete=require_complete,
    )
    if not result["ok"]:
        raise RuntimeError(f"{stage} checkpoint {result['reason']}")
    return result
