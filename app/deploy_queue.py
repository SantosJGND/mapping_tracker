"""Persistent sequential deployment queue.

A queue runs one mapping+extract job at a time, in order. Each item carries its
own deploy spec snapshot captured at enqueue time. The queue runner is detached,
so it survives Streamlit reruns. Job execution still uses the existing global
run lock and pipeline runner.
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional

APP_DIR = Path(__file__).resolve().parent
_PARENT = str(APP_DIR.parent)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

import jobs  # noqa: E402
import runs  # noqa: E402

QUEUE_FILE = "deploy_queue.json"
QUEUE_LOG = "deploy_queue.log"
QUEUE_RUNNER = APP_DIR / "queue_runner.py"

QUEUE_TERMINAL = {"done", "failed", "cancelled"}
ITEM_TERMINAL = {"done", "failed", "cancelled", "skipped"}
ITEM_RUNNING = {"running"}


def queue_path() -> Path:
    return runs.ensure_layout() / QUEUE_FILE


def queue_log_path() -> Path:
    return runs.ensure_layout() / QUEUE_LOG


def _atomic_write(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
    os.replace(tmp, path)


def read_queue() -> Optional[dict]:
    try:
        with open(queue_path()) as fh:
            data = json.load(fh)
        return data if isinstance(data, dict) else None
    except (OSError, ValueError):
        return None


def write_queue(data: dict) -> dict:
    _atomic_write(queue_path(), data)
    return data


def update_queue(**fields) -> Optional[dict]:
    queue = read_queue()
    if queue is None:
        return None
    queue.update(fields)
    return write_queue(queue)


def pid_alive(pid: Optional[int]) -> bool:
    return jobs.pid_alive(pid)


def queue_active(queue: Optional[dict] = None) -> bool:
    queue = read_queue() if queue is None else queue
    return bool(queue) and queue.get("status") in {"running", "cancelling"} and pid_alive(queue.get("pid"))


def _now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _item(run_id: str, spec: dict) -> dict:
    return {
        "item_id": uuid.uuid4().hex[:8],
        "run_id": run_id,
        "run_name": spec.get("run_name") or run_id,
        "status": "pending",
        "started": None,
        "finished": None,
        "error": "",
        "returncode": None,
        "job_status": "",
        "spec": spec,
    }


def build_spec(run_id: str, workers: int, threads: int, trim: bool,
               prefilter: bool, dedup: bool, config_path: Optional[str] = None) -> dict:
    paths = runs.run_paths(run_id)
    record = runs.read_run(run_id) or {}
    map_cmd = jobs.build_map_cmd(
        str(jobs.DEFAULT_MAP_SCRIPT), paths["references"], paths["samples"],
        paths["tsv"], paths["output"], paths["log"], int(workers), int(threads),
        bool(trim), bool(prefilter), bool(dedup), config_path=config_path)
    return {
        "run_id": run_id,
        "run_name": record.get("name", run_id),
        "map_cmd": map_cmd,
        "extract_script": str(jobs.DEFAULT_EXTRACT_SCRIPT),
        "refs_dir": paths["references"],
        "samples_dir": paths["samples"],
        "threads": int(threads),
        "workers": int(workers),
        "trim": bool(trim),
        "prefilter": bool(prefilter),
        "dedup": bool(dedup),
        "config_path": config_path,
    }


def start_queue(items: List[dict]) -> dict:
    """Start a detached sequential deployment queue.

    items: list of {"run_id": ..., "spec": ...}. Spec must be a snapshot created
    by build_spec() or equivalent.
    """
    if not items:
        raise ValueError("No runs selected for deployment")
    current = read_queue()
    if queue_active(current):
        raise RuntimeError("A deployment queue is already running")
    if runs.active_run():
        raise RuntimeError(f"A job is already running in run {runs.active_run()}")

    queue = {
        "queue_id": uuid.uuid4().hex[:10],
        "status": "running",
        "pid": None,
        "created": _now(),
        "started": _now(),
        "finished": None,
        "current_index": -1,
        "current_run_id": None,
        "continue_on_error": True,
        "cancel_requested": False,
        "items": [_item(it["run_id"], it["spec"]) for it in items],
    }
    spec_path = runs.ensure_layout() / "deploy_queue_spec.json"
    _atomic_write(spec_path, queue)
    write_queue(queue)
    with open(queue_log_path(), "w") as log:
        proc = subprocess.Popen(
            [sys.executable, str(QUEUE_RUNNER), str(spec_path)],
            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
            start_new_session=True)
    queue["pid"] = proc.pid
    return write_queue(queue)


def request_cancel() -> bool:
    queue = read_queue()
    if not queue_active(queue):
        return False
    queue["cancel_requested"] = True
    queue["status"] = "cancelling"
    write_queue(queue)
    return True


def cancel_queue_runner() -> bool:
    """Force-stop the queue runner process. Does not cancel the current job."""
    queue = read_queue()
    if not queue_active(queue):
        return False
    request_cancel()
    try:
        os.killpg(queue["pid"], signal.SIGTERM)
    except (ProcessLookupError, KeyError, TypeError):
        pass
    return True


def queue_log_tail(lines: int = 40) -> str:
    try:
        with open(queue_log_path()) as fh:
            return "".join(fh.readlines()[-lines:])
    except OSError:
        return ""


def update_item(queue: dict, index: int, **fields) -> dict:
    items = list(queue.get("items", []))
    if index < 0 or index >= len(items):
        return queue
    item = dict(items[index])
    item.update(fields)
    items[index] = item
    queue["items"] = items
    return queue
