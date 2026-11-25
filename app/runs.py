from __future__ import annotations

import json
import os
import time
import uuid
from pathlib import Path
from typing import Dict, List, Optional

APP_DIR = Path(__file__).resolve().parent
DATA_DIR = Path(os.environ.get("MAPPING_TRACKER_DATA_DIR", str(APP_DIR.parent / "data")))

INDEX_NAME = "index.json"
LOCK_NAME = ".active.lock"
TERMINAL_STATUSES = {"done", "failed", "cancelled", "invalid"}


def ensure_layout() -> Path:
    runs_dir = Path(DATA_DIR) / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)
    (Path(DATA_DIR) / "tmp").mkdir(parents=True, exist_ok=True)
    return runs_dir


def _index_path() -> Path:
    return ensure_layout() / INDEX_NAME


def _lock_path() -> Path:
    return ensure_layout() / LOCK_NAME


def _run_dir(run_id: str) -> Path:
    if not run_id or Path(run_id).name != run_id or run_id in {".", ".."}:
        raise ValueError(f"Invalid run id: {run_id}")
    return ensure_layout() / run_id


def _atomic_write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
    os.replace(tmp, path)


def _read_json(path: Path) -> Optional[dict]:
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def new_run_id() -> str:
    return time.strftime("%Y%m%d_%H%M%S") + "_" + uuid.uuid4().hex[:6]


def read_index() -> List[dict]:
    data = _read_json(_index_path())
    return data if isinstance(data, list) else []


def write_index(entries: List[dict]) -> None:
    _atomic_write_json(_index_path(), entries)


def create_run(name: str = "", reference_count: int = 0, sample_count: int = 0) -> dict:
    run_id = new_run_id()
    run_dir = _run_dir(run_id)
    paths = run_paths(run_id)
    for key in ("references", "samples", "output", "extract"):
        Path(paths[key]).mkdir(parents=True, exist_ok=True)
    run = {
        "run_id": run_id,
        "name": name or run_id,
        "status": "ready",
        "stage": "uploaded",
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "reference_count": reference_count,
        "sample_count": sample_count,
        "pid": None,
        "error": "",
    }
    write_run(run_id, run)
    return run


def read_run(run_id: str) -> Optional[dict]:
    return _read_json(_run_dir(run_id) / "run.json")


def write_run(run_id: str, data: dict) -> dict:
    run_dir = _run_dir(run_id)
    record = dict(data)
    record["run_id"] = run_id
    _atomic_write_json(run_dir / "run.json", record)

    entries = read_index()
    for i, entry in enumerate(entries):
        if entry.get("run_id") == run_id:
            entries[i] = record
            break
    else:
        entries.append(record)
    write_index(entries)
    return record


def update_run(run_id: str, **fields) -> Optional[dict]:
    record = read_run(run_id)
    if record is None:
        return None
    record.update(fields)
    return write_run(run_id, record)


def list_runs() -> List[dict]:
    entries = [e for e in read_index() if _run_exists(str(e.get("run_id", "")))]
    if not entries:
        for run_dir in sorted(_run_dir("x").parent.iterdir()):
            if not run_dir.is_dir():
                continue
            record = _read_json(run_dir / "run.json")
            if record and record.get("run_id"):
                entries.append(record)
    return sorted(entries, key=lambda r: str(r.get("created", "")), reverse=True)


def _run_exists(run_id: str) -> bool:
    try:
        return (_run_dir(run_id) / "run.json").is_file()
    except ValueError:
        return False


def latest_run() -> Optional[dict]:
    runs = list_runs()
    return runs[0] if runs else None


def run_paths(run_id: str) -> Dict[str, str]:
    run_dir = _run_dir(run_id)
    output = run_dir / "output"
    return {
        "dir": str(run_dir),
        "run_json": str(run_dir / "run.json"),
        "job_json": str(output / "job.json"),
        "job_log": str(run_dir / "job.log"),
        "references": str(run_dir / "references"),
        "samples": str(run_dir / "samples"),
        "output": str(output),
        "extract": str(output / "extract"),
        "tsv": str(output / "mapping_statistics.tsv"),
        "log": str(output / "mapping_log.txt"),
    }


def run_log_tail(run_id: str, lines: int = 80) -> str:
    try:
        with open(run_paths(run_id)["job_log"]) as fh:
            return "".join(fh.readlines()[-lines:])
    except OSError:
        return ""


def index_signature() -> tuple:
    parts = []
    for path in (_index_path(), _lock_path()):
        try:
            parts.append((str(path), path.stat().st_mtime_ns))
        except OSError:
            parts.append((str(path), None))
    return tuple(parts)


def _pid_alive(pid: object) -> bool:
    if not isinstance(pid, int) or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    try:
        with open(f"/proc/{pid}/stat") as fh:
            state = fh.read().rsplit(")", 1)[1].split()[0]
        return state != "Z"
    except (OSError, IndexError):
        return True


def active_run() -> Optional[str]:
    lock = _read_json(_lock_path())
    if not isinstance(lock, dict):
        return None
    run_id = lock.get("run_id")
    return run_id if isinstance(run_id, str) and _pid_alive(lock.get("pid")) else None


def acquire_lock(run_id: str) -> Optional[str]:
    if active_run():
        return None
    path = _lock_path()
    if path.exists():
        path.unlink()
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return None
    with os.fdopen(fd, "w") as fh:
        json.dump({"run_id": run_id, "pid": os.getpid(),
                   "started": time.strftime("%Y-%m-%d %H:%M:%S")}, fh)
    return run_id


def release_lock(run_id: str) -> bool:
    path = _lock_path()
    if not path.exists():
        return False
    lock = _read_json(path)
    if isinstance(lock, dict) and lock.get("run_id") == run_id:
        path.unlink(missing_ok=True)
        return True
    return False
