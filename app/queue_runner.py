"""Detached runner for sequential deployment queue.

Usage: python queue_runner.py <deploy_queue_spec.json>

Runs queue items in order. Each item starts a normal mapping job through
jobs.start_job() and the runner polls until terminal. On failure/cancel, it
continues to next item (continue_on_error=true by design). Queue cancel stops
before starting the next item; it does not kill the current job.
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(APP_DIR))
sys.path.insert(0, str(APP_DIR.parent))

import deploy_queue  # noqa: E402
import jobs  # noqa: E402
import runs  # noqa: E402


POLL_SECONDS = 2.0
START_WAIT_SECONDS = 300.0


def _now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _log(msg):
    print(f"[queue_runner] {msg}", flush=True)


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def _save(queue):
    deploy_queue.write_queue(queue)
    return queue


def _item_terminal(status):
    return status in deploy_queue.ITEM_TERMINAL


def _job_status(output_dir):
    job = jobs.read_job(output_dir) or {}
    return job.get("status", ""), job


def _wait_for_lock_release(run_id, timeout=120.0):
    deadline = time.time() + timeout
    while time.time() < deadline:
        active = runs.active_run()
        if active != run_id:
            return True
        time.sleep(POLL_SECONDS)
    return runs.active_run() != run_id


def _wait_for_terminal(run_id, output_dir):
    while True:
        status, job = _job_status(output_dir)
        active = runs.active_run()
        if status in {"done", "failed", "cancelled"} and active != run_id:
            return status, job
        # If job file is missing but lock released and we had started, treat as failed.
        if not job and active != run_id:
            return "failed", {"error": "job state missing"}
        time.sleep(POLL_SECONDS)


def _start_item(item):
    spec = dict(item["spec"])
    run_id = item["run_id"]
    output_dir = runs.run_paths(run_id)["output"]
    # Rebuild command/path-sensitive spec from current paths but keep captured config.
    spec["samples_dir"] = runs.run_paths(run_id)["samples"]
    spec["output_dir"] = output_dir
    jobs.start_job(spec, output_dir, run_id=run_id)
    deadline = time.time() + START_WAIT_SECONDS
    while time.time() < deadline:
        status, job = _job_status(output_dir)
        if status:
            return status, job
        if runs.active_run() == run_id:
            return "running", jobs.read_job(output_dir) or {}
        time.sleep(1.0)
    raise RuntimeError(f"timed out waiting for run {run_id} to start")


def main(argv):
    queue = _load(argv[1])
    queue["status"] = "running"
    queue["pid"] = os.getpid()
    queue["started"] = queue.get("started") or _now()
    queue["finished"] = None
    _save(queue)
    try:
        items = queue.get("items", [])
        for idx, item in enumerate(items):
            queue = deploy_queue.read_queue() or queue
            if queue.get("cancel_requested"):
                _log("cancel requested; stopping before next item")
                break
            if item.get("status") in deploy_queue.ITEM_TERMINAL:
                continue

            run_id = item["run_id"]
            queue["current_index"] = idx
            queue["current_run_id"] = run_id
            deploy_queue.update_item(queue, idx, status="running", started=_now(), error="")
            _save(queue)
            _log(f"starting {idx + 1}/{len(items)}: {run_id}")

            try:
                record = runs.read_run(run_id)
                if record is None:
                    raise RuntimeError(f"unknown run: {run_id}")
                status, job = _start_item(item)
                if status not in {"done", "failed", "cancelled"}:
                    status, job = _wait_for_terminal(run_id, runs.run_paths(run_id)["output"])
                _wait_for_lock_release(run_id)
                queue = deploy_queue.read_queue() or queue
                deploy_queue.update_item(
                    queue, idx,
                    status=status if status in deploy_queue.ITEM_TERMINAL else "failed",
                    finished=_now(),
                    returncode=job.get("returncode"),
                    error=job.get("error", ""),
                    job_status=status,
                )
                _save(queue)
                _log(f"finished {run_id}: {status}")
            except Exception as e:
                queue = deploy_queue.read_queue() or queue
                deploy_queue.update_item(queue, idx, status="failed", finished=_now(),
                                         error=f"{type(e).__name__}: {e}")
                _save(queue)
                _log(f"item failed {run_id}: {type(e).__name__}: {e}")
                _wait_for_lock_release(run_id)
                continue

        queue = deploy_queue.read_queue() or queue
        if queue.get("cancel_requested"):
            # Mark remaining pending as skipped.
            for idx, item in enumerate(queue.get("items", [])):
                if item.get("status") == "pending":
                    deploy_queue.update_item(queue, idx, status="skipped", finished=_now(),
                                             error="queue cancelled")
            queue["status"] = "cancelled"
        else:
            statuses = [i.get("status") for i in queue.get("items", [])]
            queue["status"] = "done" if statuses and all(s in deploy_queue.ITEM_TERMINAL for s in statuses) else "failed"
        queue["finished"] = _now()
        queue["current_run_id"] = None
        _save(queue)
        return 0
    except Exception as e:
        queue = deploy_queue.read_queue() or queue
        queue["status"] = "failed"
        queue["finished"] = _now()
        queue["error"] = f"{type(e).__name__}: {e}"
        _save(queue)
        _log(f"queue failed: {type(e).__name__}: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
