import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

TESTS_DIR = Path(__file__).resolve().parent
APP_DIR = TESTS_DIR.parent / "app"
sys.path.insert(0, str(APP_DIR))
sys.path.insert(0, str(APP_DIR.parent))

import deploy_queue  # noqa: E402
import jobs  # noqa: E402
import runs  # noqa: E402


class QueueTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self._orig_data_dir = runs.DATA_DIR
        runs.DATA_DIR = Path(self.tmp.name)
        runs.ensure_layout()

    def tearDown(self):
        runs.DATA_DIR = self._orig_data_dir
        self.tmp.cleanup()

    def _mk_run(self, name):
        run = runs.create_run(name)
        p = runs.run_paths(run["run_id"])
        Path(p["references"]).mkdir(parents=True, exist_ok=True)
        Path(p["samples"]).mkdir(parents=True, exist_ok=True)
        (Path(p["references"]) / "r.fa").write_text(">r\nACGT\n")
        (Path(p["samples"]) / "s_R1.fastq.gz").write_bytes(b"")
        return run

    def test_build_spec_captures_config(self):
        run = self._mk_run("r1")
        spec = deploy_queue.build_spec(run["run_id"], 2, 4, False, True, False,
                                       config_path="/cfg.json")
        self.assertEqual(spec["run_id"], run["run_id"])
        self.assertEqual(spec["workers"], 2)
        self.assertEqual(spec["threads"], 4)
        self.assertFalse(spec["trim"])
        self.assertTrue(spec["prefilter"])
        self.assertFalse(spec["dedup"])
        self.assertIn("--no-trim", spec["map_cmd"])
        self.assertNotIn("--no-prefilter", spec["map_cmd"])
        self.assertIn("--no-dedup", spec["map_cmd"])
        self.assertIn("--config", spec["map_cmd"])
        self.assertNotIn("pathogen", spec)
        self.assertNotIn("assembly_id", spec)

    def test_start_queue_rejects_when_active(self):
        run = self._mk_run("r1")
        queue = {"status": "running", "pid": os.getpid(), "items": []}
        deploy_queue.write_queue(queue)
        with mock.patch.object(deploy_queue, "pid_alive", return_value=True):
            with self.assertRaises(RuntimeError):
                deploy_queue.start_queue([{"run_id": run["run_id"], "spec": {"run_name": "r1"}}])

    def test_request_cancel(self):
        queue = {"status": "running", "pid": 123, "items": []}
        deploy_queue.write_queue(queue)
        with mock.patch.object(deploy_queue, "pid_alive", return_value=True):
            self.assertTrue(deploy_queue.request_cancel())
        got = deploy_queue.read_queue()
        self.assertTrue(got["cancel_requested"])
        self.assertEqual(got["status"], "cancelling")


class QueueRunnerTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self._orig_data_dir = runs.DATA_DIR
        runs.DATA_DIR = Path(self.tmp.name)
        runs.ensure_layout()

    def tearDown(self):
        runs.DATA_DIR = self._orig_data_dir
        self.tmp.cleanup()

    def _mk_run(self, name):
        run = runs.create_run(name)
        p = runs.run_paths(run["run_id"])
        for k in ("references", "samples", "output"):
            Path(p[k]).mkdir(parents=True, exist_ok=True)
        return run

    def _spec(self, run_id):
        p = runs.run_paths(run_id)
        return {
            "run_id": run_id,
            "run_name": run_id,
            "map_cmd": [sys.executable, "-c", "pass"],
            "extract_script": "x",
            "refs_dir": p["references"],
            "samples_dir": p["samples"],
            "threads": 1,
        }

    def test_runner_continues_after_failure(self):
        r1 = self._mk_run("r1")["run_id"]
        r2 = self._mk_run("r2")["run_id"]
        q = {
            "status": "running", "cancel_requested": False,
            "items": [
                {"run_id": r1, "status": "pending", "spec": self._spec(r1)},
                {"run_id": r2, "status": "pending", "spec": self._spec(r2)},
            ],
        }
        qpath = runs.ensure_layout() / "q.json"
        qpath.write_text(json.dumps(q))

        calls = []

        def fake_start(item):
            calls.append(item["run_id"])
            rid = item["run_id"]
            out = runs.run_paths(rid)["output"]
            status = "failed" if rid == r1 else "done"
            jobs.write_job(out, {"run_id": rid, "status": status, "pid": 1, "stage": "extract",
                                 "returncode": 1 if status == "failed" else 0})
            return status, jobs.read_job(out)

        with mock.patch("queue_runner._start_item", side_effect=fake_start), \
             mock.patch("queue_runner._wait_for_terminal", side_effect=lambda rid, out: _job_status(out)), \
             mock.patch("queue_runner._wait_for_lock_release", return_value=True), \
             mock.patch("jobs.pid_alive", return_value=False), \
             mock.patch("runs.active_run", return_value=None):
            rc = __import__("queue_runner").main(["queue_runner.py", str(qpath)])

        self.assertEqual(rc, 0)
        got = deploy_queue.read_queue()
        self.assertEqual(calls, [r1, r2])
        self.assertEqual(got["items"][0]["status"], "failed")
        self.assertEqual(got["items"][1]["status"], "done")
        self.assertEqual(got["status"], "done")

    def test_runner_stops_before_next_on_cancel(self):
        r1 = self._mk_run("r1")["run_id"]
        r2 = self._mk_run("r2")["run_id"]
        q = {
            "status": "running", "cancel_requested": False,
            "items": [
                {"run_id": r1, "status": "pending", "spec": self._spec(r1)},
                {"run_id": r2, "status": "pending", "spec": self._spec(r2)},
            ],
        }
        qpath = runs.ensure_layout() / "q.json"
        qpath.write_text(json.dumps(q))

        def fake_start(item):
            rid = item["run_id"]
            out = runs.run_paths(rid)["output"]
            jobs.write_job(out, {"run_id": rid, "status": "done", "pid": 1, "stage": "extract", "returncode": 0})
            qcur = deploy_queue.read_queue()
            qcur["cancel_requested"] = True
            deploy_queue.write_queue(qcur)
            return "done", jobs.read_job(out)

        with mock.patch("queue_runner._start_item", side_effect=fake_start), \
             mock.patch("queue_runner._wait_for_terminal", side_effect=lambda rid, out: _job_status(out)), \
             mock.patch("queue_runner._wait_for_lock_release", return_value=True), \
             mock.patch("jobs.pid_alive", return_value=False), \
             mock.patch("runs.active_run", return_value=None):
            rc = __import__("queue_runner").main(["queue_runner.py", str(qpath)])

        self.assertEqual(rc, 0)
        got = deploy_queue.read_queue()
        self.assertEqual(got["items"][0]["status"], "done")
        self.assertEqual(got["items"][1]["status"], "skipped")
        self.assertEqual(got["status"], "cancelled")


def _job_status(out):
    job = jobs.read_job(out) or {}
    return job.get("status", ""), job


if __name__ == "__main__":
    unittest.main()
