import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "app"
sys.path.insert(0, str(APP))
import runs  # noqa: E402


class RunTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="runs_test_")
        self.old_data_dir = runs.DATA_DIR
        runs.DATA_DIR = os.path.join(self.tmp, "data")

    def tearDown(self):
        runs.DATA_DIR = self.old_data_dir
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_create_update_and_list(self):
        run = runs.create_run("first", 2, 3)
        self.assertTrue(Path(runs.run_paths(run["run_id"])["dir"]).is_dir())
        runs.update_run(run["run_id"], status="running", stage="mapping")
        listed = runs.list_runs()
        self.assertEqual(len(listed), 1)
        self.assertEqual(listed[0]["status"], "running")
        self.assertEqual(runs.read_run(run["run_id"])["stage"], "mapping")

    def test_run_paths_are_isolated(self):
        a = runs.create_run("a")
        b = runs.create_run("b")
        paths_a = runs.run_paths(a["run_id"])
        paths_b = runs.run_paths(b["run_id"])
        self.assertNotEqual(paths_a["dir"], paths_b["dir"])
        for key in ("references", "samples", "output", "extract"):
            self.assertTrue(Path(paths_a[key]).is_dir())

    def test_global_lock(self):
        run = runs.create_run("locked")
        run_id = run["run_id"]
        self.assertEqual(runs.acquire_lock(run_id), run_id)
        self.assertEqual(runs.active_run(), run_id)
        self.assertIsNone(runs.acquire_lock(runs.create_run("other")["run_id"]))
        self.assertTrue(runs.release_lock(run_id))
        self.assertIsNone(runs.active_run())

    def test_stale_lock_is_replaced(self):
        proc = subprocess.Popen(["true"])
        proc.wait()
        path = runs._lock_path()
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            json.dump({"run_id": "old", "pid": proc.pid}, fh)
        self.assertIsNone(runs.active_run())
        run_id = runs.create_run("new")["run_id"]
        self.assertEqual(runs.acquire_lock(run_id), run_id)


if __name__ == "__main__":
    unittest.main()
