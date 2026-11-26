"""Tests for app/jobs.py and app/pipeline_runner.py. Run from the mapping_tracker directory:

    python -m unittest discover -s tests -v

Runner tests use stub map/extract scripts (no bwa or samtools needed).
"""
import json
import os
import shutil
import time
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "app"
sys.path.insert(0, str(APP))
import jobs  # noqa: E402
import runs  # noqa: E402

RUNNER = APP / "pipeline_runner.py"

STUB_MAP = """import os, sys
out = sys.argv[sys.argv.index('--output-dir') + 1]
open(os.path.join(os.environ['STUB_LOG']), 'a').write('map ' + ' '.join(sys.argv[1:]) + '\\n')
if os.environ.get('STUB_MAP_FAIL'):
    sys.exit(3)
for sample in ('s1',):
    d = os.path.join(out, 'PathoA')
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, sample + '_GCF_1.bam'), 'w').write('bam')
"""

STUB_EXTRACT = """import os, sys
outdir = sys.argv[sys.argv.index('--outdir') + 1]
os.makedirs(outdir, exist_ok=True)
open(os.path.join(os.environ['STUB_LOG']), 'a').write('extract ' + ' '.join(sys.argv[1:]) + '\\n')
open(os.path.join(outdir, 'subset_ref.fa.fai'), 'w').write('c\\t10\\t0\\t10\\t11\\n')
"""


class JobsTestCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="jobs_test_")
        self.old_data_dir = runs.DATA_DIR
        self.old_data_env = os.environ.get("MAPPING_TRACKER_DATA_DIR")
        runs.DATA_DIR = os.path.join(self.tmp, "data")
        os.environ["MAPPING_TRACKER_DATA_DIR"] = str(runs.DATA_DIR)
        self.samples = os.path.join(self.tmp, "samples")
        self.out = os.path.join(self.tmp, "out")
        os.makedirs(self.samples)
        for name in ("s1_1P.fastq.gz", "s1_2P.fastq.gz", "s2_1P.fastq.gz", "s2_2P.fastq.gz"):
            Path(self.samples, name).write_text("x")
        os.makedirs(self.out)

    def tearDown(self):
        runs.DATA_DIR = self.old_data_dir
        if self.old_data_env is None:
            os.environ.pop("MAPPING_TRACKER_DATA_DIR", None)
        else:
            os.environ["MAPPING_TRACKER_DATA_DIR"] = self.old_data_env
        shutil.rmtree(self.tmp, ignore_errors=True)


class CommandTests(JobsTestCase):
    def test_map_cmd_flags(self):
        cmd = jobs.build_map_cmd("/m.py", "refs", "samples", "t.tsv", "out", "log.txt",
                                 2, 8, trim=False, prefilter=True, dedup=False,
                                 config_path="/config.json")
        self.assertEqual(cmd[1], "/m.py")
        self.assertIn("--no-trim", cmd)
        self.assertNotIn("--no-prefilter", cmd)
        self.assertIn("--no-dedup", cmd)
        self.assertEqual(cmd[cmd.index("--workers") + 1], "2")
        self.assertEqual(cmd[cmd.index("--threads") + 1], "8")
        self.assertEqual(cmd[cmd.index("--config") + 1], "/config.json")

    def test_extract_cmd(self):
        cmd = jobs.build_extract_cmd("/e.py", "/r.fa", "/m.tsv", "/o", 4)
        self.assertEqual(cmd[1:], ["/e.py", "--ref", "/r.fa", "--manifest", "/m.tsv",
                                   "--outdir", "/o", "--threads", "4"])

    def test_extract_dir_layout(self):
        self.assertEqual(jobs.extract_dir("/out", "PathoA"),
                         os.path.join("/out", "extract", "PathoA"))


class RegistrationTests(JobsTestCase):
    def test_registers_only_samples_with_bam(self):
        # s1 has a mapped BAM for GCF_1; s2 does not
        os.makedirs(os.path.join(self.out, "PathoA"))
        Path(self.out, "PathoA", "s1_GCF_1.bam").write_text("bam")
        ext = jobs.extract_dir(self.out, "PathoA")
        manifest, rows = jobs.build_manifest(self.samples, self.out, "PathoA", "GCF_1", ext)
        self.assertEqual([r["name"] for r in rows], ["s1"])
        self.assertEqual(rows[0]["r2"], os.path.abspath(os.path.join(self.samples, "s1_2P.fastq.gz")))
        registered = jobs.read_registered(ext)
        self.assertEqual([r["name"] for r in registered], ["s1"])
        self.assertTrue(os.path.isabs(registered[0]["bam"]))

    def test_no_bam_raises(self):
        with self.assertRaises(ValueError):
            jobs.build_manifest(self.samples, self.out, "PathoA", "GCF_1",
                                jobs.extract_dir(self.out, "PathoA"))

    def test_single_end_uses_dash(self):
        Path(self.samples, "se_only.fastq.gz").write_text("x")
        os.makedirs(os.path.join(self.out, "PathoA"))
        Path(self.out, "PathoA", "se_only_GCF_1.bam").write_text("bam")
        _, rows = jobs.build_manifest(self.samples, self.out, "PathoA", "GCF_1",
                                      jobs.extract_dir(self.out, "PathoA"))
        self.assertEqual(rows[0]["name"], "se_only")
        self.assertEqual(rows[0]["r2"], "-")


class JobStateTests(JobsTestCase):
    def test_lock_refuses_second_job(self):
        jobs.write_job(self.out, {"status": "running", "pid": os.getpid()})
        self.assertTrue(jobs.job_active(self.out))
        with self.assertRaises(RuntimeError):
            jobs.start_job({"pathogen": "P", "assembly_id": "A"}, self.out)

    def test_dead_pid_is_not_active(self):
        proc = subprocess.Popen(["true"])
        proc.wait()
        jobs.write_job(self.out, {"status": "running", "pid": proc.pid})
        self.assertFalse(jobs.job_active(self.out))

    def test_finished_job_is_not_active(self):
        jobs.write_job(self.out, {"status": "done", "pid": os.getpid()})
        self.assertFalse(jobs.job_active(self.out))


class GlobalLockTests(JobsTestCase):
    def test_start_refuses_when_another_run_is_active(self):
        run = runs.create_run("locked")
        runs.acquire_lock(run["run_id"])
        spec = {"pathogen": "P", "assembly_id": "A"}
        try:
            with self.assertRaises(RuntimeError):
                jobs.start_job(spec, self.out)
        finally:
            runs.release_lock(run["run_id"])

    def test_start_refuses_when_lock_is_held(self):
        run = runs.create_run("locked")
        runs.acquire_lock(run["run_id"])
        spec = {"pathogen": "P", "assembly_id": "A"}
        try:
            with self.assertRaises(RuntimeError):
                jobs.start_job(spec, self.out, run_id=runs.create_run("other")["run_id"])
        finally:
            runs.release_lock(run["run_id"])


class CancelTests(JobsTestCase):
    def test_cancel_stops_job_and_children(self):
        sleeper = os.path.join(self.tmp, "sleeper.py")
        Path(sleeper).write_text("import time\ntime.sleep(120)\n")
        spec = {"map_cmd": [sys.executable, sleeper], "extract_script": "/none",
                "ref_fasta": "/none", "samples_dir": self.samples,
                "pathogen": "PathoA", "assembly_id": "GCF_1", "threads": 1}
        jobs.start_job(spec, self.out)
        for _ in range(40):
            job = jobs.read_job(self.out)
            if job.get("stage") == "mapping" and jobs.job_active(self.out):
                break
            time.sleep(0.25)
        self.assertTrue(jobs.job_active(self.out))
        self.assertTrue(jobs.cancel_job(self.out))
        job = jobs.read_job(self.out)
        self.assertEqual(job["status"], "cancelled")
        for _ in range(40):
            if not jobs.pid_alive(job["pid"]):
                break
            time.sleep(0.25)
        self.assertFalse(jobs.pid_alive(job["pid"]))
        self.assertEqual(jobs.read_job(self.out)["status"], "cancelled")
        self.assertFalse(jobs.cancel_job(self.out))


class RunnerTests(JobsTestCase):
    def _spec(self, map_script, extract_script, **extra):
        log = os.path.join(self.tmp, "calls.log")
        os.environ["STUB_LOG"] = log
        spec = {
            "output_dir": self.out,
            "map_cmd": [sys.executable, map_script, "--output-dir", self.out],
            "extract_script": extract_script,
            "ref_fasta": "/ref.fa",
            "samples_dir": self.samples,
            "pathogen": "PathoA",
            "assembly_id": "GCF_1",
            "extract_dir": jobs.extract_dir(self.out, "PathoA"),
            "threads": 2,
        }
        spec.update(extra)
        spec_path = os.path.join(self.tmp, "spec.json")
        with open(spec_path, "w") as fh:
            json.dump(spec, fh)
        return spec_path, log

    def _write_stubs(self):
        m = os.path.join(self.tmp, "stub_map.py")
        e = os.path.join(self.tmp, "stub_extract.py")
        Path(m).write_text(STUB_MAP)
        Path(e).write_text(STUB_EXTRACT)
        return m, e

    def test_success_runs_map_register_extract(self):
        m, e = self._write_stubs()
        spec_path, log = self._spec(m, e)
        rc = subprocess.run([sys.executable, str(RUNNER), spec_path]).returncode
        self.assertEqual(rc, 0)
        job = jobs.read_job(self.out)
        self.assertEqual(job["status"], "done")
        self.assertEqual(job["samples"], ["s1"])
        calls = Path(log).read_text().splitlines()
        self.assertTrue(calls[0].startswith("map "))
        self.assertTrue(calls[1].startswith("extract "))
        self.assertIn("--manifest", calls[1])
        registered = jobs.read_registered(jobs.extract_dir(self.out, "PathoA"))
        self.assertEqual([r["name"] for r in registered], ["s1"])

    def test_success_releases_run_lock_and_updates_run(self):
        m, e = self._write_stubs()
        run = runs.create_run("run with lock")
        run_id = run["run_id"]
        self.assertEqual(runs.acquire_lock(run_id), run_id)
        jobs.write_job(self.out, {"run_id": run_id})
        spec_path, log = self._spec(m, e, run_id=run_id)
        rc = subprocess.run([sys.executable, str(RUNNER), spec_path]).returncode
        self.assertEqual(rc, 0)
        self.assertIsNone(runs.active_run())
        self.assertEqual(runs.read_run(run_id)["status"], "done")

    def test_mapping_failure_skips_extract(self):
        m, e = self._write_stubs()
        spec_path, log = self._spec(m, e)
        os.environ["STUB_MAP_FAIL"] = "1"
        try:
            rc = subprocess.run([sys.executable, str(RUNNER), spec_path]).returncode
        finally:
            os.environ.pop("STUB_MAP_FAIL", None)
        self.assertEqual(rc, 3)
        job = jobs.read_job(self.out)
        self.assertEqual((job["status"], job["stage"]), ("failed", "mapping"))
        calls = Path(log).read_text().splitlines()
        self.assertFalse(any(c.startswith("extract") for c in calls))


if __name__ == "__main__":
    unittest.main()
