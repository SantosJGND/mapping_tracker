"""Detached runner for one mapping + extract job (started by jobs.start_job).

Usage: python pipeline_runner.py <job_spec.json>

Stages: mapping (map_references.py) -> registering samples (manifest) -> extract
(extract_remap.py). Progress and the final status are written to job.json in the
output dir. Child process output goes to job.log (inherited stdout/stderr).
"""
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jobs  # noqa: E402
import runs  # noqa: E402


def _now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _cancelled(out):
    job = jobs.read_job(out) or {}
    return job.get("status") == "cancelled"


def main(argv):
    with open(argv[1]) as fh:
        spec = json.load(fh)
    out = spec["output_dir"]
    run_id = spec.get("run_id")
    try:
        try:
            jobs.update_job(out, status="running", stage="mapping")
            rc = subprocess.run(spec["map_cmd"]).returncode
            if rc != 0:
                if not _cancelled(out):
                    jobs.update_job(out, status="failed", stage="mapping", returncode=rc,
                                    finished=_now())
                return rc

            jobs.update_job(out, stage="registering samples")
            manifest, rows = jobs.build_manifest(
                spec["samples_dir"], out, spec["pathogen"], spec["assembly_id"],
                spec["extract_dir"])
            names = [r["name"] for r in rows]
            print(f"[pipeline_runner] registered {len(names)} sample(s): {', '.join(names)}",
                  flush=True)

            jobs.update_job(out, stage="extract", samples=names, manifest=manifest)
            cmd = jobs.build_extract_cmd(spec["extract_script"], spec["ref_fasta"], manifest,
                                         spec["extract_dir"], spec["threads"])
            rc = subprocess.run(cmd).returncode
            if not _cancelled(out):
                jobs.update_job(out, status="done" if rc == 0 else "failed", stage="extract",
                                returncode=rc, finished=_now())
            return rc
        except Exception as e:
            jobs.update_job(out, status="failed", error=f"{type(e).__name__}: {e}",
                            finished=_now())
            print(f"[pipeline_runner] failed: {type(e).__name__}: {e}", flush=True)
            return 1
    finally:
        if run_id:
            runs.release_lock(run_id)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
