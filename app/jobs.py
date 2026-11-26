"""Launch and track map_references.py -> extract_remap.py runs from the Streamlit app.

A job is one mapping run followed by one extract run for a chosen reference:
1. map_references.py maps every sample against every reference (output: <output-dir>).
2. The samples with a mapped BAM for the chosen reference are registered in
   <output-dir>/extract/<pathogen>/manifest.tsv (this file is the record of which
   samples went into the extract).
3. extract_remap.py builds the excised reference and remaps those samples into
   <output-dir>/extract/<pathogen>/, which the Browser tab reads.

The job runs as a detached process (pipeline_runner.py), so it survives Streamlit
reruns. Its state is kept in <output-dir>/job.json; only one job runs per output dir.
"""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

APP_DIR = Path(__file__).resolve().parent
_PARENT = str(APP_DIR.parent)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from map_references import find_samples  # noqa: E402
import runs  # noqa: E402

DEFAULT_MAP_SCRIPT = APP_DIR.parent / "map_references.py"
DEFAULT_EXTRACT_SCRIPT = APP_DIR.parent.parent / "mapped_reference_extract" / "extract_remap.py"
RUNNER = APP_DIR / "pipeline_runner.py"

JOB_FILE = "job.json"
SPEC_FILE = "job_spec.json"
JOB_LOG = "job.log"
EXTRACT_DIR = "extract"
MANIFEST = "manifest.tsv"


# ------------------------------------------------------------------ commands

def extract_dir(output_dir: str, pathogen: str) -> str:
    return os.path.join(output_dir, EXTRACT_DIR, pathogen)


def build_map_cmd(map_script: str, refs_dir: str, samples_dir: str, tsv_path: str,
                  output_dir: str, log_path: str, workers: int, threads: int,
                  trim: bool = True, prefilter: bool = True, dedup: bool = True,
                  config_path: Optional[str] = None) -> List[str]:
    cmd = [sys.executable, map_script,
           "--refs-dir", refs_dir, "--samples-dir", samples_dir,
           "--output", tsv_path, "--output-dir", output_dir,
           "--log-file", log_path,
           "--workers", str(workers), "--threads", str(threads)]
    if config_path:
        cmd += ["--config", config_path]
    if not trim:
        cmd.append("--no-trim")
    if not prefilter:
        cmd.append("--no-prefilter")
    if not dedup:
        cmd.append("--no-dedup")
    return cmd


def build_extract_cmd(extract_script: str, ref_fasta: str, manifest: str,
                      outdir: str, threads: int) -> List[str]:
    return [sys.executable, extract_script,
            "--ref", ref_fasta, "--manifest", manifest,
            "--outdir", outdir, "--threads", str(threads)]


# ------------------------------------------------------------- sample registry

def mapped_bam(output_dir: str, pathogen: str, sample: str, assembly_id: str) -> str:
    return os.path.join(output_dir, pathogen, f"{sample}_{assembly_id}.bam")


def build_manifest(samples_dir: str, output_dir: str, pathogen: str, assembly_id: str,
                   out_dir: str) -> Tuple[str, List[Dict[str, str]]]:
    """Register the samples with a mapped BAM for this reference as an extract manifest.

    Writes <out_dir>/manifest.tsv (name, r1, r2, bam; '-' = no R2) and returns
    (manifest path, rows). Samples without a mapped BAM are not registered.
    Raises ValueError when no sample has a BAM.
    """
    rows = []
    for s in find_samples(samples_dir):
        bam = mapped_bam(output_dir, pathogen, s["name"], assembly_id)
        if not os.path.isfile(bam):
            continue
        rows.append({
            "name": s["name"],
            "r1": os.path.abspath(s["r1"]),
            "r2": os.path.abspath(s["r2"]) if s["r2"] else "-",
            "bam": os.path.abspath(bam),
        })
    if not rows:
        raise ValueError(f"no mapped BAM for {pathogen}/{assembly_id} in {output_dir}")
    os.makedirs(out_dir, exist_ok=True)
    manifest = os.path.join(out_dir, MANIFEST)
    with open(manifest, "w") as fh:
        fh.write("name\tr1\tr2\tbam\n")
        for r in rows:
            fh.write(f"{r['name']}\t{r['r1']}\t{r['r2']}\t{r['bam']}\n")
    return manifest, rows


def read_registered(out_dir: str) -> List[Dict[str, str]]:
    """Samples registered in an extract dir's manifest (empty if none)."""
    path = os.path.join(out_dir, MANIFEST)
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]


# ------------------------------------------------------------------- job state

def job_path(output_dir: str) -> str:
    return os.path.join(output_dir, JOB_FILE)


def read_job(output_dir: str) -> Optional[dict]:
    try:
        with open(job_path(output_dir)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def write_job(output_dir: str, data: dict) -> None:
    path = job_path(output_dir)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
    os.replace(tmp, path)
    run_id = data.get("run_id")
    if run_id:
        keys = ("status", "stage", "pid", "started", "finished", "pathogen",
                "assembly_id", "extract_dir", "samples", "manifest", "returncode",
                "error")
        fields = {key: data[key] for key in keys if key in data}
        runs.update_run(run_id, **fields)


def update_job(output_dir: str, **fields) -> None:
    data = read_job(output_dir) or {}
    data.update(fields)
    write_job(output_dir, data)


def pid_alive(pid: Optional[int]) -> bool:
    """True if pid exists and is not a zombie (a finished child not yet reaped)."""
    if not pid:
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


def job_active(output_dir: str) -> bool:
    job = read_job(output_dir)
    return bool(job) and job.get("status") == "running" and pid_alive(job.get("pid"))


def job_log_tail(output_dir: str, lines: int = 20) -> str:
    try:
        with open(os.path.join(output_dir, JOB_LOG)) as fh:
            return "".join(fh.readlines()[-lines:])
    except OSError:
        return ""


def cancel_job(output_dir: str) -> bool:
    """Stop the running job and its child processes. Returns False if none was running.

    The runner is started in its own session, so its process group holds the mapping or
    extract child too; the whole group gets SIGTERM. The status is set to 'cancelled'
    first, and the runner does not overwrite it.
    """
    output_dir = os.path.abspath(output_dir)
    if not job_active(output_dir):
        return False
    job = read_job(output_dir)
    update_job(output_dir, status="cancelled", stage=job.get("stage"),
               finished=time.strftime("%Y-%m-%d %H:%M:%S"))
    try:
        os.killpg(job["pid"], signal.SIGTERM)
    except ProcessLookupError:
        pass
    if job.get("run_id"):
        runs.release_lock(job["run_id"])
    return True


def start_job(spec: dict, output_dir: str, run_id: Optional[str] = None) -> int:
    """Start a detached pipeline job. Returns the runner pid.

    spec keys: map_cmd, extract_script, ref_fasta, samples_dir, pathogen,
    assembly_id, threads. Raises RuntimeError when a job is already running.
    """
    output_dir = os.path.abspath(output_dir)
    active = runs.active_run()
    if active:
        raise RuntimeError(f"a job is already running in run {active}")
    if job_active(output_dir):
        raise RuntimeError(f"a job is already running in {output_dir}")
    if run_id and runs.read_run(run_id) is None:
        raise RuntimeError(f"unknown run id: {run_id}")
    if run_id and not runs.acquire_lock(run_id):
        raise RuntimeError(f"could not acquire the run lock for {run_id}")
    try:
        os.makedirs(output_dir, exist_ok=True)
        spec = dict(spec, output_dir=output_dir, run_id=run_id,
                    extract_dir=extract_dir(output_dir, spec["pathogen"]))
        spec_path = os.path.join(output_dir, SPEC_FILE)
        with open(spec_path, "w") as fh:
            json.dump(spec, fh, indent=2)
        with open(os.path.join(output_dir, JOB_LOG), "w") as log:
            proc = subprocess.Popen(
                [sys.executable, str(RUNNER), spec_path],
                stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                start_new_session=True)
        write_job(output_dir, {
            "status": "running", "stage": "starting", "pid": proc.pid,
            "started": time.strftime("%Y-%m-%d %H:%M:%S"),
            "pathogen": spec["pathogen"], "assembly_id": spec["assembly_id"],
            "extract_dir": spec["extract_dir"], "run_id": run_id,
        })
        return proc.pid
    except Exception:
        if run_id:
            runs.release_lock(run_id)
        raise
