"""Detached runner for one mapping + extract job (started by jobs.start_job).

Usage: python pipeline_runner.py <job_spec.json>

Stages: mapping (map_references.py) -> per-reference extract (extract_remap.py).
Progress and final status are written to job.json in the output dir.
Child process output goes to job.log (inherited stdout/stderr).
"""
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import jobs  # noqa: E402
import runs  # noqa: E402
from map_references import find_references  # noqa: E402


def _now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _cancelled(out):
    job = jobs.read_job(out) or {}
    return job.get("status") == "cancelled"


def _result(pathogen, assembly_id, status, **extra):
    row = {
        "pathogen": pathogen, "assembly_id": assembly_id, "status": status,
        "returncode": None, "error": "", "extract_dir": "", "samples": [], "manifest": "",
    }
    row.update(extra)
    return row


def main(argv):
    with open(argv[1]) as fh:
        spec = json.load(fh)
    out = spec["output_dir"]
    run_id = spec.get("run_id")
    try:
        try:
            jobs.update_job(out, status="running", stage="mapping", extract_results=[])
            rc = subprocess.run(spec["map_cmd"]).returncode
            if rc != 0:
                if not _cancelled(out):
                    jobs.update_job(out, status="failed", stage="mapping", returncode=rc,
                                    finished=_now())
                return rc

            refs_dir = spec.get("refs_dir") or (runs.run_paths(run_id)["references"] if run_id else None)
            if not refs_dir:
                raise ValueError("refs_dir is required for extract")
            refs = find_references(refs_dir)
            results = []
            any_failed = False
            any_done = False
            for i, ref in enumerate(refs, 1):
                if _cancelled(out):
                    break
                pathogen, assembly_id = ref["pathogen"], ref["assembly_id"]
                ext_dir = jobs.extract_dir(out, pathogen, assembly_id)
                stage = f"extract ({i}/{len(refs)}): {pathogen}/{assembly_id}"
                jobs.update_job(out, stage=stage, pathogen=pathogen, assembly_id=assembly_id,
                                extract_dir=ext_dir)
                try:
                    manifest, rows = jobs.build_manifest(
                        spec["samples_dir"], out, pathogen, assembly_id, ext_dir)
                except Exception as e:
                    any_failed = True
                    msg = f"{type(e).__name__}: {e}"
                    print(f"[pipeline_runner] manifest failed for {pathogen}/{assembly_id}: {msg}", flush=True)
                    results.append(_result(pathogen, assembly_id, "failed", error=msg, extract_dir=ext_dir))
                    jobs.update_job(out, extract_results=results)
                    continue
                names = [r["name"] for r in rows]
                if not names:
                    msg = f"no mapped BAM for {pathogen}/{assembly_id}; skipped"
                    print(f"[pipeline_runner] {msg}", flush=True)
                    results.append(_result(pathogen, assembly_id, "skipped", error=msg,
                                           extract_dir=ext_dir, manifest=manifest))
                    jobs.update_job(out, extract_results=results)
                    continue

                print(f"[pipeline_runner] registered {len(names)} sample(s) for {pathogen}/{assembly_id}: {', '.join(names)}",
                      flush=True)
                jobs.update_job(out, stage=stage, samples=names, manifest=manifest, extract_results=results)
                cmd = jobs.build_extract_cmd(spec["extract_script"], ref["fasta_path"], manifest,
                                             ext_dir, spec["threads"])
                rc_ref = subprocess.run(cmd).returncode
                if rc_ref == 0 and not _cancelled(out):
                    removed = jobs.cleanup_full_ref_bam(out, pathogen, assembly_id)
                    if removed:
                        print(f"[pipeline_runner] cleaned {len(removed)} full-reference BAM file(s) for {pathogen}/{assembly_id}", flush=True)
                    any_done = True
                    results.append(_result(pathogen, assembly_id, "done", returncode=rc_ref,
                                           extract_dir=ext_dir, samples=names, manifest=manifest))
                elif _cancelled(out):
                    results.append(_result(pathogen, assembly_id, "cancelled", returncode=rc_ref,
                                           extract_dir=ext_dir, samples=names, manifest=manifest))
                else:
                    any_failed = True
                    results.append(_result(pathogen, assembly_id, "failed", returncode=rc_ref,
                                           extract_dir=ext_dir, samples=names, manifest=manifest,
                                           error=f"extract exited {rc_ref}"))
                jobs.update_job(out, extract_results=results)

            if not _cancelled(out):
                if not refs:
                    status, err = "failed", "no references found for extract"
                elif any_failed or not any_done:
                    status, err = "failed", ""
                else:
                    status, err = "done", ""
                jobs.update_job(out, status=status, stage="extract", error=err,
                                returncode=0 if status == "done" else 1,
                                finished=_now(), extract_results=results)
            return 0 if not any_failed and any_done else 1
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
