"""Read-only data access for the mapping_tracker app.

Two sources:
- TSV: mapping_statistics.tsv written by map_references.py after a run.
- live: rebuilt from completed sample x reference combinations in the output
  directory. A combination is complete when its BAM, BAI, .stats.txt and
  .aligner marker (with total_reads) exist and the .stats.txt was written
  after the BAM (summarize_bam runs last in map_sample).

Nothing in this module writes to the pipeline output directory.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

_PARENT = str(Path(__file__).resolve().parent.parent)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from map_references import (SAMTOOLS, fasta_total_length,  # noqa: E402
                            find_references, parse_flagstat, read_marker)

COLUMNS = [
    "sample", "pathogen", "assembly_id", "fasta_path", "total_reads",
    "mapped_reads", "unmapped_reads", "mapped_pct", "coverage", "avg_depth",
    "depth_covered", "ref_length", "duplicates_removed", "bam_path", "error",
]

_PROGRESS_RE = re.compile(r"\[(\d+)/(\d+)\]")
_PROCESSING_RE = re.compile(r"Processing (\d+) sample")


def _read_log(log_path):
    try:
        return Path(log_path).read_text().splitlines()
    except OSError:
        return None


def run_status(log_path):
    """'running', 'finished', 'failed' or 'unknown' for the last run in the log."""
    lines = _read_log(log_path)
    if lines is None:
        return "unknown"
    started = None
    for i, line in enumerate(lines):
        if "Starting mapping pipeline" in line:
            started = i
    if started is None:
        return "unknown"
    tail = lines[started:]
    if any("Pipeline completed" in line for line in tail):
        return "finished"
    if any("Fatal error" in line for line in tail):
        return "failed"
    return "running"


def read_progress(log_path):
    """(done, total) of the last run in the log, or None if unknown."""
    lines = _read_log(log_path)
    if lines is None:
        return None
    text = "\n".join(lines)
    totals = _PROCESSING_RE.findall(text)
    if not totals:
        return None
    total = int(totals[-1])
    done = 0
    for done_s, _ in _PROGRESS_RE.findall(text):
        done = int(done_s)
    return done, total


def parse_stats_sn(stats_path):
    """Summary numbers ('SN' lines) from a samtools stats file."""
    out = {}
    with open(stats_path) as fh:
        for line in fh:
            if not line.startswith("SN\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            key = parts[1].rstrip(":")
            value = parts[2]
            try:
                if any(c in value for c in ".eE"):
                    out[key] = float(value)
                else:
                    out[key] = int(value)
            except ValueError:
                out[key] = value
    return out


def _coverage_summary(text):
    """(coverage %, avg_depth, depth_covered) from `samtools coverage` output."""
    total_len = 0
    total_covbases = 0
    weighted_depth = 0.0
    for line in text.split("\n"):
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        try:
            contig_len = int(parts[2]) - int(parts[1]) + 1
            covbases = int(parts[4])
            meandepth = float(parts[6])
        except ValueError:
            continue
        total_len += contig_len
        total_covbases += covbases
        weighted_depth += meandepth * contig_len
    if total_len == 0:
        return None, None, None
    coverage = total_covbases / total_len * 100
    avg_depth = weighted_depth / total_len
    depth_covered = weighted_depth / total_covbases if total_covbases > 0 else 0.0
    return coverage, avg_depth, depth_covered


def combo_paths(output_dir, pathogen, sample, assembly_id):
    base = Path(output_dir) / pathogen / f"{sample}_{assembly_id}"
    return {
        "bam": str(base) + ".bam",
        "stats": str(base) + ".stats.txt",
        "aligner": str(base) + ".bam.aligner",
        "deduped": str(base) + ".bam.deduped",
    }


def _is_complete(paths):
    bam, stats = paths["bam"], paths["stats"]
    if not (os.path.exists(bam) and os.path.exists(bam + ".bai")
            and os.path.exists(stats)):
        return False
    if os.path.getmtime(stats) < os.path.getmtime(bam):
        return False  # BAM rewritten after the last stats: combo is being re-mapped
    marker = read_marker(paths["aligner"]) or {}
    try:
        int(marker["total_reads"])
    except (KeyError, ValueError):
        return False
    return True


def discover_combinations(output_dir, references):
    """Yield (sample, reference) for every .stats.txt found under output_dir."""
    for ref in references:
        pathogen_dir = Path(output_dir) / ref["pathogen"]
        if not pathogen_dir.is_dir():
            continue
        suffix = f"_{ref['assembly_id']}.stats.txt"
        for stats in sorted(pathogen_dir.glob("*" + suffix)):
            yield stats.name[:-len(suffix)], ref


def _summarize_combo(output_dir, ref, sample):
    paths = combo_paths(output_dir, ref["pathogen"], sample, ref["assembly_id"])
    bam = paths["bam"]
    marker = read_marker(paths["aligner"]) or {}
    total_reads = int(marker["total_reads"])
    row = {
        "sample": sample,
        "pathogen": ref["pathogen"],
        "assembly_id": ref["assembly_id"],
        "fasta_path": ref["fasta_path"],
        "total_reads": total_reads,
        "mapped_reads": 0,
        "unmapped_reads": 0,
        "mapped_pct": 0.0,
        "coverage": None,
        "avg_depth": None,
        "depth_covered": None,
        "ref_length": ref["ref_len"],
        "duplicates_removed": 0,
        "bam_path": bam,
        "error": "",
    }
    deduped = read_marker(paths["deduped"]) or {}
    try:
        row["duplicates_removed"] = int(deduped.get("duplicates_removed", 0))
    except ValueError:
        pass
    try:
        flagstat = subprocess.run([SAMTOOLS, "flagstat", bam], capture_output=True,
                                  text=True, check=True).stdout
        coverage_out = subprocess.run([SAMTOOLS, "coverage", bam],
                                      capture_output=True, text=True,
                                      check=True).stdout
    except (subprocess.CalledProcessError, OSError) as e:
        row["error"] = f"live summary failed: {e}"
        return row
    _, mapped = parse_flagstat(flagstat)
    row["mapped_reads"] = mapped
    row["unmapped_reads"] = total_reads - mapped
    row["mapped_pct"] = (mapped / total_reads * 100) if total_reads > 0 else 0
    row["coverage"], row["avg_depth"], row["depth_covered"] = _coverage_summary(
        coverage_out)
    return row


def load_tsv(path):
    df = pd.read_csv(path, sep="\t")
    for col in COLUMNS:
        if col not in df.columns:
            df[col] = None
    df["error"] = df["error"].fillna("").astype(str)
    df["duplicates_removed"] = df["duplicates_removed"].fillna(0)
    return df


def load_live(refs_dir, output_dir):
    """Build the statistics table from completed combinations only."""
    references = find_references(refs_dir)
    if not references:
        raise ValueError(f"No reference genomes found in {refs_dir}")
    for ref in references:
        ref["ref_len"] = fasta_total_length(ref["fasta_path"])
    rows = []
    for sample, ref in discover_combinations(output_dir, references):
        paths = combo_paths(output_dir, ref["pathogen"], sample, ref["assembly_id"])
        if _is_complete(paths):
            rows.append(_summarize_combo(output_dir, ref, sample))
    df = pd.DataFrame(rows, columns=COLUMNS)
    return df


def load_table(source, tsv_path, refs_dir, output_dir, log_path):
    """Return (dataframe, mode, run_status).

    source: 'auto', 'tsv' or 'live'. In auto mode the TSV is used when it exists
    and no run is in progress; otherwise the table is built live.
    """
    status = run_status(log_path)
    if source == "tsv":
        mode = "tsv"
    elif source == "live":
        mode = "live"
    else:
        mode = "live" if (status == "running" or not os.path.exists(tsv_path)) else "tsv"

    if mode == "tsv":
        df = load_tsv(tsv_path)
    else:
        df = load_live(refs_dir, output_dir)
    df["bam_exists"] = df["bam_path"].apply(
        lambda p: isinstance(p, str) and os.path.exists(p)).astype(bool)
    return df, mode, status


def signature(tsv_path, log_path, output_dir):
    """Cheap fingerprint of the inputs; changes whenever the data may change."""
    parts = []
    for p in (tsv_path, log_path):
        try:
            parts.append((str(p), os.stat(p).st_mtime_ns))
        except OSError:
            parts.append((str(p), None))
    for stats in sorted(Path(output_dir).glob("*/*.stats.txt")):
        parts.append((str(stats), stats.stat().st_mtime_ns))
    return tuple(parts)
