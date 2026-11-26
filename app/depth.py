"""Binned per-position depth from a BAM (opt-in, cached on disk).

Each contig is split into up to DEFAULT_BINS bins of equal width; the bin value
is the mean depth over all positions in the bin (positions with zero reads
included, via `samtools depth -a`). Results are cached under app/cache/.
"""
from __future__ import annotations

import math
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd

_PARENT = str(Path(__file__).resolve().parent.parent)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from map_references import SAMTOOLS  # noqa: E402

DEFAULT_BINS = 1000
CACHE_DIR = Path(__file__).resolve().parent / "cache"


def depth_artifact_path(bam_path, bins=DEFAULT_BINS):
    base = bam_path[:-len(".bam")] if str(bam_path).endswith(".bam") else str(bam_path)
    return f"{base}.depth.bins{bins}.tsv"


def contig_lengths(bam_path):
    """Ordered {contig: length} from the BAM header (@SQ lines)."""
    out = subprocess.run([SAMTOOLS, "view", "-H", bam_path],
                         capture_output=True, text=True, check=True).stdout
    lengths = {}
    for line in out.split("\n"):
        if not line.startswith("@SQ"):
            continue
        fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
        if "SN" in fields and "LN" in fields:
            lengths[fields["SN"]] = int(fields["LN"])
    return lengths


def _cache_path(bam_path, bins, cache_dir):
    st = os.stat(bam_path)
    name = f"{Path(bam_path).name}.{st.st_mtime_ns}.{st.st_size}.bins{bins}.tsv"
    return Path(cache_dir) / name


def _compute_binned_depth(bam_path, bins=DEFAULT_BINS):
    lengths = contig_lengths(bam_path)
    layout = {}
    sums = {}
    for contig, length in lengths.items():
        if length <= 0:
            continue
        width = max(1, math.ceil(length / bins))
        n_bins = math.ceil(length / width)
        layout[contig] = (length, width, n_bins)
        sums[contig] = [0.0] * n_bins

    proc = subprocess.Popen([SAMTOOLS, "depth", "-a", bam_path],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True)
    for line in proc.stdout:
        contig, pos, depth = line.rstrip("\n").split("\t")[:3]
        if contig not in layout:
            continue
        _, width, _ = layout[contig]
        sums[contig][(int(pos) - 1) // width] += float(depth)
    stderr = proc.stderr.read()
    if proc.wait() != 0:
        raise RuntimeError(f"samtools depth failed: {stderr[:500]}")

    rows = []
    for contig, (length, width, n_bins) in layout.items():
        for i in range(n_bins):
            start = i * width + 1
            end = min((i + 1) * width, length)
            rows.append({"contig": contig, "start": start, "end": end,
                         "mean_depth": sums[contig][i] / (end - start + 1)})
    return pd.DataFrame(rows, columns=["contig", "start", "end", "mean_depth"])


def ensure_depth_artifact(bam_path, bins=DEFAULT_BINS):
    """Compute and persist binned depth next to BAM; return artifact path."""
    artifact = depth_artifact_path(bam_path, bins=bins)
    if os.path.exists(artifact):
        return artifact
    if not os.path.exists(bam_path):
        return artifact
    df = _compute_binned_depth(bam_path, bins=bins)
    Path(artifact).parent.mkdir(parents=True, exist_ok=True)
    tmp = artifact + ".tmp"
    df.to_csv(tmp, sep="\t", index=False)
    os.replace(tmp, artifact)
    return artifact


def load_depth_artifact(bam_path, bins=DEFAULT_BINS):
    artifact = depth_artifact_path(bam_path, bins=bins)
    if os.path.exists(artifact):
        return pd.read_csv(artifact, sep="\t")
    return None


def binned_depth(bam_path, bins=DEFAULT_BINS, cache_dir=CACHE_DIR):
    """DataFrame(contig, start, end, mean_depth), 1-based inclusive bin coords."""
    artifact_df = load_depth_artifact(bam_path, bins=bins)
    if artifact_df is not None:
        return artifact_df
    cache_file = _cache_path(bam_path, bins, cache_dir)
    if cache_file.exists():
        return pd.read_csv(cache_file, sep="\t")
    if not os.path.exists(bam_path):
        raise FileNotFoundError(f"BAM and depth artifact not found for {bam_path}")

    df = _compute_binned_depth(bam_path, bins=bins)
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    df.to_csv(cache_file, sep="\t", index=False)
    return df
