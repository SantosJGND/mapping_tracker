#!/usr/bin/env python3

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

import mapping_config

REFS_DIR = "references"
SAMPLES_DIR = "."
OUTPUT_FILE = "mapping_statistics.tsv"
LOG_FILE = "mapping_log.txt"


def resolve_tool(env_name, default):
    value = os.environ.get(env_name) or default
    return shutil.which(value) or value


BWA = resolve_tool("MAPPING_BWA", "bwa")
MINIMAP2 = resolve_tool("MAPPING_MINIMAP2", "minimap2")
BWA_ID = os.path.basename(shutil.which(BWA) or BWA)
MINIMAP2_ID = os.path.basename(shutil.which(MINIMAP2) or MINIMAP2)
SAMTOOLS = resolve_tool("MAPPING_SAMTOOLS", "samtools")
MSAMTOOLS = resolve_tool("MAPPING_MSAMTOOLS", "msamtools")
PRINSEQ = resolve_tool("MAPPING_PRINSEQ", "prinseq++")
PICARD = resolve_tool("MAPPING_PICARD", "picard")
TRIMMOMATIC = resolve_tool("MAPPING_TRIMMOMATIC", "trimmomatic")

MIN_LEN = 100
MIN_IDENTITY = 0.98
MIN_QUERY_ALIGNED = 0.50

MIN_ENTROPY = 0.5
MIN_DUST = 0.5

TRIM_PHRED = 33
TRIM_LEADING = 3
TRIM_TRAILING = 3
TRIM_MINLEN = 36
TRIM_SLIDING_WINDOW = "5:20"

NANOFILT = resolve_tool("MAPPING_NANOFILT", "NanoFilt")
NANOFILT_QUALITY = 7
NANOFILT_MINLEN = 200
NANOFILT_HEADCROP = 0
NANOFILT_TAILCROP = 0


def resolve_aligner(value):
    """Return (aligner_id, aligner_bin) for an aligner name or full path."""
    if value == BWA_ID:
        return BWA_ID, BWA
    if value == MINIMAP2_ID:
        return MINIMAP2_ID, MINIMAP2
    if os.path.exists(value) and os.access(value, os.X_OK):
        base = os.path.basename(value)
        if "minimap2" in base:
            return MINIMAP2_ID, value
        if "bwa" in base:
            return BWA_ID, value
        raise ValueError(
            f"Unrecognized aligner binary (basename has neither 'bwa' nor "
            f"'minimap2'): {value}")
    raise ValueError(f"Aligner not found or not executable: {value}")


def setup_logging(log_file, verbose=False):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file, mode='w'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


def find_references(refs_dir):
    """Find reference genomes, auto-detecting nested NCBI vs flat layouts."""
    refs_path = Path(refs_dir)
    if not refs_path.is_dir():
        raise ValueError(f"References directory not found: {refs_dir}")

    references = []

    nested = any(
        (p / "ncbi_dataset" / "data").is_dir()
        for p in refs_path.iterdir()
        if p.is_dir()
    )

    if nested:
        for pathogen_dir in sorted(refs_path.iterdir()):
            if not pathogen_dir.is_dir():
                continue
            data_dir = pathogen_dir / "ncbi_dataset" / "data"
            if not data_dir.is_dir():
                continue

            best_by_accession = {}
            for assembly_dir in sorted(data_dir.iterdir()):
                if not assembly_dir.is_dir():
                    continue
                assembly_id = assembly_dir.name
                if not (assembly_id.startswith("GCF_") or assembly_id.startswith("GCA_")):
                    continue
                fasta_files = sorted(assembly_dir.glob("*.fna"))
                if not fasta_files:
                    continue
                accession = assembly_id.split("_", 1)[1]
                current = best_by_accession.get(accession)
                if current is None or (
                        assembly_id.startswith("GCF_") and current[0].startswith("GCA_")):
                    best_by_accession[accession] = (assembly_id, str(fasta_files[0]))

            for assembly_id, fasta_path in best_by_accession.values():
                references.append({
                    "pathogen": pathogen_dir.name,
                    "assembly_id": assembly_id,
                    "fasta_path": fasta_path
                })
    else:
        for fasta_file in sorted(list(refs_path.glob("*.fa"))
                                 + list(refs_path.glob("*.fasta"))):
            stem = fasta_file.stem
            if stem.endswith(".ref"):
                stem = stem[:-4]
            parts = stem.split("_")
            pathogen = parts[0] if parts else stem
            references.append({
                "pathogen": pathogen,
                "assembly_id": stem,
                "fasta_path": str(fasta_file)
            })

    return references


def _sample_name(fastq):
    """Sample name for an unpaired FASTQ: drop the whole .fastq.gz suffix."""
    name = fastq.name
    return name[:-len(".fastq.gz")] if name.endswith(".fastq.gz") else fastq.stem


def find_samples(samples_dir, nanopore=False):
    """Find samples, auto-detecting nested (per-subdir) vs flat layouts.

    In nanopore mode reads are single-end long reads: every FASTQ becomes its
    own sample and no read-pairing is attempted.
    """
    base_path = Path(samples_dir)
    if not base_path.is_dir():
        raise ValueError(f"Samples directory not found: {samples_dir}")

    samples = []
    fastqs = sorted(base_path.glob("*.fastq.gz"))
    fastqs = [f for f in fastqs if not f.name.startswith(".")]

    def _se(fastq):
        return {"name": _sample_name(fastq), "r1": str(fastq), "r2": None,
                "is_paired": False}

    if fastqs:
        if nanopore:
            samples = [_se(f) for f in fastqs]
            return samples
        used = set()
        pair_markers = [("_R1", "_R2"), ("_R1_001", "_R2_001"), ("_1P", "_2P"),
                        ("_1", "_2"), (".R1", ".R2")]
        for fastq in fastqs:
            if fastq.name in used:
                continue
            r2_name = None
            fwd_marker = None
            for fwd, rev in pair_markers:
                idx = fastq.name.find(fwd)
                if idx == -1:
                    continue
                after = fastq.name[idx + len(fwd):]
                if after and after[0].isdigit():
                    continue
                candidate = fastq.name.replace(fwd, rev, 1)
                if (fastq.parent / candidate).exists():
                    r2_name = candidate
                    fwd_marker = fwd
                    break
            r2 = None
            if r2_name:
                r2_path = fastq.parent / r2_name
                if r2_path.exists():
                    r2 = str(r2_path)
                    used.add(r2_path.name)
            if r2 and fwd_marker:
                idx = fastq.stem.find(fwd_marker)
                name = fastq.stem[:idx] if idx != -1 else fastq.stem
            else:
                name = _sample_name(fastq)
            samples.append({
                "name": name,
                "r1": str(fastq),
                "r2": r2,
                "is_paired": r2 is not None
            })
            used.add(fastq.name)
    else:
        for sample_dir in sorted(base_path.iterdir()):
            if not sample_dir.is_dir():
                continue
            if sample_dir.name.startswith(".") or sample_dir.name == "references":
                continue
            all_files = sorted(sample_dir.glob("*.fastq.gz"))
            if nanopore:
                samples.extend(_se(f) for f in all_files)
                continue
            r1_files = [f for f in all_files if "_R1" in f.name]
            if r1_files:
                sample_name = sample_dir.name
                r1 = str(r1_files[0])
                r2_pattern = r1.replace("_R1", "_R2")
                r2 = r2_pattern if os.path.exists(r2_pattern) else None
                samples.append({
                    "name": sample_name,
                    "r1": r1,
                    "r2": r2,
                    "is_paired": r2 is not None
                })
            else:
                se_files = sorted(sample_dir.glob("*.fastq.gz"))
                if se_files:
                    samples.append({
                        "name": sample_dir.name,
                        "r1": str(se_files[0]),
                        "r2": None,
                        "is_paired": False
                    })

    return samples


def tool_available(path):
    if not path:
        return False
    if os.path.isabs(path):
        return os.path.exists(path) and os.access(path, os.X_OK)
    return shutil.which(path) is not None


def check_dependencies(msamtools_path=None, need_msamtools=False,
                       prinseq_path=None, need_prinseq=False,
                       picard_path=None, need_picard=False, aligner=BWA_ID,
                       aligner_bin=None, trimmomatic_path=None,
                       need_trimmomatic=False, nanofilt_path=None,
                       need_nanofilt=False):
    """Verify required software is available."""
    missing = []
    if aligner_bin and os.path.isabs(aligner_bin):
        if not tool_available(aligner_bin):
            missing.append(f"aligner ({aligner_bin})")
    elif not tool_available(aligner):
        missing.append(aligner)
    if not tool_available(SAMTOOLS):
        missing.append(SAMTOOLS)
    if need_trimmomatic and not tool_available(trimmomatic_path):
        missing.append(f"trimmomatic ({trimmomatic_path})")
    if need_msamtools and not tool_available(msamtools_path):
        missing.append(f"msamtools ({msamtools_path})")
    if need_prinseq and not tool_available(prinseq_path):
        missing.append(f"prinseq++ ({prinseq_path})")
    if need_picard and not tool_available(picard_path):
        missing.append(f"picard ({picard_path})")
    if need_nanofilt and not tool_available(nanofilt_path):
        missing.append(f"NanoFilt ({nanofilt_path})")
    if missing:
        raise RuntimeError(f"Required tools not found: {', '.join(missing)}")


def fasta_total_length(fasta_path):
    """Sum of all sequence lengths in a FASTA file (streaming)."""
    total = 0
    current = 0
    with open(fasta_path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                total += current
                current = 0
            else:
                current += len(line.strip())
        total += current
    return total


def index_path(temp_dir, assembly_id):
    """BWA index prefix (suffix files are added by `bwa index`) in temp_dir."""
    return os.path.join(temp_dir, assembly_id)


def build_index(ref_info, temp_dir, aligner=BWA_ID, aligner_bin=None):
    """Prepare whichever reference artifacts the aligner needs.

    minimap2 aligns against the reference FASTA directly (no index step);
    BWA builds a `.bwt`-family index into temp_dir. Returns the string to use
    as the reference when mapping.
    """
    fasta = ref_info["fasta_path"]
    bin_path = aligner_bin or aligner

    if aligner == MINIMAP2_ID:
        return fasta

    index_base = index_path(temp_dir, ref_info["assembly_id"])
    if not os.path.exists(index_base + ".bwt"):
        try:
            subprocess.run(
                [bin_path, "index", "-p", index_base, fasta],
                check=True, capture_output=True, text=True
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"BWA index failed for {ref_info['pathogen']}/{ref_info['assembly_id']}: "
                f"stderr: {e.stderr}"
            ) from e

    return index_base


def parse_flagstat(stdout):
    """Return (primary_total, primary_mapped) from samtools flagstat output.

    Only primary alignments are counted: secondary (0x100) and supplementary
    (0x800) records are excluded, so the total is a stable per-sample read
    count that does not depend on the reference or mapping ambiguity.
    """
    total_reads = 0
    mapped_reads = 0
    for line in stdout.split("\n"):
        if line.strip().endswith("primary"):
            total_reads = int(line.split()[0])
        elif "primary mapped (" in line:
            mapped_reads = int(line.split()[0])
    return total_reads, mapped_reads


def summarize_bam(bam_path, stats_path, sample_name, pathogen, assembly_id,
                  fasta_path, ref_len, threads, skipped, total_override=None):
    """Run flagstat/stats/coverage on a BAM and return the full statistics row.

    `total_override` sets total_reads independently of the BAM. It is the
    sample's primary read count captured before dedup/filter, so it is a fixed
    library size: dedup and filter physically drop reads and flagstat would
    otherwise report a reference-dependent total.
    """
    flagstat_result = subprocess.run(
        [SAMTOOLS, "flagstat", "-@", str(threads), bam_path],
        capture_output=True, text=True, check=True)
    stats_result = subprocess.run(
        [SAMTOOLS, "stats", "-@", str(threads), bam_path],
        capture_output=True, text=True, check=True)
    coverage_result = subprocess.run(
        [SAMTOOLS, "coverage", bam_path],
        capture_output=True, text=True, check=True)

    with open(stats_path, "w") as f:
        f.write(stats_result.stdout)

    flag_total, mapped_reads = parse_flagstat(flagstat_result.stdout)
    total_reads = total_override if total_override is not None else flag_total

    coverage = None
    avg_depth = None
    depth_covered = None
    total_len = 0
    total_covbases = 0
    weighted_depth = 0.0
    for line in coverage_result.stdout.split("\n"):
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
    if total_len > 0:
        coverage = (total_covbases / total_len) * 100
        avg_depth = weighted_depth / total_len
        depth_covered = (weighted_depth / total_covbases) if total_covbases > 0 else 0.0

    mapped_pct = (mapped_reads / total_reads * 100) if total_reads > 0 else 0

    return {
        "sample": sample_name,
        "pathogen": pathogen,
        "assembly_id": assembly_id,
        "fasta_path": fasta_path,
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "unmapped_reads": total_reads - mapped_reads,
        "mapped_pct": mapped_pct,
        "coverage": coverage,
        "avg_depth": avg_depth,
        "depth_covered": depth_covered,
        "ref_length": ref_len,
        "bam_path": bam_path,
        "skipped": skipped
    }


def filter_bam(bam_path, filter_config, total_reads):
    """Filter a BAM in place with msamtools (length/identity/query-aligned)."""
    tmp_path = bam_path + ".filtered.tmp"
    cmd = [
        filter_config["msamtools"], "filter", "-b",
        "-l", str(filter_config["min_len"]),
        "-p", str(filter_config["min_identity_pct"]),
        "-z", str(filter_config["min_query_pct"]),
        bam_path
    ]
    try:
        with open(tmp_path, "wb") as out:
            result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE,
                                    text=True)
        if result.returncode != 0:
            raise RuntimeError(f"msamtools filter (rc={result.returncode}): "
                               f"{result.stderr[:500]}")

        subprocess.run([SAMTOOLS, "index", tmp_path],
                       check=True, capture_output=True)

        os.remove(bam_path)
        os.remove(bam_path + ".bai")
        os.replace(tmp_path, bam_path)
        os.replace(tmp_path + ".bai", bam_path + ".bai")

        marker = bam_path + ".filtered"
        with open(marker, "w") as f:
            f.write(f"min_len={filter_config['min_len']}\n"
                    f"min_identity={filter_config['min_identity']}\n"
                    f"min_query_aligned={filter_config['min_query_aligned']}\n"
                    f"total_reads={total_reads}\n")
    except Exception:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        raise


def dedup_bam(bam_path, dedup_config, temp_dir):
    """Remove duplicate reads in place with Picard MarkDuplicates.

    Returns the number of duplicate reads removed (pre - post dedup flagstat
    totals). Writes `<bam>.dedup.metrics.txt` and a `<bam>.deduped` marker.
    """
    metrics_path = bam_path + ".dedup.metrics.txt"
    tmp_path = bam_path + ".dedup.tmp"
    picard_tmp = os.path.join(temp_dir, f"picard_{os.path.basename(bam_path)}")
    os.makedirs(picard_tmp, exist_ok=True)

    pre_total, _ = parse_flagstat(
        subprocess.run([SAMTOOLS, "flagstat", bam_path],
                       capture_output=True, text=True, check=True).stdout)

    cmd = [dedup_config["picard"], f"-Xmx{dedup_config['heap']}",
           "MarkDuplicates",
           f"INPUT={bam_path}",
           f"OUTPUT={tmp_path}",
           f"METRICS_FILE={metrics_path}",
           "REMOVE_DUPLICATES=true",
           "ASSUME_SORT_ORDER=coordinate",
           "CREATE_INDEX=false",
           "VALIDATION_STRINGENCY=SILENT",
           f"TMP_DIR={picard_tmp}"]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"picard MarkDuplicates (rc={result.returncode}): "
                               f"{result.stderr[:500]}")

        subprocess.run([SAMTOOLS, "index", tmp_path],
                       check=True, capture_output=True)

        os.remove(bam_path)
        os.remove(bam_path + ".bai")
        os.replace(tmp_path, bam_path)
        os.replace(tmp_path + ".bai", bam_path + ".bai")

        post_total, _ = parse_flagstat(
            subprocess.run([SAMTOOLS, "flagstat", bam_path],
                           capture_output=True, text=True, check=True).stdout)
        removed = max(0, pre_total - post_total)

        marker = bam_path + ".deduped"
        with open(marker, "w") as f:
            f.write(f"remove_duplicates=true\n"
                    f"pre_total={pre_total}\n"
                    f"post_total={post_total}\n"
                    f"duplicates_removed={removed}\n")
        return removed
    except Exception:
        for p in (tmp_path, tmp_path + ".bai"):
            if os.path.exists(p):
                os.remove(p)
        raise


def count_fastq_reads(path):
    """Count reads in a FASTQ file (gzip-aware).

    Counts a read only at the record boundary (line 0 of each 4-line record).
    Only lines starting with '@' at that position are headers; quality lines
    may also start with '@' (Phred+33 = 31) and must not be counted.
    """
    if not path or not os.path.exists(path):
        return 0
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0 and line.startswith("@"):
                n += 1
    return n


def read_marker(marker_path):
    """Read a key=value marker file into a dict (None if unreadable)."""
    d = {}
    try:
        with open(marker_path) as fh:
            for line in fh:
                if "=" in line:
                    key, value = line.rstrip("\n").split("=", 1)
                    d[key] = value
    except OSError:
        return None
    return d


def _nonempty(*paths):
    """Existing, non-empty paths as strings (None/empty entries dropped)."""
    return [str(p) for p in paths
            if p is not None and os.path.exists(p) and os.path.getsize(p) > 0]


def trim_sample(sample, config, threads):
    """Trim a sample's reads with Trimmomatic, caching results.

    Cached outputs live in `<output-dir>/trimmed/` and are reused when the
    stored parameters match the requested ones. PE orphan survivors are kept
    as single-end reads in `sample["singles"]`. The sample dict is updated in
    place and returned.
    """
    sample_name = sample["name"]
    out_dir = Path(config["dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    r1, r2 = sample["r1"], sample["r2"]
    out_r1 = out_dir / f"{sample_name}_R1.fq"
    out_r2 = out_dir / f"{sample_name}_R2.fq" if r2 else None
    unpaired_r1 = out_dir / f"{sample_name}_R1.single.fq" if r2 else None
    unpaired_r2 = out_dir / f"{sample_name}_R2.single.fq" if r2 else None
    marker_path = out_dir / f"{sample_name}.trimmed"

    params = {
        "phred": str(config["phred"]),
        "leading": str(config["leading"]),
        "trailing": str(config["trailing"]),
        "minlen": str(config["minlen"]),
        "sliding_window": config["sliding_window"],
        "source_r1": r1,
        "source_r2": r2 or "",
    }

    expected_outputs = [out_r1]
    if out_r2:
        expected_outputs.append(out_r2)

    cached = (read_marker(marker_path) == params
              and all(os.path.exists(p) for p in expected_outputs))
    if cached:
        sample["r1"] = str(out_r1)
        sample["r2"] = str(out_r2) if out_r2 else None
        sample["singles"] = _nonempty(unpaired_r1, unpaired_r2)
        sample["trimmed"] = True
        sample["trim_fingerprint"] = "|".join(
            f"{k}={v}" for k, v in sorted(params.items()))
        return sample

    for path in [out_r1, out_r2, unpaired_r1, unpaired_r2]:
        if path is not None and path.exists():
            path.unlink()

    trimmers = [f"LEADING:{config['leading']}",
                f"TRAILING:{config['trailing']}",
                f"SLIDINGWINDOW:{config['sliding_window']}",
                f"MINLEN:{config['minlen']}"]

    if r2:
        cmd = [config["trimmomatic"], "PE", "-threads", str(threads),
               f"-phred{config['phred']}", r1, r2,
               str(out_r1), str(unpaired_r1), str(out_r2), str(unpaired_r2)]
    else:
        cmd = [config["trimmomatic"], "SE", "-threads", str(threads),
               f"-phred{config['phred']}", r1, str(out_r1)]
    cmd += trimmers

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"Trimmomatic failed for {sample_name} (rc={result.returncode}): "
            f"{result.stderr[:500]}"
        )

    with open(marker_path, "w") as f:
        for key, value in params.items():
            f.write(f"{key}={value}\n")

    input_reads = (count_fastq_reads(r1)
                   + (count_fastq_reads(r2) if r2 else 0))
    singles = _nonempty(unpaired_r1, unpaired_r2)
    passed_reads = (count_fastq_reads(str(out_r1))
                    + (count_fastq_reads(str(out_r2)) if out_r2 else 0)
                    + sum(count_fastq_reads(p) for p in singles))
    sample["r1"] = str(out_r1)
    sample["r2"] = str(out_r2) if out_r2 else None
    sample["singles"] = singles
    sample["trimmed"] = True
    sample["trim_fingerprint"] = "|".join(
        f"{k}={v}" for k, v in sorted(params.items()))
    sample["trim_log"] = {
        "input_reads": input_reads,
        "passed_reads": passed_reads,
        "removed_reads": input_reads - passed_reads,
    }
    return sample


def prefilter_sample(sample, config, threads):
    """Filter a sample's reads with prinseq++ (entropy/dust), caching results.

    Cached outputs live in `<output-dir>/prefiltered/` and are reused when the
    stored parameters match the requested ones. Single-end reads carried over
    from trimming (`sample["singles"]`) are filtered as single-end too, and
    the trim fingerprint is part of the cache key so that re-trimming
    invalidates this cache. The sample dict is updated in place with the
    filtered paths and returned.
    """
    sample_name = sample["name"]
    out_dir = Path(config["dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    ext = ".fq.gz" if config["gz"] else ".fq"
    good_r1 = out_dir / f"{sample_name}_R1{ext}"
    good_r2 = out_dir / f"{sample_name}_R2{ext}" if sample["r2"] else None
    single_r1 = out_dir / f"{sample_name}_R1.single{ext}"
    single_r2 = out_dir / f"{sample_name}_R2.single{ext}" if sample["r2"] else None
    marker_path = out_dir / f"{sample_name}.prefiltered"

    input_singles = list(sample.get("singles", []))
    single_outputs = [out_dir / f"{sample_name}_s{i}{ext}"
                      for i in range(len(input_singles))]

    params = {
        "min_entropy": str(config["min_entropy"]),
        "min_dust": str(config["min_dust"]),
        "gz": str(config["gz"]),
        "source_r1": sample["r1"],
        "source_r2": sample["r2"] or "",
        "source_singles": ";".join(input_singles),
        "source_trimmed": sample.get("trim_fingerprint", ""),
        "source_nanofiltered": sample.get("nanofilt_fingerprint", ""),
    }

    expected_outputs = [good_r1] + single_outputs
    if good_r2:
        expected_outputs.append(good_r2)

    cached = (read_marker(marker_path) == params
              and all(os.path.exists(p) for p in expected_outputs))
    if cached:
        sample["r1"] = str(good_r1)
        sample["r2"] = str(good_r2) if good_r2 else None
        sample["singles"] = _nonempty(single_r1, single_r2, *single_outputs)
        sample["prefiltered"] = True
        return sample

    for path in ([good_r1, good_r2, single_r1, single_r2]
                 + single_outputs):
        if path is not None and path.exists():
            path.unlink()

    cmd = [config["prinseq"], "-threads", str(threads),
           "-fastq", sample["r1"],
           "-lc_entropy", str(config["min_entropy"]),
           "-lc_dust", str(config["min_dust"])]
    if sample["r2"]:
        cmd += ["-fastq2", sample["r2"],
                "-out_good", str(good_r1),
                "-out_good2", str(good_r2),
                "-out_single", str(single_r1),
                "-out_single2", str(single_r2),
                "-out_bad", "/dev/null",
                "-out_bad2", "/dev/null"]
    else:
        cmd += ["-out_good", str(good_r1),
                "-out_single", str(single_r1),
                "-out_bad", "/dev/null"]
    if config["gz"]:
        cmd.append("-out_gz")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"prinseq++ prefilter failed for {sample_name} (rc={result.returncode}): "
            f"{result.stderr[:500]}"
        )

    for src, dst in zip(input_singles, single_outputs):
        se_cmd = [config["prinseq"], "-threads", str(threads),
                  "-fastq", src,
                  "-lc_entropy", str(config["min_entropy"]),
                  "-lc_dust", str(config["min_dust"]),
                  "-out_good", str(dst),
                  "-out_bad", "/dev/null"]
        if config["gz"]:
            se_cmd.append("-out_gz")
        result = subprocess.run(se_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"prinseq++ prefilter failed for {sample_name} single-end reads "
                f"(rc={result.returncode}): {result.stderr[:500]}"
            )

    with open(marker_path, "w") as f:
        for key in params:
            f.write(f"{key}={params[key]}\n")

    input_reads = (count_fastq_reads(sample["r1"])
                   + (count_fastq_reads(sample["r2"]) if sample["r2"] else 0)
                   + sum(count_fastq_reads(p) for p in input_singles))
    passed_reads = (count_fastq_reads(good_r1)
                    + (count_fastq_reads(good_r2) if good_r2 else 0)
                    + count_fastq_reads(single_r1)
                    + (count_fastq_reads(single_r2) if single_r2 else 0)
                    + sum(count_fastq_reads(str(p)) for p in single_outputs))
    sample["r1"] = str(good_r1)
    sample["r2"] = str(good_r2) if good_r2 else None
    sample["singles"] = _nonempty(single_r1, single_r2, *single_outputs)
    sample["prefiltered"] = True
    sample["prefilter_log"] = {
        "input_reads": input_reads,
        "passed_reads": passed_reads,
        "removed_reads": input_reads - passed_reads,
    }
    return sample


def nanofilt_sample(sample, config, threads):
    """Filter a nanopore sample's reads with NanoFilt, caching results.

    Reads are single-end long reads. NanoFilt requires decompressed FASTQ on
    stdin (its documented `gunzip -c ... | NanoFilt` interface), so gzip inputs
    are piped through `gzip -dc` (no temporary files). Filtered reads are
    written to `<output-dir>/nanofiltered/<sample>.fq`. A
    `<sample>.nanofiltered` marker records the parameters used, so a re-run with
    unchanged settings reuses the cached reads. The sample dict is updated in
    place and returned.
    """
    sample_name = sample["name"]
    out_dir = Path(config["dir"])
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / f"{sample_name}.fq"
    marker_path = out_dir / f"{sample_name}.nanofiltered"

    input_paths = [sample["r1"]]
    if sample["r2"]:
        input_paths.append(sample["r2"])
    input_paths += list(sample.get("singles", []))

    params = {
        "quality": str(config["quality"]),
        "minlen": str(config["minlen"]),
        "headcrop": str(config["headcrop"]),
        "tailcrop": str(config["tailcrop"]),
        "source": ";".join(input_paths),
    }

    sample["nanofilt_fingerprint"] = "|".join(
        f"{k}={v}" for k, v in sorted(params.items()))

    cached = (read_marker(marker_path) == params
              and out_path.exists() and os.path.getsize(out_path) > 0)
    if cached:
        sample["r1"] = str(out_path)
        sample["r2"] = None
        sample["singles"] = []
        sample["nanofiltered"] = True
        return sample

    if out_path.exists():
        out_path.unlink()

    input_reads = 0
    passed_reads = 0
    try:
        for i, path in enumerate(input_paths):
            cmd = [config["nanofilt"], "-q", str(config["quality"]),
                   "-l", str(config["minlen"]),
                   "--logfile", str(out_path) + ".log"]
            if config["headcrop"]:
                cmd += ["--headcrop", str(config["headcrop"])]
            if config["tailcrop"]:
                cmd += ["--tailcrop", str(config["tailcrop"])]
            src_proc = None
            if str(path).endswith(".gz"):
                src_proc = subprocess.Popen(
                    ["gzip", "-dc", path], stdout=subprocess.PIPE)
                stdin = src_proc.stdout
            else:
                stdin = open(path)
            try:
                with open(out_path, "w" if i == 0 else "a") as dst:
                    result = subprocess.run(cmd, stdin=stdin, stdout=dst,
                                            stderr=subprocess.PIPE, text=True)
            finally:
                if src_proc is not None:
                    src_proc.stdout.close()
                    src_proc.wait()
                else:
                    stdin.close()
            if result.returncode != 0:
                raise RuntimeError(
                    f"NanoFilt failed for {sample_name} (rc={result.returncode}): "
                    f"{result.stderr[:500]}")
            input_reads += count_fastq_reads(path)
            passed_reads = count_fastq_reads(str(out_path))
    except Exception:
        if out_path.exists():
            out_path.unlink()
        raise

    with open(marker_path, "w") as f:
        for key, value in params.items():
            f.write(f"{key}={value}\n")

    sample["r1"] = str(out_path)
    sample["r2"] = None
    sample["singles"] = []
    sample["nanofiltered"] = True
    sample["nanofilt_log"] = {
        "input_reads": input_reads,
        "passed_reads": passed_reads,
        "removed_reads": input_reads - passed_reads,
    }
    return sample


def map_sample(sample, ref_info, temp_dir, output_dir, threads, filter_config,
               prefilter_config, dedup_config, aligner=BWA_ID, trim_state="",
               aligner_bin=None, nanopore=False, sort_memory=None):
    sample_name = sample["name"]
    pathogen = ref_info["pathogen"]
    assembly_id = ref_info["assembly_id"]

    output_subdir = os.path.join(output_dir, pathogen)
    os.makedirs(output_subdir, exist_ok=True)

    bam_path = os.path.join(output_subdir, f"{sample_name}_{assembly_id}.bam")
    stats_path = os.path.join(output_subdir, f"{sample_name}_{assembly_id}.stats.txt")

    ref_len = ref_info.get("ref_len") or fasta_total_length(ref_info["fasta_path"])

    ready = (os.path.exists(bam_path)
             and os.path.exists(bam_path + ".bai")
             and os.path.exists(stats_path))
    if filter_config["enabled"]:
        ready = ready and os.path.exists(bam_path + ".filtered")
    if prefilter_config["enabled"]:
        ready = ready and os.path.exists(bam_path + ".prefiltered")
    if dedup_config["enabled"]:
        ready = ready and os.path.exists(bam_path + ".deduped")
    recorded = read_marker(bam_path + ".aligner") or {}
    if recorded.get("aligner", BWA_ID) != aligner:
        ready = False
    if recorded.get("trimmed", "") != trim_state:
        ready = False
    requested_platform = "nanopore" if nanopore else "illumina"
    if recorded.get("platform", "illumina") != requested_platform:
        ready = False
    if "total_reads" not in recorded:
        ready = False

    if ready:
        try:
            total_override = None
            try:
                total_override = int(recorded["total_reads"])
            except (KeyError, ValueError):
                if filter_config["enabled"]:
                    marker = bam_path + ".filtered"
                    try:
                        with open(marker) as f:
                            for line in f:
                                if line.startswith("total_reads="):
                                    total_override = int(line.strip().split("=", 1)[1])
                    except (OSError, ValueError):
                        total_override = None
            result = summarize_bam(bam_path, stats_path, sample_name, pathogen,
                                   assembly_id, ref_info["fasta_path"], ref_len,
                                   threads, skipped=True,
                                   total_override=total_override)
            if dedup_config["enabled"]:
                marker = read_marker(bam_path + ".deduped")
                result["duplicates_removed"] = (
                    int(marker.get("duplicates_removed")) if marker else 0)
            return result
        except subprocess.CalledProcessError as e:
            return {
                "sample": sample_name,
                "pathogen": pathogen,
                "assembly_id": assembly_id,
                "fasta_path": ref_info["fasta_path"],
                "error": f"Resume stats failed (rc={e.returncode}): "
                         f"{e.stderr[:500]}",
                "skipped": True
            }
        except Exception as e:
            return {
                "sample": sample_name,
                "pathogen": pathogen,
                "assembly_id": assembly_id,
                "fasta_path": ref_info["fasta_path"],
                "error": f"Unexpected error: {type(e).__name__}: {str(e)}",
                "skipped": True
            }

    if aligner == MINIMAP2_ID:
        index_base = ref_info["fasta_path"]
    else:
        index_base = index_path(temp_dir, assembly_id)
    rg = f"@RG\\tID:{sample_name}\\tSM:{sample_name}"

    groups = []
    if sample["r1"]:
        groups.append((sample["r1"], sample.get("r2")))
    for path in sample.get("singles", []):
        groups.append((path, None))

    if not groups:
        raise RuntimeError("No reads available to map after prefiltering")

    tmp_bams = []
    try:
        bin_path = aligner_bin or aligner
        for g_r1, g_r2 in groups:
            if aligner == MINIMAP2_ID:
                preset = "map-ont" if nanopore else "sr"
                align_cmd = [bin_path, "-ax", preset, "-t", str(threads),
                             "--secondary=no", "-R", rg, index_base, g_r1]
                if g_r2:
                    align_cmd.append(g_r2)
            elif g_r2:
                align_cmd = [bin_path, "mem", "-t", str(threads), "-R", rg,
                             index_base, g_r1, g_r2]
            else:
                align_cmd = [bin_path, "mem", "-t", str(threads), "-R", rg,
                             index_base, g_r1]
            part_bam = bam_path + f".part{len(tmp_bams)}"
            sort_cmd = [SAMTOOLS, "sort", "-@", str(threads)]
            if sort_memory:
                sort_cmd += ["-m", sort_memory]
            sort_cmd += ["-o", part_bam]
            with subprocess.Popen(align_cmd, stdout=subprocess.PIPE,
                                  stderr=subprocess.DEVNULL) as align_proc, \
                    subprocess.Popen(sort_cmd, stdin=align_proc.stdout,
                                     stderr=subprocess.PIPE,
                                     text=True) as sort_proc:
                align_proc.stdout.close()
                _, sort_stderr = sort_proc.communicate()
            if align_proc.returncode != 0 or sort_proc.returncode != 0:
                raise RuntimeError(
                    f"{aligner} align (rc={align_proc.returncode}) / samtools sort "
                    f"(rc={sort_proc.returncode}) failed: {sort_stderr[:500]}"
                )
            tmp_bams.append(part_bam)

        if len(tmp_bams) == 1:
            os.replace(tmp_bams[0], bam_path)
        elif tmp_bams:
            subprocess.run(
                [SAMTOOLS, "merge", "-f", "-@", str(threads), bam_path]
                + tmp_bams,
                check=True, capture_output=True)
        else:
            raise RuntimeError("No reads available to map after prefiltering")

        subprocess.run([SAMTOOLS, "index", "-@", str(threads), bam_path],
                       check=True, capture_output=True)

        flagstat_result = subprocess.run(
            [SAMTOOLS, "flagstat", "-@", str(threads), bam_path],
            capture_output=True, text=True, check=True)
        primary_total, _ = parse_flagstat(flagstat_result.stdout)

        with open(bam_path + ".aligner", "w") as f:
            f.write(f"aligner={aligner}\ntrimmed={trim_state}\n"
                    f"platform={'nanopore' if nanopore else 'illumina'}\n"
                    f"total_reads={primary_total}\n")

        if prefilter_config["enabled"]:
            marker = bam_path + ".prefiltered"
            with open(marker, "w") as f:
                f.write(f"min_entropy={prefilter_config['min_entropy']}\n"
                        f"min_dust={prefilter_config['min_dust']}\n")

        duplicates_removed = 0
        if dedup_config["enabled"]:
            duplicates_removed = dedup_bam(bam_path, dedup_config, temp_dir)

        if filter_config["enabled"]:
            filter_bam(bam_path, filter_config, primary_total)

        result = summarize_bam(bam_path, stats_path, sample_name, pathogen,
                               assembly_id, ref_info["fasta_path"], ref_len,
                               threads, skipped=False,
                               total_override=primary_total)
        if dedup_config["enabled"]:
            result["duplicates_removed"] = duplicates_removed
        return result

    except subprocess.CalledProcessError as e:
        error_msg = f"Command failed with return code {e.returncode}"
        if e.stdout:
            error_msg += f"\nstdout: {e.stdout[:500]}"
        if e.stderr:
            error_msg += f"\nstderr: {e.stderr[:500]}"
        return {
            "sample": sample_name,
            "pathogen": pathogen,
            "assembly_id": assembly_id,
            "fasta_path": ref_info["fasta_path"],
            "error": error_msg,
            "skipped": False
        }

    except Exception as e:
        return {
            "sample": sample_name,
            "pathogen": pathogen,
            "assembly_id": assembly_id,
            "fasta_path": ref_info["fasta_path"],
            "error": f"Unexpected error: {type(e).__name__}: {str(e)}",
            "skipped": False
        }

    finally:
        for tb in tmp_bams:
            if os.path.exists(tb):
                os.remove(tb)


def process_sample_ref_combination(args):
    (sample, ref, temp_dir, output_dir, threads,
     filter_config, prefilter_config, dedup_config, aligner,
     trim_state, aligner_bin, nanopore, sort_memory) = args
    return map_sample(sample, ref, temp_dir, output_dir, threads, filter_config,
                      prefilter_config, dedup_config, aligner, trim_state,
                      aligner_bin, nanopore, sort_memory)


def main():
    config_parser = argparse.ArgumentParser(add_help=False)
    config_parser.add_argument("--config", default=None,
                               help="JSON configuration file path")
    config_args, _ = config_parser.parse_known_args()
    try:
        run_config, config_file = mapping_config.load_config(config_args.config)
    except (OSError, ValueError) as e:
        raise SystemExit(f"Invalid mapping configuration: {e}")

    parser = argparse.ArgumentParser(description="Map samples against reference genomes")
    parser.add_argument("--config", default=None,
                        help="JSON configuration file path [default: %(default)s]")
    parser.add_argument("--refs-dir", default=REFS_DIR, help="References directory")
    parser.add_argument("--samples-dir", default=SAMPLES_DIR, help="Samples base directory")
    parser.add_argument("--output", default=OUTPUT_FILE, help="Output TSV file")
    parser.add_argument("--output-dir", default="mapping_results", help="BAM output directory")
    parser.add_argument("--workers", type=int, default=run_config["workers"],
                        help="Number of parallel workers [default: %(default)s]")
    parser.add_argument("--threads", type=int, default=run_config["threads"],
                        help="Total threads across all workers [default: %(default)s]")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--log-file", default=LOG_FILE, help="Log file path")
    parser.add_argument("--msamtools", default=MSAMTOOLS,
                        help="Path to msamtools binary")
    parser.add_argument("--min-len", type=int, default=MIN_LEN,
                        help="Minimum alignment length (bp) [default: %(default)s]")
    parser.add_argument("--min-identity", type=float, default=MIN_IDENTITY,
                        help="Minimum sequence identity, 0-1 [default: %(default)s]")
    parser.add_argument("--min-query-aligned", type=float, default=MIN_QUERY_ALIGNED,
                        help="Minimum fraction of the query aligned, 0-1 [default: %(default)s]")
    parser.add_argument("--no-filter", action="store_true", default=None,
                        help="Disable msamtools length/identity filtering")
    parser.add_argument("--prinseq", default=PRINSEQ,
                        help="Path to prinseq++ binary")
    parser.add_argument("--no-prefilter", action="store_true", default=None,
                        help="Disable prinseq++ entropy/dust read prefiltering")
    parser.add_argument("--min-entropy", type=float, default=MIN_ENTROPY,
                        help="Minimum sequence entropy for prefilter, 0-1 "
                             "[default: %(default)s]")
    parser.add_argument("--min-dust", type=float, default=MIN_DUST,
                        help="Minimum dust score for prefilter, 0-1 "
                             "[default: %(default)s]")
    parser.add_argument("--prefilter-gz", action="store_true",
                        help="Write prefiltered FASTQs gzip-compressed")
    parser.add_argument("--picard", default=PICARD,
                        help="Path to picard binary")
    parser.add_argument("--no-dedup", action="store_true", default=None,
                        help="Disable Picard MarkDuplicates duplicate removal")
    parser.add_argument("--dedup-heap", default=run_config["dedup_heap"],
                        help="JVM heap for MarkDuplicates, e.g. 4g [default: %(default)s]")
    parser.add_argument("--sort-memory", default=run_config["sort_memory"],
                        help="Max memory per samtools sort thread, e.g. 256M")
    parser.add_argument("--aligner", default=BWA_ID,
                        help="Aligner for mapping: 'bwa', 'minimap2', or a full "
                             "path to the aligner binary (detected by basename) "
                             "[default: %(default)s]")
    parser.add_argument("--trimmomatic", default=TRIMMOMATIC,
                        help="Path to trimmomatic binary")
    parser.add_argument("--no-trim", action="store_true", default=None,
                        help="Disable Trimmomatic quality trimming")
    parser.add_argument("--trim-phred", type=int, choices=[33, 64],
                        default=TRIM_PHRED,
                        help="Phred offset for trimming [default: %(default)s]")
    parser.add_argument("--trim-leading", type=int, default=TRIM_LEADING,
                        help="Trim bases below this quality from the start "
                             "[default: %(default)s]")
    parser.add_argument("--trim-trailing", type=int, default=TRIM_TRAILING,
                        help="Trim bases below this quality from the end "
                             "[default: %(default)s]")
    parser.add_argument("--trim-minlen", type=int, default=TRIM_MINLEN,
                        help="Drop reads shorter than this after trimming "
                             "[default: %(default)s]")
    parser.add_argument("--trim-sliding-window", type=str, default=TRIM_SLIDING_WINDOW,
                        help="Trimmomatic sliding window (WINDOW:QUAL, e.g. 5:20) "
                             "[default: %(default)s]")
    parser.add_argument("--nanopore", action="store_true",
                        help="Nanopore long-read mode: NanoFilt QC, single-end "
                             "reads, minimap2 -ax map-ont")
    parser.add_argument("--nanofilt", default=NANOFILT,
                        help="Path to NanoFilt binary")
    parser.add_argument("--nanofilt-quality", type=int, default=NANOFILT_QUALITY,
                        help="Minimum average read quality for NanoFilt "
                             "[default: %(default)s]")
    parser.add_argument("--nanofilt-minlen", type=int, default=NANOFILT_MINLEN,
                        help="Minimum read length for NanoFilt [default: %(default)s]")
    parser.add_argument("--nanofilt-headcrop", type=int, default=NANOFILT_HEADCROP,
                        help="Trim n nucleotides from the start of each read "
                             "(0 disables) [default: %(default)s]")
    parser.add_argument("--nanofilt-tailcrop", type=int, default=NANOFILT_TAILCROP,
                        help="Trim n nucleotides from the end of each read "
                             "(0 disables) [default: %(default)s]")
    args = parser.parse_args()

    no_filter = run_config["no_filter"] if args.no_filter is None else args.no_filter
    no_prefilter = (run_config["no_prefilter"] if args.no_prefilter is None
                    else args.no_prefilter)
    no_dedup = run_config["no_dedup"] if args.no_dedup is None else args.no_dedup
    no_trim = run_config["no_trim"] if args.no_trim is None else args.no_trim

    if args.workers < 1:
        parser.error("--workers must be >= 1")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if not 0 < args.min_identity <= 1:
        parser.error("--min-identity must be in (0, 1]")
    if not 0 < args.min_query_aligned <= 1:
        parser.error("--min-query-aligned must be in (0, 1]")
    if args.min_len < 0:
        parser.error("--min-len must be >= 0")
    if not 0 < args.min_entropy <= 1:
        parser.error("--min-entropy must be in (0, 1]")
    if not 0 < args.min_dust <= 1:
        parser.error("--min-dust must be in (0, 1]")
    if args.trim_minlen < 1:
        parser.error("--trim-minlen must be >= 1")
    if args.trim_leading < 0 or args.trim_trailing < 0:
        parser.error("--trim-leading/--trim-trailing must be >= 0")
    try:
        sw_win, sw_qual = (int(p) for p in args.trim_sliding_window.split(":", 1))
    except ValueError:
        parser.error("--trim-sliding-window must be in WINDOW:QUAL form (e.g. 5:20)")
    if sw_win < 1 or sw_qual < 1:
        parser.error("--trim-sliding-window WINDOW and QUAL must be >= 1")
    if args.nanopore:
        if args.nanofilt_quality < 0:
            parser.error("--nanofilt-quality must be >= 0")
        if args.nanofilt_minlen < 1:
            parser.error("--nanofilt-minlen must be >= 1")
        if args.nanofilt_headcrop < 0:
            parser.error("--nanofilt-headcrop must be >= 0")
        if args.nanofilt_tailcrop < 0:
            parser.error("--nanofilt-tailcrop must be >= 0")

    try:
        aligner_id, aligner_bin = resolve_aligner(args.aligner)
    except ValueError as e:
        parser.error(str(e))
    if args.nanopore and aligner_id != MINIMAP2_ID:
        aligner_id, aligner_bin = MINIMAP2_ID, MINIMAP2

    filter_config = {
        "enabled": not no_filter,
        "msamtools": args.msamtools,
        "min_len": args.min_len,
        "min_identity": args.min_identity,
        "min_identity_pct": round(args.min_identity * 100),
        "min_query_aligned": args.min_query_aligned,
        "min_query_pct": round(args.min_query_aligned * 100),
    }

    prefilter_config = {
        "enabled": not no_prefilter,
        "prinseq": args.prinseq,
        "min_entropy": args.min_entropy,
        "min_dust": args.min_dust,
        "gz": args.prefilter_gz,
        "dir": os.path.join(args.output_dir, "prefiltered"),
    }

    trim_config = {
        "enabled": not no_trim,
        "trimmomatic": args.trimmomatic,
        "phred": args.trim_phred,
        "leading": args.trim_leading,
        "trailing": args.trim_trailing,
        "minlen": args.trim_minlen,
        "sliding_window": args.trim_sliding_window,
        "dir": os.path.join(args.output_dir, "trimmed"),
    }

    if args.nanopore:
        trim_config["enabled"] = False

    nanofilt_config = {
        "enabled": args.nanopore,
        "nanofilt": args.nanofilt,
        "quality": args.nanofilt_quality,
        "minlen": args.nanofilt_minlen,
        "headcrop": args.nanofilt_headcrop,
        "tailcrop": args.nanofilt_tailcrop,
        "dir": os.path.join(args.output_dir, "nanofiltered"),
    }

    dedup_config = {
        "enabled": not no_dedup,
        "picard": args.picard,
        "heap": args.dedup_heap,
    }

    per_combo_threads = max(1, args.threads // args.workers)

    logger = setup_logging(args.log_file, args.verbose)
    logger.info("Starting mapping pipeline")
    logger.info(f"Configuration file: {config_file}")

    logger.debug(
        f"Config: refs={args.refs_dir} samples={args.samples_dir} "
        f"workers={args.workers} total_threads={args.threads} "
        f"per_combo_threads={per_combo_threads} aligner={aligner_id} "
        f"sort_memory={args.sort_memory} "
        f"aligner_bin={aligner_bin} "
        f"trim={trim_config} filter={filter_config} "
        f"prefilter={prefilter_config} dedup={dedup_config} "
        f"nanofilt={nanofilt_config}"
    )

    logger.info(f"Checking dependencies (aligner={aligner_id})...")
    check_dependencies(msamtools_path=args.msamtools,
                       need_msamtools=filter_config["enabled"],
                       prinseq_path=args.prinseq,
                       need_prinseq=prefilter_config["enabled"],
                       picard_path=args.picard,
                       need_picard=dedup_config["enabled"],
                       aligner=aligner_id,
                       aligner_bin=aligner_bin,
                       trimmomatic_path=trim_config["trimmomatic"],
                       need_trimmomatic=trim_config["enabled"],
                       nanofilt_path=nanofilt_config["nanofilt"],
                       need_nanofilt=nanofilt_config["enabled"])
    logger.info("  All required tools found.")
    if trim_config["enabled"]:
        logger.info(f"Trimming reads (LEADING:{trim_config['leading']} "
                    f"TRAILING:{trim_config['trailing']} "
                    f"SLIDINGWINDOW:{trim_config['sliding_window']} "
                    f"MINLEN:{trim_config['minlen']}, "
                    f"phred{trim_config['phred']})")
    if filter_config["enabled"]:
        logger.info(f"Filtering alignments (min_len={filter_config['min_len']}, "
                    f"min_identity={filter_config['min_identity_pct']}%, "
                    f"min_query_aligned={filter_config['min_query_pct']}%)")
    if prefilter_config["enabled"]:
        logger.info(f"Prefiltering reads (min_entropy={prefilter_config['min_entropy']}, "
                    f"min_dust={prefilter_config['min_dust']}, "
                    f"gzip={prefilter_config['gz']})")
    if dedup_config["enabled"]:
        logger.info(f"Removing duplicates (MarkDuplicates, heap={dedup_config['heap']})")
    if nanofilt_config["enabled"]:
        logger.info(f"NanoFilt filtering reads (quality={nanofilt_config['quality']}, "
                    f"minlen={nanofilt_config['minlen']}, "
                    f"headcrop={nanofilt_config['headcrop']}, "
                    f"tailcrop={nanofilt_config['tailcrop']})")
    if args.nanopore:
        logger.info(f"Nanopore mode: minimap2 -ax map-ont (aligner={aligner_id})")

    logger.info("Finding references...")
    references = find_references(args.refs_dir)
    logger.info(f"Found {len(references)} reference genomes:")
    for r in references:
        logger.info(f"  - {r['pathogen']} ({r['assembly_id']})")

    logger.info("Finding samples...")
    samples = find_samples(args.samples_dir, nanopore=args.nanopore)
    logger.info(f"Found {len(samples)} samples:")
    for s in samples:
        logger.info(f"  - {s['name']} ({'paired-end' if s['is_paired'] else 'single-end'})")

    if not references:
        raise ValueError(f"No reference genomes found in {args.refs_dir}")

    logger.info("Validating references...")
    for ref in references:
        if not os.path.exists(ref["fasta_path"]):
            raise FileNotFoundError(f"Reference FASTA not found: {ref['fasta_path']}")
        ref["ref_len"] = fasta_total_length(ref["fasta_path"])
    logger.info("  All reference files valid.")

    if not samples:
        raise ValueError(f"No samples found in {args.samples_dir}")

    logger.info("Validating samples...")
    for sample in samples:
        if not os.path.exists(sample["r1"]):
            raise FileNotFoundError(f"R1 FASTQ not found: {sample['r1']}")
        if sample["is_paired"] and sample["r2"] and not os.path.exists(sample["r2"]):
            raise FileNotFoundError(f"R2 FASTQ not found: {sample['r2']}")
    logger.info("  All sample files valid.")

    os.makedirs(args.output_dir, exist_ok=True)

    temp_dir = tempfile.mkdtemp(prefix=f"index_{aligner_id}_")
    logger.info(f"Using temp directory: {temp_dir}")
    logger.info(f"Threads per combination: {per_combo_threads}")

    try:
        logger.info(f"Building {aligner_id} indexes...")
        for ref in references:
            try:
                build_index(ref, temp_dir, aligner_id, aligner_bin)
                logger.info(f"  Indexed: {ref['pathogen']} ({ref['assembly_id']})")
            except RuntimeError as e:
                logger.error(f"  ERROR indexing {ref['pathogen']}/{ref['assembly_id']}: {e}")
                raise

        combinations = [(s, r, temp_dir, args.output_dir, per_combo_threads,
                         filter_config, prefilter_config, dedup_config,
                         aligner_id, s.get("trim_fingerprint", ""), aligner_bin,
                         args.nanopore, args.sort_memory)
                        for s in samples for r in references]
        total_combinations = len(combinations)

        logger.info(f"Processing {total_combinations} sample × reference combinations...")

        results = []
        completed = 0

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            if nanofilt_config["enabled"]:
                logger.info("Filtering nanopore reads with NanoFilt...")
                filtered = []
                nf_futures = {executor.submit(nanofilt_sample, s,
                                              nanofilt_config, per_combo_threads): s
                              for s in samples}
                for future in as_completed(nf_futures):
                    sample = future.result()
                    filtered.append(sample)
                    info = sample.get("nanofilt_log")
                    if info:
                        logger.info(
                            f"  {sample['name']}: {info['input_reads']} input reads -> "
                            f"{info['passed_reads']} passed NanoFilt "
                            f"({info['removed_reads']} removed)")
                    else:
                        logger.info(f"  {sample['name']}: reused cached nanofiltered reads")
                samples = filtered
                combinations = [(s, r, temp_dir, args.output_dir,
                                 per_combo_threads, filter_config,
                                 prefilter_config, dedup_config,
                                 aligner_id, s.get("trim_fingerprint", ""),
                                 aligner_bin, args.nanopore, args.sort_memory)
                                for s in samples for r in references]

            if trim_config["enabled"]:
                logger.info("Trimming reads with Trimmomatic...")
                trimmed = []
                trim_futures = {executor.submit(trim_sample, s, trim_config,
                                                per_combo_threads): s
                                for s in samples}
                for future in as_completed(trim_futures):
                    sample = future.result()
                    trimmed.append(sample)
                    info = sample.get("trim_log")
                    if info:
                        logger.info(
                            f"  {sample['name']}: {info['input_reads']} input reads -> "
                            f"{info['passed_reads']} passed trimming "
                            f"({info['removed_reads']} removed)")
                    else:
                        logger.info(f"  {sample['name']}: reused cached trimmed reads")
                samples = trimmed
                combinations = [(s, r, temp_dir, args.output_dir,
                                 per_combo_threads, filter_config,
                                 prefilter_config, dedup_config,
                                  aligner_id, s.get("trim_fingerprint", ""), aligner_bin,
                                  args.nanopore, args.sort_memory)
                                 for s in samples for r in references]


            if prefilter_config["enabled"]:
                logger.info("Prefiltering reads with prinseq++ (entropy/dust)...")
                prefiltered = []
                pf_futures = {executor.submit(prefilter_sample, s,
                                              prefilter_config, per_combo_threads): s
                              for s in samples}
                for future in as_completed(pf_futures):
                    sample = future.result()
                    prefiltered.append(sample)
                    info = sample.get("prefilter_log")
                    if info:
                        logger.info(
                            f"  {sample['name']}: {info['input_reads']} input reads -> "
                            f"{info['passed_reads']} passed prefilter "
                            f"({info['removed_reads']} removed)")
                    else:
                        logger.info(f"  {sample['name']}: reused cached prefiltered reads")
                samples = prefiltered
                combinations = [(s, r, temp_dir, args.output_dir,
                                 per_combo_threads, filter_config,
                                 prefilter_config, dedup_config,
                                  aligner_id, s.get("trim_fingerprint", ""), aligner_bin,
                                  args.nanopore, args.sort_memory)
                                 for s in samples for r in references]


            futures = {executor.submit(process_sample_ref_combination, combo): combo
                       for combo in combinations}

            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                completed += 1

                if result.get("error"):
                    logger.error(
                        f"MAPPING ERROR - Sample: {result['sample']} | "
                        f"Pathogen: {result['pathogen']} | "
                        f"Assembly: {result['assembly_id']} | "
                        f"Error: {result['error']}"
                    )

                status = "SKIPPED" if result.get("skipped") else "DONE" if "error" not in result else "ERROR"
                dup_info = ""
                if "duplicates_removed" in result:
                    dup_info = f" (duplicates removed: {result['duplicates_removed']})"
                logger.info(f"[{completed}/{total_combinations}] {result['sample']} vs {result['pathogen']}: {status}{dup_info}")

        logger.info("Collecting statistics...")

        stats_data = []
        for r in results:
            stats_data.append({
                "sample": r.get("sample", ""),
                "pathogen": r.get("pathogen", ""),
                "assembly_id": r.get("assembly_id", ""),
                "fasta_path": r.get("fasta_path", ""),
                "total_reads": r.get("total_reads", 0),
                "mapped_reads": r.get("mapped_reads", 0),
                "unmapped_reads": r.get("unmapped_reads", 0),
                "mapped_pct": r.get("mapped_pct", 0),
                "coverage": r.get("coverage"),
                "avg_depth": r.get("avg_depth"),
                "depth_covered": r.get("depth_covered"),
                "ref_length": r.get("ref_length"),
                "duplicates_removed": r.get("duplicates_removed", 0),
                "bam_path": r.get("bam_path", ""),
                "error": r.get("error", "")
            })

        df = pd.DataFrame(stats_data)

        column_order = ["sample", "pathogen", "assembly_id", "fasta_path", "total_reads",
                        "mapped_reads", "unmapped_reads", "mapped_pct", "coverage",
                        "avg_depth", "depth_covered", "ref_length", "duplicates_removed",
                        "bam_path", "error"]

        existing_cols = [c for c in column_order if c in df.columns]
        df = df[existing_cols]

        df.to_csv(args.output, sep="\t", index=False)

        logger.info(f"Results saved to: {args.output}")
        logger.info(f"BAM files saved to: {args.output_dir}")

        logger.info("=== Summary ===")
        logger.info(f"Total combinations: {total_combinations}")
        logger.info(f"Successful mappings: {len([r for r in results if not r.get('skipped') and 'error' not in r])}")
        logger.info(f"Resumed from existing BAMs (skipped): {len([r for r in results if r.get('skipped') and 'error' not in r])}")

        errors = [r for r in results if r.get("error")]
        if errors:
            logger.error(f"=== {len(errors)} MAPPING ERRORS ===")
            for err in errors:
                logger.error(f"  {err['sample']} vs {err['pathogen']}: {err['error'][:500]}")

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
        logger.info(f"Temp directory cleaned up: {temp_dir}")
        logger.info("Pipeline completed")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(f"Fatal error: {type(e).__name__}: {str(e)}")
        sys.exit(1)
