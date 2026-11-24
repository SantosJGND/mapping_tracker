"""Streamlit companion app for map_references.py results.

Run from the mapping_tracker directory (relative paths resolve from there):
    streamlit run app/app.py

Reads mapping_statistics.tsv when available, or builds the table live from
completed combinations while the pipeline is still running.

The Run tab can start a detached job (map_references.py, then registration of the
samples used, then extract_remap.py). The job writes into the output dir
(job.json, job.log, extract/<pathogen>/). The other tabs only read outputs.
"""
import glob
import os
import subprocess
import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
import browser  # noqa: E402
import depth  # noqa: E402
import jobs  # noqa: E402
import loader  # noqa: E402
import views  # noqa: E402
from map_references import find_references  # noqa: E402

st.set_page_config(page_title="mapping_tracker", layout="wide")

with st.sidebar:
    st.header("Inputs")
    refs_dir = st.text_input("References dir", "references")
    output_dir = st.text_input("Output dir", "mapping_results")
    tsv_path = st.text_input("Statistics TSV", "mapping_statistics.tsv")
    log_path = st.text_input("Log file", "mapping_log.txt")
    source = st.selectbox(
        "Source", ["auto", "tsv", "live"],
        help="auto: TSV if it exists and no run is in progress, otherwise live.")
    live = st.checkbox("Auto-refresh", value=True)
    interval = st.number_input("Refresh interval (s)", min_value=2,
                               max_value=300, value=5)


@st.cache_data(show_spinner="Loading results...")
def _load(source, tsv_path, refs_dir, output_dir, log_path, signature):
    # `signature` is part of the cache key (no leading underscore on purpose):
    # it changes whenever inputs change on disk, invalidating the cache.
    return loader.load_table(source, tsv_path, refs_dir, output_dir, log_path)


def _show_table(frame, **kwargs):
    st.dataframe(frame, hide_index=True, **kwargs)


def _run_tab(refs_dir, output_dir, tsv_path, log_path):
    st.subheader("Run mapping + extract")
    st.caption("Maps all samples against all references, then excises and remaps the "
               "samples that mapped to the chosen reference. One job runs at a time.")
    samples_dir = st.text_input("Samples dir", "samples", key="run_samples")
    map_script = st.text_input("map_references.py", str(jobs.DEFAULT_MAP_SCRIPT),
                               key="run_map_script")
    extract_script = st.text_input("extract_remap.py", str(jobs.DEFAULT_EXTRACT_SCRIPT),
                                   key="run_extract_script")
    c1, c2 = st.columns(2)
    workers = c1.number_input("Workers", min_value=1, max_value=64, value=4, key="run_workers")
    threads = c2.number_input("Threads", min_value=1, max_value=256,
                              value=os.cpu_count() or 4, key="run_threads")
    t1, t2, t3 = st.columns(3)
    trim = t1.checkbox("Trimmomatic trimming", True, key="run_trim")
    prefilter = t2.checkbox("prinseq++ prefilter", True, key="run_prefilter")
    dedup = t3.checkbox("Picard dedup", True, key="run_dedup")

    try:
        refs = find_references(refs_dir)
    except ValueError as e:
        st.error(str(e))
        refs = []
    if not refs:
        st.info("No references found in the References dir (sidebar).")
        return
    labels = [f"{r['pathogen']} / {r['assembly_id']}" for r in refs]
    choice = st.selectbox("Reference for extract", labels, key="run_ref")
    ref = refs[labels.index(choice)]

    job = jobs.read_job(output_dir)
    active = jobs.job_active(output_dir)
    if job:
        st.markdown(f"**Last job:** {job.get('status')} - stage: {job.get('stage')} "
                    f"(started {job.get('started')}, pid {job.get('pid')})")
        if job.get("error"):
            st.error(job["error"])
        tail = jobs.job_log_tail(output_dir)
        if tail:
            st.code(tail, language="text")

    if st.button("Cancel job", disabled=not active, key="run_cancel"):
        if jobs.cancel_job(output_dir):
            st.warning("Job cancelled. Partial outputs are left in place; "
                       "a re-run resumes from existing BAMs.")
            active = False

    if st.button("Start job", disabled=active, key="run_start"):
        try:
            spec = {
                "map_cmd": jobs.build_map_cmd(
                    os.path.abspath(map_script), os.path.abspath(refs_dir),
                    os.path.abspath(samples_dir), os.path.abspath(tsv_path),
                    os.path.abspath(output_dir), os.path.abspath(log_path),
                    int(workers), int(threads), trim, prefilter, dedup),
                "extract_script": os.path.abspath(extract_script),
                "ref_fasta": os.path.abspath(ref["fasta_path"]),
                "samples_dir": os.path.abspath(samples_dir),
                "pathogen": ref["pathogen"],
                "assembly_id": ref["assembly_id"],
                "threads": int(threads),
            }
            jobs.start_job(spec, output_dir)
            st.success("Job started. Status above updates on refresh.")
        except RuntimeError as e:
            st.error(str(e))

    ext = jobs.extract_dir(os.path.abspath(output_dir), ref["pathogen"])
    registered = jobs.read_registered(ext)
    st.markdown(f"**Samples registered for extract ({ref['pathogen']}):**")
    if registered:
        _show_table([{k: r[k] for k in ("name", "r1", "r2", "bam")} for r in registered])
    else:
        st.caption("None yet.")


@st.cache_resource(show_spinner=False)
def _browser_server(outdir):
    # One range server per extract_remap output dir, started on first use.
    return browser.start_server(outdir)


def _browser_signature(outdir, picked):
    paths = [
        os.path.join(outdir, "subset_ref.fa"),
        os.path.join(outdir, "subset_ref.fa.fai"),
        os.path.join(outdir, "coverage_intervals.bed"),
    ]
    for sample in picked:
        bam = os.path.join(outdir, sample, "remapped.sorted.bam")
        paths.extend([bam, bam + ".bai"])
    return tuple(
        (path, os.stat(path).st_mtime_ns, os.stat(path).st_size)
        for path in paths if os.path.isfile(path)
    )


@st.cache_data(show_spinner=False)
def _browser_view(outdir, picked, base, signature):
    intervals = browser.read_intervals(outdir)
    contigs = browser.read_contigs(outdir)
    junction_path = browser.write_junctions(outdir, intervals)
    flags_error = None
    try:
        flags_path = browser.write_read_flags(outdir, list(picked), intervals)
    except subprocess.CalledProcessError as e:
        flags_path = None
        flags_error = e.stderr or str(e)
    config = browser.build_config(base, list(picked), contigs[0],
                                  junction_path, flags_path)
    return browser.viewer_html(base, config), flags_error


def _extract_outputs(output_dir):
    """Extract output dirs (<output_dir>/extract/<pathogen>) that hold a subset reference."""
    found = glob.glob(os.path.join(output_dir, jobs.EXTRACT_DIR, "*", "subset_ref.fa.fai"))
    return sorted(os.path.dirname(p) for p in found)


def _browser_tab(output_dir):
    options = _extract_outputs(output_dir)
    manual = st.text_input("extract_remap output dir (overrides the choice below)", "",
                           key="br_outdir")
    if manual:
        outdir = manual
    elif options:
        outdir = st.selectbox("Extract output", options, key="br_pick")
    else:
        st.info("No extract output yet. Run a job in the Run tab, or enter a directory above.")
        return
    if not os.path.isfile(os.path.join(outdir, "subset_ref.fa.fai")):
        st.info("No subset_ref.fa.fai in this directory yet.")
        return
    samples = browser.list_samples(outdir)
    if not samples:
        st.info("No <sample>/remapped.sorted.bam found in this directory.")
        return
    picked = st.multiselect(f"Samples (max {browser.MAX_TRACKS})", samples,
                            default=samples[:browser.MAX_TRACKS], key="br_samples")
    if len(picked) > browser.MAX_TRACKS:
        st.warning(f"Only the first {browser.MAX_TRACKS} selected samples are shown.")
        picked = picked[:browser.MAX_TRACKS]
    if not picked:
        st.info("Select at least one sample.")
        return

    server = _browser_server(outdir)
    base = f"http://127.0.0.1:{server.server_address[1]}"
    intervals = browser.read_intervals(outdir)
    contigs = browser.read_contigs(outdir)
    picked = tuple(picked)
    html, flags_error = _browser_view(outdir, picked, base,
                                      _browser_signature(outdir, picked))
    if flags_error:
        st.error(f"Junction read scan failed: {flags_error}")
    st.iframe(html, height=640, alt="Genome browser for excised reference and selected samples")
    st.caption("Coordinates are subset coordinates (subset_ref.fa). "
               "Red track: junctions between joined intervals. "
               "Orange track: reads whose alignment spans a junction.")

    with st.expander("Subset -> original coordinate lookup"):
        contig = st.selectbox("Subset contig", contigs, key="br_contig")
        pos = st.number_input("Subset position (0-based)", min_value=0, value=0,
                              step=1, key="br_pos")
        hit = browser.subset_to_original(intervals, contig, int(pos))
        if hit is None:
            st.write("Position is not inside any excised interval.")
        else:
            st.write(f"Original position: `{hit[0]}:{hit[1]}` (0-based)")


def render():
    signature = loader.signature(tsv_path, log_path, output_dir)
    try:
        df, mode, status = _load(source, tsv_path, refs_dir, output_dir,
                                 log_path, signature)
    except Exception as e:  # shown to the user, app keeps running
        st.error(f"Could not load results: {type(e).__name__}: {e}")
        return

    c1, c2, c3 = st.columns(3)
    c1.metric("Data source", mode.upper())
    c2.metric("Run status", status)
    c3.metric("Completed combinations", len(df))

    progress = loader.read_progress(log_path)
    if progress and status == "running":
        done, total = progress
        if total:
            st.progress(min(done / total, 1.0),
                        text=f"Mapping progress: {done}/{total} combinations")

    ok = df[df["error"].fillna("") == ""]
    if ok.empty:
        st.info("No completed combinations with results yet.")

    tab_refs, tab_heat, tab_sample, tab_pair, tab_browser, tab_run, tab_warn = st.tabs(
        ["References", "Heatmap", "Sample", "Pair", "Browser", "Run", "Warnings"],
        key="main_tabs", on_change="rerun")

    if tab_refs.open:
        with tab_refs:
            refs = views.references_table(df)
            _show_table(refs)
            st.plotly_chart(views.references_figure(refs))

    if tab_heat.open:
        with tab_heat:
            metric = st.selectbox("Metric", list(views.METRICS),
                                  format_func=views.METRICS.get, key="heat_metric")
            fig = views.heatmap_figure(ok, metric)
            if fig is None:
                st.info("Nothing to show yet.")
            else:
                st.plotly_chart(fig)

    if tab_sample.open:
        with tab_sample:
            if not ok.empty:
                sample = st.selectbox("Sample", sorted(ok["sample"].unique()),
                                      key="sample_pick")
                sdf = ok[ok["sample"] == sample]
                metric_s = st.selectbox("Metric", list(views.METRICS),
                                        format_func=views.METRICS.get,
                                        key="sample_metric")
                st.plotly_chart(views.sample_figure(sdf, metric_s))
                _show_table(sdf[["pathogen", "assembly_id", "total_reads",
                                 "mapped_reads", "mapped_pct", "coverage",
                                 "avg_depth", "depth_covered", "duplicates_removed"]])

    if tab_pair.open:
        with tab_pair:
            if not ok.empty:
                sample = st.selectbox("Sample", sorted(ok["sample"].unique()),
                                      key="pair_sample")
                pairs = ok[ok["sample"] == sample].reset_index(drop=True)
                labels = views.ref_label(pairs).tolist()
                choice = st.selectbox("Reference", labels, key="pair_ref")
                row = pairs.iloc[labels.index(choice)]

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Mapped (%)", f"{row['mapped_pct']:.2f}")
                m2.metric("Coverage (%)", _fmt(row["coverage"]))
                m3.metric("Average depth", _fmt(row["avg_depth"]))
                m4.metric("Duplicates removed", int(row["duplicates_removed"]))
                st.caption(f"BAM: {row['bam_path']}")

                bam = str(row["bam_path"])
                stats_file = (bam[:-len(".bam")] + ".stats.txt"
                              if bam.endswith(".bam") else "")
                if stats_file and os.path.exists(stats_file):
                    _show_table(views.stats_table(loader.parse_stats_sn(stats_file)))
                else:
                    st.warning("samtools stats file not found for this pair.")

                if st.checkbox("Show per-position coverage (binned, cached)",
                               key="pair_depth"):
                    if not row["bam_exists"]:
                        st.warning("BAM not found; per-position coverage unavailable.")
                    else:
                        try:
                            with st.spinner("Computing binned depth..."):
                                binned = depth.binned_depth(bam, depth.DEFAULT_BINS)
                            st.plotly_chart(views.depth_figure(
                                binned, f"{sample} vs {choice}"))
                        except Exception as e:
                            st.error(f"Depth computation failed: {e}")

    if tab_browser.open:
        with tab_browser:
            _browser_tab(output_dir)

    if tab_run.open:
        with tab_run:
            _run_tab(refs_dir, output_dir, tsv_path, log_path)

    if tab_warn.open:
        with tab_warn:
            missing = df[~df["bam_exists"]]
            if missing.empty:
                st.success("All BAM paths exist.")
            else:
                st.warning(f"{len(missing)} row(s) reference a missing BAM.")
                _show_table(missing[["sample", "pathogen", "assembly_id", "bam_path"]])
            errored = df[df["error"].fillna("") != ""]
            if not errored.empty:
                st.error(f"{len(errored)} row(s) with errors.")
                _show_table(errored[["sample", "pathogen", "assembly_id", "error"]])


def _fmt(value):
    return "n/a" if value is None or value != value else f"{value:.2f}"


if live:
    st.fragment(run_every=int(interval))(render)()
else:
    render()
