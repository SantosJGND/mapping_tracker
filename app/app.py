"""Streamlit app for uploading data, running mapping_tracker, and exploring results.

Run from the mapping_tracker directory:
    streamlit run app/app.py
"""
import glob
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

APP_DIR = Path(__file__).resolve().parent
ROOT_DIR = APP_DIR.parent
for _path in (str(ROOT_DIR), str(APP_DIR)):
    if _path not in sys.path:
        sys.path.insert(0, _path)
import browser  # noqa: E402
import deploy_queue  # noqa: E402
import depth  # noqa: E402
import jobs  # noqa: E402
import loader  # noqa: E402
import runs  # noqa: E402
import uploads  # noqa: E402
import views  # noqa: E402
from map_references import find_references, find_samples  # noqa: E402
from mapping_config import load_config  # noqa: E402

MAX_UPLOAD_MB = 5120
RUN_CONFIG, RUN_CONFIG_PATH = load_config()

st.set_page_config(page_title="mapping_tracker", layout="wide")


def _run_label(run_id):
    record = runs.read_run(run_id) or {}
    return f"{record.get('name', run_id)} [{record.get('status', 'unknown')}]"


def _default_index(options, wanted):
    if wanted in options:
        return options.index(wanted)
    return 0


with st.sidebar:
    st.header("Results source")
    source_mode = st.segmented_control("Source", ["Run", "Manual paths"],
                                       default="Run", key="source_mode")
    if source_mode == "Run":
        run_options = [r["run_id"] for r in runs.list_runs()]
        results_run = st.selectbox(
            "Run", run_options, key="results_run",
            format_func=_run_label,
            help="Use data/runs/<run_id> from the Data and Runs tabs.")
        if results_run:
            paths = runs.run_paths(results_run)
            refs_dir = paths["references"]
            output_dir = paths["output"]
            tsv_path = paths["tsv"]
            log_path = paths["log"]
        else:
            refs_dir = "references"
            output_dir = "mapping_results"
            tsv_path = "mapping_statistics.tsv"
            log_path = "mapping_log.txt"
    else:
        results_run = None
        refs_dir = st.text_input("References dir", "references")
        output_dir = st.text_input("Output dir", "mapping_results")
        tsv_path = st.text_input("Statistics TSV", "mapping_statistics.tsv")
        log_path = st.text_input("Log file", "mapping_log.txt")
    source = st.selectbox(
        "Table source", ["auto", "tsv", "live"],
        help="auto: TSV if it exists and no run is in progress, otherwise live.")
    live = st.checkbox("Auto-refresh", value=True)
    interval = st.number_input("Refresh interval (s)", min_value=2,
                               max_value=300, value=5)


@st.cache_data(show_spinner="Loading results...")
def _load(source, tsv_path, refs_dir, output_dir, log_path, signature):
    return loader.load_table(source, tsv_path, refs_dir, output_dir, log_path)


def _show_table(frame, alt, **kwargs):
    st.dataframe(frame, hide_index=True, alt=alt, **kwargs)


def _data_tab():
    st.subheader("Upload data")
    st.caption("Reference genomes as ZIP files; reads as individual FASTQ files. "
               "Plain FASTQ files are gzipped automatically.")
    with st.form("upload_form", border=True):
        name = st.text_input("Run name", key="upload_name",
                             placeholder="Panel run 1")
        ref_uploads = st.file_uploader(
            "Reference ZIP files", type="zip", accept_multiple_files=True,
            max_upload_size=MAX_UPLOAD_MB, key="upload_refs",
            help="Flat FASTA files or an NCBI dataset directory.")
        fastq_uploads = st.file_uploader(
            "FASTQ files", type=["fastq.gz", "fq.gz", "fastq", "fq"],
            accept_multiple_files=True, max_upload_size=MAX_UPLOAD_MB,
            key="upload_fastqs")
        submitted = st.form_submit_button("Upload data", icon=":material/upload:")

    if submitted:
        if not ref_uploads or not fastq_uploads:
            st.error("Upload at least one reference ZIP and one FASTQ file.")
        else:
            try:
                run = uploads.create_run(name, ref_uploads, fastq_uploads)
                st.session_state["last_run_id"] = run["run_id"]
                st.success(f"Uploaded run `{run['run_id']}` with "
                           f"{run['reference_count']} reference(s) and "
                           f"{run['sample_count']} sample(s).")
                st.rerun()
            except Exception as e:
                st.error(f"Upload failed: {type(e).__name__}: {e}")

    st.divider()
    st.subheader("Uploaded datasets")
    records = runs.list_runs()
    if not records:
        st.info("No datasets uploaded yet.")
        return
    options = [r["run_id"] for r in records]
    default = _default_index(options, st.session_state.get("last_run_id"))
    selected = st.selectbox("Dataset", options, index=default, key="data_run",
                            format_func=_run_label)
    record = runs.read_run(selected) or {}
    selected_paths = runs.run_paths(selected)
    c1, c2, c3 = st.columns(3)
    c1.metric("References", record.get("reference_count", 0))
    c2.metric("Samples", record.get("sample_count", 0))
    c3.metric("Status", record.get("status", "unknown"))
    if record.get("error"):
        st.error(record["error"])
    st.code("\n".join(f"{key}: {value}" for key, value in selected_paths.items()),
            language="text")
    if st.checkbox("Show detected files", key="data_detected"):
        try:
            detected_refs = find_references(selected_paths["references"])
            detected_samples = find_samples(selected_paths["samples"])
            _show_table(
                [{"pathogen": r["pathogen"], "assembly_id": r["assembly_id"],
                  "fasta": r["fasta_path"]} for r in detected_refs],
                alt="Detected reference genomes")
            _show_table(
                [{"sample": s["name"], "r1": s["r1"], "r2": s["r2"] or "-"}
                 for s in detected_samples],
                alt="Detected FASTQ samples")
        except Exception as e:
            st.error(f"Could not inspect uploaded files: {type(e).__name__}: {e}")


def _ref_options(refs):
    labels = [f"{r['pathogen']} / {r['assembly_id']}" for r in refs]
    return labels, {label: r for label, r in zip(labels, refs)}


def _run_deploy_form(run_id, disabled):
    prefix = f"queue_{run_id}"
    with st.expander(f"Config: {_run_label(run_id)}", expanded=False):
        c1, c2 = st.columns(2)
        workers = c1.number_input("Workers", min_value=1, max_value=64,
                                  value=int(RUN_CONFIG["workers"]), key=f"{prefix}_workers")
        threads = c2.number_input("Threads", min_value=1, max_value=256,
                                  value=int(RUN_CONFIG["threads"]), key=f"{prefix}_threads")
        t1, t2, t3 = st.columns(3)
        trim = t1.checkbox("Trimmomatic trimming", not RUN_CONFIG["no_trim"], key=f"{prefix}_trim")
        prefilter = t2.checkbox("prinseq++ prefilter", not RUN_CONFIG["no_prefilter"],
                                key=f"{prefix}_prefilter")
        dedup = t3.checkbox("Picard dedup", not RUN_CONFIG["no_dedup"], key=f"{prefix}_dedup")
        st.caption(f"Run defaults source: `{RUN_CONFIG_PATH}`. Extract runs for all references.")
    return {
        "workers": int(workers),
        "threads": int(threads),
        "trim": bool(trim),
        "prefilter": bool(prefilter),
        "dedup": bool(dedup),
    }


def _runs_tab():
    st.subheader("Run mapping + extract")
    st.caption("One job runs at a time across all datasets. Sequential deployment can queue "
               "multiple runs; each run keeps its own config.")
    active = runs.active_run()
    queue = deploy_queue.read_queue()
    queue_active = deploy_queue.queue_active(queue)
    if active:
        st.warning(f"Job running in `{active}`.")
    if queue_active:
        st.info(f"Deployment queue running ({queue.get('status')}, {len(queue.get('items', []))} item(s)).")
    records = runs.list_runs()
    if not records:
        st.info("Upload a dataset in the Data tab first.")
        return
    options = [r["run_id"] for r in records]
    default = _default_index(options, active or st.session_state.get("last_run_id"))
    selected = st.selectbox("Run", options, index=default, key="runs_run",
                            format_func=_run_label)
    paths = runs.run_paths(selected)
    record = runs.read_run(selected) or {}
    c1, c2 = st.columns(2)
    c1.metric("Status", record.get("status", "unknown"))
    c2.metric("Stage", record.get("stage", "-"))
    if record.get("error"):
        st.error(record["error"])

    tail = runs.run_log_tail(selected) or jobs.job_log_tail(paths["output"])
    if tail:
        with st.expander("Job log"):
            st.code(tail, language="text")

    if st.button("Cancel job", disabled=not active or active != selected,
                 key="runs_cancel"):
        if jobs.cancel_job(paths["output"]):
            st.warning("Job cancelled. Partial outputs are left in place.")
            st.rerun()

    try:
        refs = find_references(paths["references"])
    except ValueError as e:
        st.error(str(e))
        refs = []
    if refs:
        st.caption(f"Run defaults: `{RUN_CONFIG_PATH}`. Extract runs for all references.")
        with st.form("start_form", border=True):
            c1, c2 = st.columns(2)
            workers = c1.number_input("Workers", min_value=1, max_value=64,
                                      value=int(RUN_CONFIG["workers"]),
                                      key="runs_workers")
            threads = c2.number_input("Threads", min_value=1, max_value=256,
                                      value=int(RUN_CONFIG["threads"]),
                                      key="runs_threads")
            t1, t2, t3 = st.columns(3)
            trim = t1.checkbox("Trimmomatic trimming", not RUN_CONFIG["no_trim"],
                               key="runs_trim")
            prefilter = t2.checkbox("prinseq++ prefilter", not RUN_CONFIG["no_prefilter"],
                                    key="runs_prefilter")
            dedup = t3.checkbox("Picard dedup", not RUN_CONFIG["no_dedup"],
                                key="runs_dedup")
            start = st.form_submit_button("Start mapping + extract (all refs)",
                                          disabled=bool(active or queue_active),
                                          icon=":material/play_arrow:")

        if start:
            try:
                spec = {
                    "map_cmd": jobs.build_map_cmd(
                        str(jobs.DEFAULT_MAP_SCRIPT), paths["references"],
                        paths["samples"], paths["tsv"], paths["output"], paths["log"],
                        int(workers), int(threads), trim, prefilter, dedup,
                        config_path=RUN_CONFIG_PATH),
                    "extract_script": str(jobs.DEFAULT_EXTRACT_SCRIPT),
                    "refs_dir": paths["references"],
                    "samples_dir": paths["samples"],
                    "threads": int(threads),
                }
                jobs.start_job(spec, paths["output"], run_id=selected)
                st.success("Job started. Status updates on refresh.")
                st.rerun()
            except RuntimeError as e:
                st.error(str(e))

        results = (jobs.read_job(paths["output"]) or {}).get("extract_results") or []
        if results:
            st.markdown("**Extract results:**")
            _show_table(results, alt="Extract results per reference")
        else:
            st.caption("No extract results yet.")
    else:
        st.info("No references found in this dataset.")

    st.divider()
    st.subheader("Sequential deployment")
    st.caption("Select runs to deploy one by one. Each run keeps its own config. "
               "Queue continues on job failure/cancel. Redeploying uses existing output and resume rules.")
    deployable = [r["run_id"] for r in records]
    selected_queue = st.multiselect("Runs to deploy (in selected order)", deployable,
                                    default=[], key="deploy_queue_select",
                                    format_func=_run_label)
    queue_items = []
    invalid_queue = False
    for run_id in selected_queue:
        try:
            run_refs = find_references(runs.run_paths(run_id)["references"])
        except ValueError as e:
            st.error(f"{_run_label(run_id)}: {e}")
            invalid_queue = True
            continue
        if not run_refs:
            st.error(f"{_run_label(run_id)}: no references found.")
            invalid_queue = True
            continue
        cfg = _run_deploy_form(run_id, disabled=queue_active)
        queue_items.append({
            "run_id": run_id,
            "spec": deploy_queue.build_spec(
                run_id, cfg["workers"], cfg["threads"],
                cfg["trim"], cfg["prefilter"], cfg["dedup"], config_path=RUN_CONFIG_PATH),
        })

    q1, q2 = st.columns(2)
    if q1.button("Start sequential deployment", disabled=bool(active or queue_active or invalid_queue or not queue_items),
                 key="deploy_queue_start", icon=":material/playlist_play:"):
        try:
            deploy_queue.start_queue(queue_items)
            st.success("Sequential deployment started. Status updates on refresh.")
            st.rerun()
        except Exception as e:
            st.error(f"Could not start queue: {type(e).__name__}: {e}")
    if q2.button("Cancel queue (after current job)", disabled=not queue_active,
                 key="deploy_queue_cancel", icon=":material/stop:"):
        if deploy_queue.request_cancel():
            st.warning("Queue cancel requested. Current job will finish unless you cancel it.")
            st.rerun()

    if queue:
        st.markdown("**Queue status:**")
        qmeta = {k: queue.get(k) for k in ("status", "current_index", "current_run_id", "started", "finished", "cancel_requested")}
        st.code("\n".join(f"{k}: {v}" for k, v in qmeta.items()), language="text")
        items = queue.get("items", [])
        if items:
            rows = []
            for i, it in enumerate(items, 1):
                spec = it.get("spec", {})
                rows.append({
                    "order": i,
                    "run": _run_label(it.get("run_id", "")),
                    "extract": "all references",
                    "status": it.get("status", ""),
                    "job_status": it.get("job_status", ""),
                    "error": it.get("error", ""),
                })
            _show_table(rows, alt="Sequential deployment queue")
        qtail = deploy_queue.queue_log_tail()
        if qtail:
            with st.expander("Queue log"):
                st.code(qtail, language="text")


@st.cache_resource(show_spinner=False)
def _browser_server(outdir):
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
    height = browser.browser_height(config)
    return browser.viewer_html(base, config, height=height), flags_error, height


def _extract_outputs(output_dir):
    patterns = [
        os.path.join(output_dir, jobs.EXTRACT_DIR, "*", "subset_ref.fa.fai"),
        os.path.join(output_dir, jobs.EXTRACT_DIR, "*", "*", "subset_ref.fa.fai"),
    ]
    found = []
    for pat in patterns:
        found.extend(glob.glob(pat))
    return sorted({os.path.dirname(p) for p in found})


def _browser_tab(output_dir):
    options = _extract_outputs(output_dir)
    manual = st.text_input("extract_remap output dir (overrides the choice below)", "",
                           key="br_outdir")
    if manual:
        outdir = manual
    elif options:
        outdir = st.selectbox("Extract output", options, key="br_pick")
    else:
        st.info("No extract output yet. Run a job in the Runs tab, or enter a directory above.")
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
    html, flags_error, height = _browser_view(outdir, picked, base,
                                               _browser_signature(outdir, picked))
    if flags_error:
        st.error(f"Junction read scan failed: {flags_error}")
    st.iframe(html, height=height, width="stretch", alt="Genome browser for excised reference and selected samples")
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
    tabs = st.tabs(
        ["Data", "Runs", "References", "Heatmap", "Sample", "Pair", "Browser",
         "Warnings"],
        key="main_tabs", on_change="rerun")
    (data_tab, runs_tab, refs_tab, heat_tab, sample_tab, pair_tab,
     browser_tab, warn_tab) = tabs

    if data_tab.open:
        with data_tab:
            _data_tab()
    if runs_tab.open:
        with runs_tab:
            _runs_tab()

    result_tabs = (refs_tab, heat_tab, sample_tab, pair_tab, warn_tab)
    df = pd.DataFrame(columns=loader.COLUMNS + ["bam_exists"])
    mode, status = "n/a", "unknown"
    load_error = None
    if any(tab.open for tab in result_tabs):
        signature = loader.signature(tsv_path, log_path, output_dir) + (refs_dir,)
        try:
            df, mode, status = _load(source, tsv_path, refs_dir, output_dir,
                                     log_path, signature)
        except Exception as e:
            load_error = f"{type(e).__name__}: {e}"
        c1, c2, c3 = st.columns(3)
        c1.metric("Data source", mode.upper())
        c2.metric("Run status", status)
        c3.metric("Completed combinations", len(df))
        if load_error:
            st.error(f"Could not load results: {load_error}")
        progress = loader.read_progress(log_path)
        if progress and status == "running":
            done, total = progress
            if total:
                st.progress(min(done / total, 1.0),
                            text=f"Mapping progress: {done}/{total} combinations")
        ok = df[df["error"].fillna("") == ""]
        if ok.empty and not load_error:
            st.info("No completed combinations with results yet.")
    else:
        ok = df[df["error"].fillna("") == ""]

    if refs_tab.open:
        with refs_tab:
            refs = views.references_table(df)
            _show_table(refs, alt="Mapping results by reference")
            st.plotly_chart(views.references_figure(refs),
                            alt="Mapped reads by reference")

    if heat_tab.open:
        with heat_tab:
            metric = st.selectbox("Metric", list(views.METRICS),
                                  format_func=views.METRICS.get, key="heat_metric")
            fig = views.heatmap_figure(ok, metric)
            if fig is None:
                st.info("Nothing to show yet.")
            else:
                st.plotly_chart(fig, alt="Heatmap of mapping metrics")

    if sample_tab.open:
        with sample_tab:
            if not ok.empty:
                sample = st.selectbox("Sample", sorted(ok["sample"].unique()),
                                      key="sample_pick")
                sdf = ok[ok["sample"] == sample]
                metric_s = st.selectbox("Metric", list(views.METRICS),
                                        format_func=views.METRICS.get,
                                        key="sample_metric")
                st.plotly_chart(views.sample_figure(sdf, metric_s),
                                alt=f"Mapping metrics for sample {sample}")
                _show_table(sdf[["pathogen", "assembly_id", "total_reads",
                                 "mapped_reads", "mapped_pct", "coverage",
                                 "avg_depth", "depth_covered", "duplicates_removed"]],
                            alt=f"Details for sample {sample}")

    if pair_tab.open:
        with pair_tab:
            if not ok.empty:
                sample = st.selectbox("Sample", sorted(ok["sample"].unique()),
                                      key="pair_sample")
                pairs = ok[ok["sample"] == sample].reset_index(drop=True)
                labels = views.ref_label(pairs).tolist()
                choice = st.selectbox("Reference", labels, key="pair_ref")
                row = pairs.iloc[labels.index(choice)]

                m1, m2, m3, m4 = st.columns(4)
                m1.metric("Mapped (%)", _fmt(row["mapped_pct"]))
                m2.metric("Coverage (%)", _fmt(row["coverage"]))
                m3.metric("Average depth", _fmt(row["avg_depth"]))
                m4.metric("Duplicates removed", _fmt(row["duplicates_removed"]))
                st.caption(f"BAM: {row['bam_path']}")

                bam = str(row["bam_path"])
                stats_file = (bam[:-len(".bam")] + ".stats.txt"
                              if bam.endswith(".bam") else "")
                if stats_file and os.path.exists(stats_file):
                    _show_table(views.stats_table(loader.parse_stats_sn(stats_file)),
                                alt="samtools stats summary")
                else:
                    st.warning("samtools stats file not found for this pair.")

                if st.checkbox("Show per-position coverage (binned, cached)",
                               key="pair_depth"):
                    try:
                        with st.spinner("Loading binned depth..."):
                            binned = depth.binned_depth(bam, depth.DEFAULT_BINS)
                        st.plotly_chart(views.depth_figure(
                            binned, f"{sample} vs {choice}"),
                            alt=f"Binned coverage for {sample} versus {choice}")
                    except FileNotFoundError:
                        st.warning("BAM and depth artifact not found for this pair.")
                    except Exception as e:
                        st.error(f"Depth computation failed: {e}")

    if browser_tab.open:
        with browser_tab:
            _browser_tab(output_dir)

    if warn_tab.open:
        with warn_tab:
            missing = df[~df["bam_exists"]]
            if missing.empty:
                st.success("All BAM paths exist.")
            else:
                st.warning(f"{len(missing)} row(s) have no BAM on disk (cleaned up or missing). Depth remains available when a depth artifact exists.")
                _show_table(missing[["sample", "pathogen", "assembly_id", "bam_path"]],
                            alt="Rows with missing BAM files")
            errored = df[df["error"].fillna("") != ""]
            if not errored.empty:
                st.error(f"{len(errored)} row(s) with errors.")
                _show_table(errored[["sample", "pathogen", "assembly_id", "error"]],
                            alt="Rows with mapping errors")


def _fmt(value):
    return "n/a" if value is None or value != value else f"{value:.2f}"


if live:
    st.fragment(run_every=int(interval))(render)()
else:
    render()
