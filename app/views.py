"""Plotly figures and tables for the mapping_tracker app."""
from __future__ import annotations

import pandas as pd
import plotly.graph_objects as go

STATS_KEYS = [
    "raw total sequences", "reads mapped", "reads mapped and paired",
    "reads duplicated", "reads MQ0", "average length", "average quality",
    "insert size average", "insert size standard deviation", "error rate",
    "bases mapped (cigar)",
]

METRICS = {
    "mapped_pct": "Mapped (%)",
    "coverage": "Coverage (%)",
    "avg_depth": "Average depth",
    "depth_covered": "Depth over covered bases",
}


def ref_label(df):
    return df["pathogen"].astype(str) + " / " + df["assembly_id"].astype(str)


def references_table(df):
    refs = (df[["pathogen", "assembly_id", "ref_length", "fasta_path"]]
            .drop_duplicates(subset=["pathogen", "assembly_id"])
            .sort_values(["pathogen", "assembly_id"])
            .reset_index(drop=True))
    return refs


def references_figure(refs):
    fig = go.Figure()
    for pathogen, grp in refs.groupby("pathogen", sort=True):
        fig.add_bar(x=grp["assembly_id"], y=grp["ref_length"], name=str(pathogen))
    fig.update_layout(title="Reference length per assembly", barmode="group",
                      xaxis_title="assembly", yaxis_title="length (bp)",
                      height=420, margin=dict(t=60, b=40))
    return fig


def heatmap_figure(df, metric):
    if df.empty:
        return None
    work = df.assign(ref=ref_label(df))
    pivot = work.pivot_table(index="sample", columns="ref", values=metric,
                             aggfunc="first")
    fig = go.Figure(go.Heatmap(
        z=pivot.values, x=list(pivot.columns), y=list(pivot.index),
        colorscale="Viridis", colorbar=dict(title=METRICS.get(metric, metric)),
        hovertemplate="sample=%{y}<br>reference=%{x}<br>value=%{z:.3f}<extra></extra>",
    ))
    fig.update_layout(title=METRICS.get(metric, metric),
                      height=max(300, 40 * len(pivot.index) + 150),
                      xaxis_title="reference", yaxis_title="sample",
                      margin=dict(t=60, b=40))
    return fig


def sample_figure(sample_df, metric):
    fig = go.Figure(go.Bar(x=ref_label(sample_df), y=sample_df[metric]))
    fig.update_layout(title=f"{METRICS.get(metric, metric)} per reference",
                      yaxis_title=METRICS.get(metric, metric), height=380,
                      margin=dict(t=60, b=40))
    return fig


def stats_table(sn):
    rows = [{"metric": k, "value": sn[k]} for k in STATS_KEYS if k in sn]
    return pd.DataFrame(rows, columns=["metric", "value"])


def depth_figure(depth_df, title):
    """Binned depth along the concatenated contigs, with contig boundaries."""
    fig = go.Figure()
    offset = 0
    boundaries = []
    xs, ys = [], []
    for index, (contig, grp) in enumerate(depth_df.groupby("contig", sort=False)):
        length = int(grp["end"].max())
        if index > 0:
            boundaries.append((offset, str(contig)))
        centers = offset + (grp["start"] + grp["end"]) / 2
        xs.extend(centers.tolist() + [None])
        ys.extend(grp["mean_depth"].tolist() + [None])
        offset += length
    fig.add_scatter(x=xs, y=ys, mode="lines", name="mean depth",
                    hovertemplate="pos=%{x:.0f}<br>depth=%{y:.2f}<extra></extra>")
    show_names = len(boundaries) <= 20  # avoid clutter on fragmented drafts
    for x, name in boundaries:
        fig.add_vline(x=x, line_dash="dot", line_color="grey",
                      annotation_text=name if show_names else None,
                      annotation_position="top")
    fig.update_layout(title=title, xaxis_title="position (bp, contigs concatenated)",
                      yaxis_title="mean depth", height=420, margin=dict(t=60, b=40))
    return fig
