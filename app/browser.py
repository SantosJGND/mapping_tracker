"""igv.js genome browser support for the extract_remap output.

Reads an extract_remap output directory (subset_ref.fa(.fai), coverage_intervals.bed,
<sample>/remapped.sorted.bam) and builds:
- an igv.js config (subset reference, one BAM track per selected sample,
  junction overlay track),
- a junction BED overlay (one 2 bp feature per junction between joined intervals),
- a subset -> original coordinate lookup.

Served by a small range-capable HTTP server on 127.0.0.1 (igv.js reads BAM and FASTA
through HTTP byte-range requests). The server exposes the extract output dir at /,
the vendored igv.js at /static/, and generated files at /cache/.
Nothing here writes into the extract_remap output directory.
"""
from __future__ import annotations

import hashlib
import json
import mimetypes
import os
import re
import subprocess
import threading
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from urllib.parse import unquote, urlsplit

APP_DIR = Path(__file__).resolve().parent
STATIC_DIR = APP_DIR / "static"
CACHE_DIR = APP_DIR / "cache" / "browser"
SAMTOOLS = "samtools"
MAX_TRACKS = 10
CHUNK = 1024 * 1024


@dataclass(frozen=True)
class Interval:
    subset_contig: str
    subset_start: int  # 0-based, half-open
    subset_end: int
    orig_contig: str
    orig_start: int
    orig_end: int


# ------------------------------------------------------------- extract outputs

def read_contigs(outdir: str) -> List[str]:
    """Subset contig names in subset_ref.fa.fai order."""
    with open(os.path.join(outdir, "subset_ref.fa.fai")) as fh:
        return [line.split("\t", 1)[0] for line in fh if line.strip()]


def read_intervals(outdir: str) -> List[Interval]:
    """Rows of coverage_intervals.bed ('#' header skipped)."""
    rows: List[Interval] = []
    with open(os.path.join(outdir, "coverage_intervals.bed")) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append(Interval(f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])))
    return rows


def list_samples(outdir: str) -> List[str]:
    """Sample subdirectories holding a remapped, indexed BAM."""
    root = Path(outdir)
    if not root.is_dir():
        return []
    return sorted(
        d.name for d in root.iterdir()
        if d.is_dir()
        and (d / "remapped.sorted.bam").is_file()
        and (d / "remapped.sorted.bam.bai").is_file()
    )


def subset_to_original(intervals: List[Interval], contig: str, pos: int) -> Optional[Tuple[str, int]]:
    """Map a 0-based subset position to (original contig, 0-based original position)."""
    for iv in intervals:
        if iv.subset_contig == contig and iv.subset_start <= pos < iv.subset_end:
            return iv.orig_contig, iv.orig_start + (pos - iv.subset_start)
    return None


def junction_gaps(intervals: List[Interval]) -> Dict[str, List[Tuple[int, int]]]:
    """Subset contig -> list of (start, end) 0-based half-open junction gaps.

    A junction sits between consecutive intervals on the same subset contig. With
    spacer 0 the gap is empty (start == end == next interval's start). With a spacer
    S, the gap is the N run [prev subset_end, next subset_start), of length S.
    """
    out: Dict[str, List[Tuple[int, int]]] = {}
    for a, b in zip(intervals, intervals[1:]):
        if a.subset_contig == b.subset_contig:
            out.setdefault(a.subset_contig, []).append((a.subset_end, b.subset_start))
    return out


def junction_bed_lines(intervals: List[Interval]) -> List[str]:
    """Overlay BED of junctions.

    Spacer 0: a 2 bp feature over the two bases flanking the join.
    Spacer S: the N run itself. The feature name records the original coordinates on
    either side of the join.
    """
    lines = []
    for a, b in zip(intervals, intervals[1:]):
        if a.subset_contig != b.subset_contig:
            continue
        g0, g1 = a.subset_end, b.subset_start
        if g1 > g0:
            start, end = g0, g1
        else:
            start, end = g0 - 1, g0 + 1
        name = f"{a.orig_contig}:{a.orig_end}|{b.orig_contig}:{b.orig_start}"
        lines.append(f"{a.subset_contig}\t{start}\t{end}\t{name}\n")
    return lines


def write_junctions(outdir: str, intervals: List[Interval]) -> Optional[str]:
    """Write the junction overlay into the app cache. Returns the relative URL path, or None."""
    lines = junction_bed_lines(intervals)
    if not lines:
        return None
    key = hashlib.sha1(os.path.abspath(outdir).encode()).hexdigest()[:12]
    folder = CACHE_DIR / key
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "junctions.bed").write_text("".join(lines))
    return f"/cache/{key}/junctions.bed"


# ------------------------------------------------- per-read junction flags

_CIGAR_REF_RE = re.compile(r"(\d+)([MDN=X])")


def _overlaps_junction(start0: int, end0: int, gap: Tuple[int, int]) -> bool:
    g0, g1 = gap
    if g1 > g0:  # spacer: any overlap with the N run
        return start0 < g1 and end0 > g0
    return start0 < g0 < end0  # no spacer: covers both flanking bases


def spanning_read_lines(sample: str, bam: str,
                        gaps: Dict[str, List[Tuple[int, int]]]) -> List[Tuple[str, int, str]]:
    """BED rows (contig, start, line) for primary reads whose aligned span crosses a junction.

    Spacer 0: flagged when the read covers both joined bases. Spacer S: flagged when
    the read overlaps the N run. Reads crossing several junctions are reported once.
    """
    regions = [f"{c}:{g0 + 1}-{max(g1, g0 + 1)}" for c, gs in gaps.items() for g0, g1 in gs]
    if not regions:
        return []
    out = subprocess.run([SAMTOOLS, "view", "-F", "0x904", bam, *regions],
                         capture_output=True, text=True, check=True).stdout
    seen = set()
    rows: List[Tuple[str, int, str]] = []
    for line in out.splitlines():
        f = line.split("\t", 6)
        qname, flag, contig, cigar = f[0], int(f[1]), f[2], f[5]
        start0 = int(f[3]) - 1
        end0 = start0 + sum(int(n) for n, _ in _CIGAR_REF_RE.findall(cigar))
        if not any(_overlaps_junction(start0, end0, g) for g in gaps.get(contig, ())):
            continue
        key = (qname, flag, contig, start0)
        if key in seen:
            continue
        seen.add(key)
        strand = "-" if flag & 0x10 else "+"
        rows.append((contig, start0, f"{contig}\t{start0}\t{end0}\t{sample}:{qname}\t0\t{strand}\n"))
    return rows


def write_read_flags(outdir: str, samples: List[str],
                     intervals: List[Interval]) -> Optional[str]:
    """Write a BED track of junction-spanning reads for the selected samples.

    Returns the relative URL path, or None when there are no junctions or no flagged reads.
    The file name includes a hash of the sample selection, so the browser never reuses
    a stale track after the selection changes.
    """
    gaps = junction_gaps(intervals)
    if not gaps:
        return None
    rows: List[Tuple[str, int, str]] = []
    for s in samples:
        bam = os.path.join(outdir, s, "remapped.sorted.bam")
        rows.extend(spanning_read_lines(s, bam, gaps))
    if not rows:
        return None
    order = {c: i for i, c in enumerate(gaps)}
    rows.sort(key=lambda r: (order.get(r[0], len(order)), r[1]))
    key = hashlib.sha1(os.path.abspath(outdir).encode()).hexdigest()[:12]
    sel = hashlib.sha1("|".join(sorted(samples)).encode()).hexdigest()[:8]
    folder = CACHE_DIR / key
    folder.mkdir(parents=True, exist_ok=True)
    name = f"junction_reads_{sel}.bed"
    (folder / name).write_text("".join(line for _, _, line in rows))
    return f"/cache/{key}/{name}"


# ------------------------------------------------------------- igv.js config

def build_config(base_url: str, samples: List[str], locus_contig: str,
                 junction_path: Optional[str], read_flags_path: Optional[str] = None) -> dict:
    tracks = []
    if junction_path:
        tracks.append({
            "name": "junctions (joined intervals)",
            "type": "annotation",
            "format": "bed",
            "url": base_url + junction_path,
            "displayMode": "EXPANDED",
            "color": "rgb(220,20,60)",
        })
    if read_flags_path:
        tracks.append({
            "name": "junction-spanning reads",
            "type": "annotation",
            "format": "bed",
            "url": base_url + read_flags_path,
            "displayMode": "EXPANDED",
            "color": "rgb(255,140,0)",
        })
    for s in samples[:MAX_TRACKS]:
        tracks.append({
            "name": s,
            "type": "alignment",
            "format": "bam",
            "url": f"{base_url}/{s}/remapped.sorted.bam",
            "indexURL": f"{base_url}/{s}/remapped.sorted.bam.bai",
        })
    return {
        "genome": {
            "id": "subset_ref",
            "name": "subset reference",
            "fastaURL": f"{base_url}/subset_ref.fa",
            "indexURL": f"{base_url}/subset_ref.fa.fai",
        },
        "locus": locus_contig,
        "tracks": tracks,
    }


def _script_json(obj) -> str:
    """JSON safe to embed in a <script> block (no '<', '>' or '&' can close it)."""
    return (json.dumps(obj).replace("<", "\\u003c")
            .replace(">", "\\u003e").replace("&", "\\u0026"))


def browser_height(config: dict, min_height: int = 900, per_track: int = 140, base: int = 260) -> int:
    tracks = len(config.get("tracks", [])) if isinstance(config, dict) else 0
    return max(min_height, base + per_track * max(1, tracks))


def viewer_html(base_url: str, config: dict, height: int = 900) -> str:
    return (
        "<style>html,body{margin:0;padding:0;background:#fff;width:100%;}"
        f"#igv{{width:100%;min-height:{height}px;height:{height}px;}}</style>\n"
        f'<div id="igv"></div>\n'
        f'<script src="{base_url}/static/igv.min.js"></script>\n'
        "<script>\n"
        f"igv.createBrowser(document.getElementById('igv'), {_script_json(config)})\n"
        "  .catch(function (e) { document.getElementById('igv').innerText = 'igv.js error: ' + e; });\n"
        "</script>\n"
    )


# ------------------------------------------------------------------ server

_ROUTES = (("/static/", STATIC_DIR), ("/cache/", CACHE_DIR))
_RANGE_RE = re.compile(r"bytes=(\d*)-(\d*)$")


class _RangeHandler(BaseHTTPRequestHandler):
    server_version = "mapping_tracker_browser"

    def log_message(self, fmt, *args):  # keep the Streamlit console quiet
        pass

    def _cors(self):
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Range")
        self.send_header("Access-Control-Expose-Headers", "Content-Range, Content-Length, Accept-Ranges")

    def _target(self) -> Optional[Path]:
        path = unquote(urlsplit(self.path).path)
        root = Path(self.server.outdir)
        rel = path
        for prefix, base in _ROUTES:
            if path.startswith(prefix):
                root, rel = base, path[len(prefix):]
                break
        root = root.resolve()
        target = (root / rel.lstrip("/")).resolve()
        if not str(target).startswith(str(root) + os.sep) or not target.is_file():
            return None
        return target

    def do_OPTIONS(self):
        self.send_response(204)
        self._cors()
        self.send_header("Access-Control-Allow-Methods", "GET, HEAD")
        self.end_headers()

    def do_HEAD(self):
        self._send(head_only=True)

    def do_GET(self):
        self._send(head_only=False)

    def _send(self, head_only: bool):
        target = self._target()
        if target is None:
            self.send_error(404)
            return
        size = target.stat().st_size
        ctype = mimetypes.guess_type(str(target))[0] or "application/octet-stream"
        header = self.headers.get("Range")
        if size == 0:
            self.send_response(200)
            self._cors()
            self.send_header("Accept-Ranges", "bytes")
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        start, end, status = 0, size - 1, 200
        if header:
            m = _RANGE_RE.match(header.strip())
            if not m or (m.group(1) == "" and m.group(2) == ""):
                self.send_error(416)
                return
            if m.group(1) == "":  # suffix range: last N bytes
                start = max(0, size - int(m.group(2)))
            else:
                start = int(m.group(1))
                if m.group(2):
                    end = min(int(m.group(2)), size - 1)
            if start >= size or end < start:
                self.send_response(416)
                self._cors()
                self.send_header("Content-Range", f"bytes */{size}")
                self.end_headers()
                return
            status = 206
        length = end - start + 1
        self.send_response(status)
        self._cors()
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(length))
        if status == 206:
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.end_headers()
        if head_only:
            return
        with open(target, "rb") as fh:
            fh.seek(start)
            remaining = length
            while remaining > 0:
                chunk = fh.read(min(CHUNK, remaining))
                if not chunk:
                    break
                self.wfile.write(chunk)
                remaining -= len(chunk)


def start_server(outdir: str) -> ThreadingHTTPServer:
    """Start the range server on a free localhost port, serving outdir at /. Daemon thread."""
    server = ThreadingHTTPServer(("127.0.0.1", 0), _RangeHandler)
    server.outdir = str(Path(outdir).resolve())
    threading.Thread(target=server.serve_forever, daemon=True,
                     name="browser-range-server").start()
    return server
