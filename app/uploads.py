from __future__ import annotations

import gzip
import io
import os
import re
import stat
import sys
import zipfile
from pathlib import Path
from typing import Any, Dict, List

APP_DIR = Path(__file__).resolve().parent
_PARENT = str(APP_DIR.parent)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

import runs  # noqa: E402
from map_references import find_references, find_samples  # noqa: E402

FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
GZIP_MAGIC = b"\x1f\x8b"


def _content(uploaded: Any) -> bytes:
    if hasattr(uploaded, "getvalue"):
        return bytes(uploaded.getvalue())
    if hasattr(uploaded, "read"):
        try:
            uploaded.seek(0)
        except (AttributeError, OSError):
            pass
        return bytes(uploaded.read())
    return bytes(uploaded)


def _safe_name(name: Any) -> str:
    base = os.path.basename(str(name).replace("\\", "/"))
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", base).strip("._")
    return cleaned or "upload"


def _unique_path(folder: Path, name: str) -> Path:
    target = folder / name
    if not target.exists():
        return target
    stem, suffix = os.path.splitext(name)
    i = 2
    while True:
        target = folder / f"{stem}_{i}{suffix}"
        if not target.exists():
            return target
        i += 1


def save_fastq(uploaded: Any, samples_dir: str) -> str:
    data = _content(uploaded)
    name = _safe_name(getattr(uploaded, "name", "sample.fastq"))
    lower = name.lower()
    is_gzip = data.startswith(GZIP_MAGIC)

    if lower.endswith((".fastq.gz", ".fq.gz")):
        if not is_gzip:
            data = gzip.compress(data, mtime=0)
    elif lower.endswith((".fastq", ".fq")):
        stem = name.rsplit(".", 1)[0]
        name = stem + ".fastq.gz"
        if not is_gzip:
            data = gzip.compress(data, mtime=0)
    else:
        data = gzip.compress(data, mtime=0)
        name = name + ".fastq.gz"

    target_dir = Path(samples_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    target = _unique_path(target_dir, name)
    target.write_bytes(data)
    return str(target)


def _zip_members(zf: zipfile.ZipFile) -> List[zipfile.ZipInfo]:
    members = []
    for info in zf.infolist():
        if info.is_dir():
            continue
        name = info.filename.replace("\\", "/")
        parts = [p for p in name.split("/") if p not in {"", "."}]
        if name.startswith("/") or ".." in parts:
            raise ValueError(f"Unsafe ZIP entry: {info.filename}")
        if stat.S_ISLNK(info.external_attr >> 16):
            raise ValueError(f"Symlink ZIP entry: {info.filename}")
        members.append(info)
    return members


def _write_member(zf: zipfile.ZipFile, info: zipfile.ZipInfo, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(zf.read(info))


def _extract_nested(zf: zipfile.ZipFile, members: List[zipfile.ZipInfo],
                    dest_dir: Path, label: str) -> List[str]:
    roots = []
    for info in members:
        parts = info.filename.replace("\\", "/").split("/")
        for i, part in enumerate(parts[:-1]):
            if part == "ncbi_dataset" and i + 1 < len(parts) and parts[i + 1] == "data":
                root = "/".join(parts[:i])
                if root not in roots:
                    roots.append(root)
    if not roots:
        return []

    written = []
    used_names = set()
    dest_dir.mkdir(parents=True, exist_ok=True)
    for root in roots:
        prefix = root + "/" if root else ""
        base = root.split("/")[-1] if root else label
        base = _safe_name(base)
        if base in used_names:
            i = 2
            while f"{base}_{i}" in used_names:
                i += 1
            base = f"{base}_{i}"
        used_names.add(base)
        for info in members:
            name = info.filename.replace("\\", "/")
            if not name.startswith(prefix + "ncbi_dataset/data/"):
                continue
            rel = name[len(prefix):]
            target = _unique_path(dest_dir / base, _safe_member_path(rel))
            _write_member(zf, info, target)
            written.append(str(target))
    return written


def _safe_member_path(path: str) -> str:
    parts = [_safe_name(p) for p in path.split("/") if p not in {"", "."}]
    return "/".join(parts)


def _extract_flat(zf: zipfile.ZipFile, members: List[zipfile.ZipInfo],
                  dest_dir: Path) -> List[str]:
    written = []
    dest_dir.mkdir(parents=True, exist_ok=True)
    for info in members:
        name = info.filename.replace("\\", "/")
        base = _safe_name(name.split("/")[-1])
        if not base.lower().endswith(FASTA_SUFFIXES):
            continue
        if base.lower().endswith(".fna"):
            base = base.rsplit(".", 1)[0] + ".fa"
        target = _unique_path(dest_dir, base)
        _write_member(zf, info, target)
        written.append(str(target))
    return written


def save_reference_zip(uploaded: Any, references_dir: str) -> List[str]:
    data = _content(uploaded)
    label = _safe_name(Path(str(getattr(uploaded, "name", "references.zip"))).stem)
    dest_dir = Path(references_dir)
    with zipfile.ZipFile(io.BytesIO(data)) as zf:
        members = _zip_members(zf)
        written = _extract_nested(zf, members, dest_dir, label)
        if not written:
            written = _extract_flat(zf, members, dest_dir)
    if not written:
        raise ValueError(f"No reference FASTA found in {getattr(uploaded, 'name', 'ZIP')}")
    return written


def create_run(name: str, reference_uploads: List[Any],
               fastq_uploads: List[Any]) -> Dict[str, Any]:
    run = runs.create_run(name)
    paths = runs.run_paths(run["run_id"])
    try:
        for uploaded in reference_uploads or []:
            save_reference_zip(uploaded, paths["references"])
        for uploaded in fastq_uploads or []:
            save_fastq(uploaded, paths["samples"])

        references = find_references(paths["references"])
        samples = find_samples(paths["samples"])
        if not references:
            raise ValueError("No reference genomes found in uploaded ZIPs")
        if not samples:
            raise ValueError("No FASTQ samples found in uploaded files")

        return runs.update_run(
            run["run_id"],
            reference_count=len(references),
            sample_count=len(samples),
            references=[{"pathogen": r["pathogen"], "assembly_id": r["assembly_id"]}
                        for r in references],
            samples=[s["name"] for s in samples],
        )
    except Exception as e:
        runs.update_run(run["run_id"], status="invalid", error=f"{type(e).__name__}: {e}")
        raise
