import gzip
import io
import os
import shutil
import sys
import tempfile
import unittest
import zipfile
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "app"
sys.path.insert(0, str(APP))
import runs  # noqa: E402
import uploads  # noqa: E402


class Upload:
    def __init__(self, name, data):
        self.name = name
        self._data = data

    def getvalue(self):
        return self._data


def zip_bytes(entries):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for name, data in entries:
            zf.writestr(name, data)
    return buf.getvalue()


class UploadTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="uploads_test_")
        self.old_data_dir = runs.DATA_DIR
        runs.DATA_DIR = os.path.join(self.tmp, "data")

    def tearDown(self):
        runs.DATA_DIR = self.old_data_dir
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_flat_reference_zip(self):
        refs = os.path.join(self.tmp, "refs")
        written = uploads.save_reference_zip(
            Upload("panel.zip", zip_bytes([("refA.fna", b">chr1\nACGT\n")])), refs)
        self.assertEqual([Path(p).name for p in written], ["refA.fa"])
        found = uploads.find_references(refs)
        self.assertEqual([r["assembly_id"] for r in found], ["refA"])

    def test_ncbi_reference_zip(self):
        refs = os.path.join(self.tmp, "refs")
        uploads.save_reference_zip(Upload("panel.zip", zip_bytes([
            ("PathoA/ncbi_dataset/data/GCF_1/GCF_1_genomic.fna", b">c\nAAAA\n"),
        ])), refs)
        found = uploads.find_references(refs)
        self.assertEqual([(r["pathogen"], r["assembly_id"]) for r in found],
                         [("PathoA", "GCF_1")])

    def test_unsafe_zip_entry_rejected(self):
        with self.assertRaises(ValueError):
            uploads.save_reference_zip(
                Upload("bad.zip", zip_bytes([("../escape.fa", b">c\nA\n")])),
                os.path.join(self.tmp, "refs"))

    def test_plain_fastq_is_gzipped(self):
        samples = os.path.join(self.tmp, "samples")
        target = uploads.save_fastq(Upload("s_R1.fastq", b"@r\nACGT\n+\nIIII\n"), samples)
        self.assertTrue(target.endswith("s_R1.fastq.gz"))
        with gzip.open(target, "rb") as fh:
            self.assertEqual(fh.read(), b"@r\nACGT\n+\nIIII\n")

    def test_gzipped_fastq_is_kept(self):
        samples = os.path.join(self.tmp, "samples")
        data = gzip.compress(b"@r\nACGT\n+\nIIII\n", mtime=0)
        target = uploads.save_fastq(Upload("s_R2.fastq.gz", data), samples)
        self.assertEqual(Path(target).read_bytes(), data)

    def test_create_run_detects_uploads(self):
        run = uploads.create_run(
            "test run",
            [Upload("panel.zip", zip_bytes([("refA.fa", b">c\nACGT\n")]))],
            [Upload("s_R1.fastq", b"@r\nACGT\n+\nIIII\n")])
        paths = runs.run_paths(run["run_id"])
        self.assertEqual(run["status"], "ready")
        self.assertEqual(run["reference_count"], 1)
        self.assertEqual(run["sample_count"], 1)
        self.assertTrue(Path(paths["references"], "refA.fa").is_file())
        self.assertTrue(Path(paths["samples"], "s_R1.fastq.gz").is_file())


if __name__ == "__main__":
    unittest.main()
