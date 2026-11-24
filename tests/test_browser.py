"""Tests for app/browser.py. Run from the mapping_tracker directory:

    python -m unittest discover -s tests -v

Needs samtools on PATH. Builds a small extract_remap-style output in a temp dir:
one subset contig with one junction, and a sample BAM with reads that do and do not
cross the join.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import urllib.error
import urllib.request
from pathlib import Path

APP = Path(__file__).resolve().parent.parent / "app"
sys.path.insert(0, str(APP))
import browser  # noqa: E402

CONTIG = "NZ_TEST.1"
SUBSET_LEN = 300
JUNCTION = 150  # 0-based subset boundary: bases 149 and 150 are joined
SEQ = "ACGT" * 25  # 100 bp


def _write_sam(path, reads, ref_len):
    lines = ["@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:{CONTIG}\tLN:{ref_len}"]
    lines += ["\t".join(r) for r in reads]
    sam = path + ".sam"
    Path(sam).write_text("\n".join(lines) + "\n")
    subprocess.run(["samtools", "view", "-b", "-o", path + ".unsorted.bam", sam], check=True)
    subprocess.run(["samtools", "sort", "-o", path, path + ".unsorted.bam"], check=True)
    subprocess.run(["samtools", "index", path], check=True)


def _read(name, flag, pos1, cigar):
    # pos1 is 1-based SAM POS
    return [name, str(flag), CONTIG, str(pos1), "60", cigar, "*", "0", "0", SEQ, "*"]


class BrowserFixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not shutil.which("samtools"):
            raise unittest.SkipTest("samtools not on PATH")
        cls.tmp = tempfile.mkdtemp(prefix="browser_test_")
        out = cls.tmp
        Path(out, "subset_ref.fa.fai").write_text(f"{CONTIG}\t{SUBSET_LEN}\t0\t60\t61\n")
        # Two intervals joined at subset 150: original 1000-1149 and 5000-5149 (1-based)
        Path(out, "coverage_intervals.bed").write_text(
            "#subset_contig\tsubset_start\tsubset_end\torig_contig\torig_start\torig_end\n"
            f"{CONTIG}\t0\t150\t{CONTIG}\t999\t1149\n"
            f"{CONTIG}\t150\t300\t{CONTIG}\t4999\t5149\n")
        sdir = Path(out, "sampleA")
        sdir.mkdir()
        _write_sam(str(sdir / "remapped.sorted.bam"), [
            _read("span_fwd", 0, JUNCTION - 49, "100M"),     # covers 99..149 -> crosses 150
            _read("span_rev", 16, JUNCTION - 49, "100M"),    # reverse strand, crosses
            _read("span_del", 0, JUNCTION - 49, "50M10D50M"),  # deletion over join
            _read("inside", 0, JUNCTION - 149, "100M"),      # ends at 50, no cross
            _read("starts_at_join", 0, JUNCTION + 1, "100M"),  # starts after join
            _read("secondary", 256, JUNCTION - 49, "100M"),  # excluded (secondary)
        ], SUBSET_LEN)
        cls.outdir = out

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)


class IntervalTests(BrowserFixture):
    def test_read_intervals_and_contigs(self):
        ivs = browser.read_intervals(self.outdir)
        self.assertEqual(len(ivs), 2)
        self.assertEqual(browser.read_contigs(self.outdir), [CONTIG])

    def test_junction_gaps_no_spacer(self):
        ivs = browser.read_intervals(self.outdir)
        self.assertEqual(browser.junction_gaps(ivs), {CONTIG: [(JUNCTION, JUNCTION)]})

    def test_junction_bed_no_spacer(self):
        ivs = browser.read_intervals(self.outdir)
        lines = browser.junction_bed_lines(ivs)
        self.assertEqual(len(lines), 1)
        self.assertEqual(lines[0].split("\t")[1:3], [str(JUNCTION - 1), str(JUNCTION + 1)])

    def test_spacer_gap_and_bed(self):
        # extract_remap --spacer 10: next interval starts 10 bases after the previous end
        ivs = [browser.Interval(CONTIG, 0, 150, CONTIG, 999, 1149),
               browser.Interval(CONTIG, 160, 310, CONTIG, 4999, 5149)]
        self.assertEqual(browser.junction_gaps(ivs), {CONTIG: [(150, 160)]})
        self.assertEqual(browser.junction_bed_lines(ivs)[0].split("\t")[1:3], ["150", "160"])

    def test_spacer_flags_reads_over_n_run(self):
        # read at 0-based 100..200 overlaps the N run [150,160); read at 0..100 does not
        rows = browser.spanning_read_lines("sampleA",
                                           os.path.join(self.outdir, "sampleA", "remapped.sorted.bam"),
                                           {CONTIG: [(150, 160)]})
        names = {r[2].split("\t")[3] for r in rows}
        self.assertIn("sampleA:span_fwd", names)
        self.assertNotIn("sampleA:inside", names)

    def test_subset_to_original(self):
        ivs = browser.read_intervals(self.outdir)
        self.assertEqual(browser.subset_to_original(ivs, CONTIG, 0), (CONTIG, 999))
        self.assertEqual(browser.subset_to_original(ivs, CONTIG, 150), (CONTIG, 4999))
        self.assertEqual(browser.subset_to_original(ivs, CONTIG, 299), (CONTIG, 5148))
        self.assertIsNone(browser.subset_to_original(ivs, CONTIG, 300))
        self.assertIsNone(browser.subset_to_original(ivs, "other", 0))


class SpanningReadTests(BrowserFixture):
    def test_flags_only_reads_crossing_join(self):
        ivs = browser.read_intervals(self.outdir)
        gaps = browser.junction_gaps(ivs)
        rows = browser.spanning_read_lines("sampleA",
                                           os.path.join(self.outdir, "sampleA", "remapped.sorted.bam"),
                                           gaps)
        names = sorted(r[2].split("\t")[3] for r in rows)
        self.assertEqual(names, ["sampleA:span_del", "sampleA:span_fwd", "sampleA:span_rev"])

    def test_strand_column(self):
        ivs = browser.read_intervals(self.outdir)
        rows = browser.spanning_read_lines("sampleA",
                                           os.path.join(self.outdir, "sampleA", "remapped.sorted.bam"),
                                           browser.junction_gaps(ivs))
        strands = {r[2].split("\t")[3]: r[2].split("\t")[5].strip() for r in rows}
        self.assertEqual(strands["sampleA:span_rev"], "-")
        self.assertEqual(strands["sampleA:span_fwd"], "+")

    def test_write_read_flags_writes_track(self):
        ivs = browser.read_intervals(self.outdir)
        url = browser.write_read_flags(self.outdir, ["sampleA"], ivs)
        self.assertIsNotNone(url)
        path = browser.CACHE_DIR / url.split("/")[2] / url.split("/")[3]
        try:
            self.assertEqual(len(path.read_text().splitlines()), 3)
        finally:
            path.unlink(missing_ok=True)

    def test_no_junctions_no_track(self):
        self.assertIsNone(browser.write_read_flags(self.outdir, ["sampleA"], [
            browser.Interval(CONTIG, 0, 300, CONTIG, 999, 1298)]))


class ConfigTests(unittest.TestCase):
    def test_script_json_escapes_markup(self):
        out = browser._script_json({"name": "</script><b>&"})
        self.assertNotIn("</script>", out)
        self.assertNotIn("<", out)
        self.assertEqual(json.loads(out)["name"], "</script><b>&")

    def test_config_caps_tracks(self):
        samples = [f"s{i}" for i in range(15)]
        cfg = browser.build_config("http://x", samples, CONTIG, None, None)
        self.assertEqual(len(cfg["tracks"]), browser.MAX_TRACKS)

    def test_config_includes_flag_track(self):
        cfg = browser.build_config("http://x", ["s"], CONTIG, "/cache/a/j.bed", "/cache/a/r.bed")
        names = [t["name"] for t in cfg["tracks"]]
        self.assertIn("junction-spanning reads", names)


class ServerTests(BrowserFixture):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.srv = browser.start_server(cls.outdir)
        cls.base = f"http://127.0.0.1:{cls.srv.server_address[1]}"

    @classmethod
    def tearDownClass(cls):
        cls.srv.shutdown()
        cls.srv.server_close()
        super().tearDownClass()

    def _get(self, path, headers=None):
        req = urllib.request.Request(self.base + path, headers=headers or {})
        try:
            with urllib.request.urlopen(req) as r:
                return r.status, r.read()
        except urllib.error.HTTPError as e:
            return e.code, b""

    def test_full_and_range(self):
        full = Path(self.outdir, "sampleA", "remapped.sorted.bam").read_bytes()
        status, body = self._get("/sampleA/remapped.sorted.bam")
        self.assertEqual((status, body), (200, full))
        status, body = self._get("/sampleA/remapped.sorted.bam", {"Range": "bytes=0-9"})
        self.assertEqual((status, body), (206, full[:10]))

    def test_traversal_and_missing_blocked(self):
        self.assertEqual(self._get("/../../etc/passwd")[0], 404)
        self.assertEqual(self._get("/nope.bam")[0], 404)


if __name__ == "__main__":
    unittest.main()
