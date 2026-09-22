"""Biological controls and an independent oracle from minimap2's real index."""
from collections import Counter
import gzip
import json
import math
from pathlib import Path
import random
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest

from reference_hardness import minimap2_minimizers, quantify_reference_hardness, reference_hardness_score


def random_dna(length, seed=7):
    return "".join(random.Random(seed).choices("ACGT", k=length))


def read_mmi_seeds(path):
    """Decode keys/positions from a native index, only for the parity oracle.

    Layout: minimap2 v2.31 index.c mm_idx_dump / worker_post. The production
    scorer does not depend on this internal binary format.
    """
    with path.open("rb") as handle:
        def read(fmt):
            return struct.unpack("=" + fmt, handle.read(struct.calcsize("=" + fmt)))

        assert handle.read(4) == b"MMI\x02"
        w, k, bits, n_seq, flags = read("5I")
        for _ in range(n_seq):
            size, = read("B")
            handle.read(size)
            read("I")
        seeds = Counter()
        for bucket in range(1 << bits):
            n, = read("I")
            positions = read(f"{n}Q")
            size, = read("I")
            for _ in range(size):
                encoded_key, value = read("2Q")
                key = ((encoded_key >> 1) << bits) | bucket
                if encoded_key & 1:
                    occurrences = (value,)
                else:
                    offset, count = value >> 32, value & 0xFFFFFFFF
                    occurrences = positions[offset:offset + count]
                for position in occurrences:
                    rid = position >> 32
                    end = ((position & 0xFFFFFFFF) >> 1) + 1
                    strand = position & 1
                    seeds[rid, key, end - k, end, strand] += 1
        assert not handle.read(1), "oracle expects a single --idx-no-seq index"
        return seeds


class ReferenceHardnessTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "reference.fa"

    def score(self, records, **kwargs):
        self.path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))
        return quantify_reference_hardness(self.path, **kwargs)

    def test_random_duplicated_and_low_complexity_controls(self):
        sequence = random_dna(10000)
        unique = self.score([("a", sequence)])["summary"]
        duplicated = self.score([("a", sequence), ("b", sequence)])["summary"]
        low_complexity = self.score([("a", "A" * 20000)])["summary"]
        self.assertLess(unique["hardness_score"], 1)
        self.assertGreater(unique["unique_minimizers_per_kb"], 100)
        self.assertEqual(duplicated["repetitive_minimizer_fraction"], 1)
        self.assertEqual(duplicated["longest_unique_anchor_gap_bp"], 10000)
        self.assertGreater(duplicated["seed_hardness_score"], 70)
        self.assertGreater(duplicated["hardness_score"], 20)
        self.assertGreater(low_complexity["hardness_score"], 95)
        self.assertGreater(low_complexity["mean_seed_copy_number"], 10000)

    def test_reverse_complement_and_within_contig_repeats(self):
        sequence = random_dna(12000)
        rc = sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        reverse = self.score([("forward", sequence), ("reverse", rc)])["summary"]
        same_contig = self.score([("a", sequence + "N" * 100 + sequence)])["summary"]
        self.assertGreater(reverse["repetitive_minimizer_fraction"], 0.999)
        self.assertGreater(same_contig["repetitive_minimizer_fraction"], 0.999)
        self.assertGreater(reverse["seed_hardness_score"], 70)
        self.assertGreater(reverse["hardness_score"], 20)

    def test_long_repeats_harder_than_fragmented_repeats(self):
        repeat = random_dna(6000, 13)
        flank = random_dna(10000, 27)
        long = flank[:5000] + repeat + flank[5000:] + repeat
        fragmented = "".join(flank[i * 500:(i + 1) * 500] + repeat[i * 300:(i + 1) * 300] * 2
                             for i in range(20))
        a = self.score([("a", long)])["summary"]
        b = self.score([("a", fragmented)])["summary"]
        self.assertGreater(a["longest_unique_anchor_gap_bp"], 5900)
        self.assertLess(b["longest_unique_anchor_gap_bp"], 610)
        self.assertGreater(a["anchor_gap_burden"], b["anchor_gap_burden"] + 0.15)
        self.assertGreater(a["hardness_score"], b["hardness_score"] + 5)

    def test_ns_do_not_become_repeats_or_bridge_gaps(self):
        report = self.score([("a", "A" * 2000 + "N" * 30000 + "A" * 2000)])
        summary = report["summary"]
        self.assertEqual(summary["callable_bases"], 4000)
        self.assertEqual(summary["ambiguous_bases"], 30000)
        self.assertEqual(summary["longest_unique_anchor_gap_bp"], 2000)
        self.assertEqual([(g["start"], g["end"]) for g in report["long_unique_anchor_gaps"]],
                         [(0, 2000), (32000, 34000)])
        self.assertEqual([(g["start"], g["end"]) for g in report["low_complexity_tracts"]],
                         [(0, 2000), (32000, 34000)])
        self.assertEqual(summary["longest_low_complexity_tract_bp"], 2000)
        self.assertTrue(report["warnings"])

    def test_no_data_is_not_an_easy_score(self):
        n = self.score([("a", "N" * 50)])["summary"]
        short = self.score([("a", "ACGT")])["summary"]
        self.assertIsNone(n["hardness_score"])
        self.assertEqual(n["status"], "no_callable_bases")
        self.assertEqual(n["longest_unique_anchor_gap_bp"], 0)
        self.assertIsNone(short["hardness_score"])
        self.assertEqual(short["status"], "no_seeds")
        self.assertEqual(short["seedless_acgt_run_fraction"], 1)

    def test_pooling_uses_occurrences_and_bases(self):
        report = self.score([("a", random_dna(10000)), ("b", "A" * 500), ("n", "N" * 20)])
        summary = report["summary"]
        contigs = report["contigs"][:2]
        r = sum(c["minimizer_occurrences"] * c["seed_ambiguity"] for c in contigs)
        r /= summary["minimizer_occurrences"]
        g = sum(c["callable_bases"] * c["anchor_gap_burden"] for c in contigs)
        g /= summary["callable_bases"]
        self.assertAlmostEqual(summary["seed_hardness_score"], 50 * (r + g))
        mean = sum(c["callable_bases"] * c["low_complexity_mean_risk"] for c in contigs)
        mean /= summary["callable_bases"]
        peak = max(c["low_complexity_peak_risk"] for c in contigs)
        self.assertAlmostEqual(summary["low_complexity_mean_risk"], mean)
        self.assertAlmostEqual(summary["low_complexity_peak_risk"], peak)
        self.assertAlmostEqual(summary["hardness_score"], 15 * (r + g) + 70 * (0.25 * mean + 0.75 * peak))
        self.assertEqual(summary["longest_low_complexity_tract_bp"], 500)

    def test_tunable_length_scale_and_weight(self):
        small_scale = self.score([("a", "A" * 2000)], length_scale=1000)["summary"]
        large_scale = self.score([("a", "A" * 2000)], length_scale=10000)["summary"]
        seed_only = self.score([("a", "A" * 2000)], seed_weight=1, local_weight=0)["summary"]
        self.assertAlmostEqual(small_scale["anchor_gap_burden"], 2 / 3)
        self.assertAlmostEqual(large_scale["anchor_gap_burden"], 1 / 6)
        self.assertAlmostEqual(seed_only["hardness_score"], 100 * seed_only["seed_ambiguity"])

    def test_extent_still_distinguishes_repeats_above_length_scale(self):
        short = self.score([("a", "A" * 2000)])["summary"]
        long = self.score([("a", "A" * 20000)])["summary"]
        self.assertGreater(long["anchor_gap_burden"], short["anchor_gap_burden"] + 0.25)
        self.assertGreater(long["hardness_score"], short["hardness_score"] + 10)

    def test_wrapping_case_gzip_and_u(self):
        sequence = random_dna(3000)
        expected = self.score([("a", sequence)])
        compressed = self.path.with_suffix(".fa.gz")
        with gzip.open(compressed, "wt") as handle:
            handle.write(">a description\n" + "\n".join(sequence[i:i + 73].lower().replace("t", "u")
                                                       for i in range(0, len(sequence), 73)) + "\n")
        actual = quantify_reference_hardness(compressed)
        self.assertEqual(actual["summary"], expected["summary"])

    def test_invalid_inputs(self):
        for content in ("", "ACGT\n", ">\nACGT\n", ">a\n", ">a\nACGT\n>a\nTT\n", ">a\nACT!\n"):
            with self.subTest(content=content):
                self.path.write_text(content)
                with self.assertRaises(ValueError):
                    quantify_reference_hardness(self.path)
        for options in ({"k": 0}, {"k": 29}, {"w": 256}, {"w": True},
                        {"length_scale": 0}, {"long_gap": -1}, {"high_copy": 0},
                        {"seed_weight": float("nan")}, {"seed_weight": 1.1},
                        {"local_weight": -0.1}, {"local_weight": True},
                        {"hotspot_weight": float("nan")}, {"hotspot_weight": 1.1},
                        {"complexity_window": 31}, {"complexity_window": 200.5},
                        {"entropy_cutoff": 0}, {"entropy_cutoff": 7}, {"entropy_cutoff": float("nan")}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.score([("a", "A" * 50)], **options)

    def test_cli_json_and_bed(self):
        self.score([("a", "A" * 2000)])
        script = Path(__file__).resolve().parents[1] / "reference_hardness.py"
        output = self.path.with_suffix(".json")
        bed = self.path.with_suffix(".bed")
        subprocess.run([sys.executable, str(script), str(self.path), "-o", str(output),
                        "--bed", str(bed)], check=True, capture_output=True, text=True)
        self.assertGreater(json.loads(output.read_text())["summary"]["hardness_score"], 70)
        self.assertGreater(json.loads(output.read_text())["summary"]["seed_hardness_score"], 80)
        self.assertEqual(bed.read_text(), "a\t0\t2000\tunique_anchor_gap\n")
        original = self.path.read_text()
        result = subprocess.run([sys.executable, str(script), str(self.path), "-o", str(self.path)],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.path.read_text(), original)

    def test_default_cli_prints_only_one_overall_float(self):
        expected = self.score([("a", random_dna(10000)), ("b", "A" * 3000)])["summary"]["hardness_score"]
        script = Path(__file__).resolve().parents[1] / "reference_hardness.py"
        result = subprocess.run([sys.executable, str(script), str(self.path)],
                                check=True, capture_output=True, text=True)
        self.assertEqual(result.stderr, "")
        self.assertEqual(result.stdout, f"{expected}\n")
        self.assertEqual(float(result.stdout), expected)
        self.assertEqual(reference_hardness_score(self.path), expected)
        self.assertIsInstance(reference_hardness_score(self.path), float)

    def test_unscorable_reference_fails_without_printing_a_number(self):
        self.score([("a", "NNNN")])
        with self.assertRaisesRegex(ValueError, "no_callable_bases"):
            reference_hardness_score(self.path)
        script = Path(__file__).resolve().parents[1] / "reference_hardness.py"
        result = subprocess.run([sys.executable, str(script), str(self.path)],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertIn("no_callable_bases", result.stderr)
        detailed = subprocess.run([sys.executable, str(script), str(self.path), "--format", "json"],
                                  check=True, capture_output=True, text=True)
        self.assertIsNone(json.loads(detailed.stdout)["summary"]["hardness_score"])

    def test_output_cannot_overwrite_input_through_a_hardlink(self):
        self.score([("a", "A" * 1000)])
        alias = self.path.with_name("alias.txt")
        alias.hardlink_to(self.path)
        original = self.path.read_text()
        script = Path(__file__).resolve().parents[1] / "reference_hardness.py"
        result = subprocess.run([sys.executable, str(script), str(self.path), "-o", str(alias)],
                                capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(self.path.read_text(), original)

    def test_interrupted_low_complexity_can_outrank_more_repeated_complex_sequence(self):
        rng = random.Random(123)
        repeated_complex = "".join(rng.choices("ACGT", k=80)) * 100
        interrupted_ct = "".join(rng.choices("CT", weights=[0.75, 0.25], k=8000))
        complex_stats = self.score([("r", repeated_complex)])["summary"]
        simple_stats = self.score([("r", interrupted_ct)])["summary"]
        self.assertGreater(complex_stats["seed_ambiguity"], simple_stats["seed_ambiguity"])
        self.assertEqual(complex_stats["local_complexity_risk"], 0)
        self.assertGreater(simple_stats["local_complexity_risk"], 0.5)
        self.assertGreater(simple_stats["hardness_score"], complex_stats["hardness_score"] + 20)

    def test_complexity_severity_and_length_both_matter(self):
        rng = random.Random(31)
        interrupted = "".join(rng.choices("CT", k=4000))
        short = self.score([("r", "A" * 500)])["summary"]
        long = self.score([("r", "A" * 4000)])["summary"]
        varied = self.score([("r", interrupted)])["summary"]
        self.assertGreater(long["low_complexity_peak_risk"], short["low_complexity_peak_risk"])
        self.assertGreater(long["low_complexity_peak_risk"], varied["low_complexity_peak_risk"])
        self.assertAlmostEqual(long["low_complexity_peak_risk"], 4000 / 5000)

    def test_unique_flanks_do_not_dilute_the_hardest_tract(self):
        # Same focal repeat; padding comes from one nested, deterministic random sequence.
        flank = random_dna(30000, 39)
        focal = flank[-1000:] + "A" * 3000 + flank[:1000]
        padded = flank + "A" * 3000 + flank
        small = self.score([("r", focal)])["summary"]
        large = self.score([("r", padded)])["summary"]
        self.assertAlmostEqual(large["low_complexity_peak_risk"], small["low_complexity_peak_risk"])
        self.assertLess(large["low_complexity_mean_risk"], small["low_complexity_mean_risk"] / 5)
        self.assertGreater(large["local_complexity_risk"], 0.75 * small["local_complexity_risk"])

    def test_complexity_windows_match_direct_entropy_counts_and_reverse_complement(self):
        seq = random_dna(350, 15) + "CCTC" * 200 + random_dna(350, 29)
        report = self.score([("r", seq)])
        flags = []
        for pos in range(len(seq) - 199):
            counts = Counter(seq[i:i + 3] for i in range(pos, pos + 198))
            entropy = -math.fsum(n / 198 * math.log2(n / 198) for n in counts.values())
            if entropy < 4:
                flags.append((pos, (4 - entropy) / 4))
        coverage = set()
        for pos, _severity in flags:
            coverage.update(range(pos, pos + 200))
        self.assertEqual(report["summary"]["low_complexity_base_fraction"], len(coverage) / len(seq))
        for tract in report["low_complexity_tracts"]:
            severities = [severity for pos, severity in flags if tract["start"] <= pos and pos + 200 <= tract["end"]]
            severity = math.fsum(severities) / len(severities)
            effective = severity * tract["length_bp"]
            self.assertAlmostEqual(tract["mean_severity"], severity)
            self.assertAlmostEqual(tract["risk"], effective / (effective + 1000))
        rc = seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        reverse = self.score([("r", rc)])
        self.assertAlmostEqual(report["summary"]["local_complexity_risk"], reverse["summary"]["local_complexity_risk"])
        expected_intervals = [(len(seq) - t["end"], len(seq) - t["start"])
                              for t in reversed(report["low_complexity_tracts"])]
        self.assertEqual([(t["start"], t["end"]) for t in reverse["low_complexity_tracts"]], expected_intervals)

    def test_short_complexity_runs_and_contig_boundaries(self):
        report = self.score([("short", "A" * 31), ("assessed", "A" * 100), ("other", "A" * 100)])
        self.assertEqual(report["summary"]["complexity_evaluated_bases"], 200)
        self.assertEqual(report["summary"]["longest_low_complexity_tract_bp"], 100)
        self.assertEqual([(t["contig"], t["start"], t["end"], t["window_bp"])
                          for t in report["low_complexity_tracts"]],
                         [("assessed", 0, 100, 100), ("other", 0, 100, 100)])

    def test_legacy_mode_retains_original_numeric_scale(self):
        sequence = random_dna(10000)
        current = self.score([("a", sequence), ("b", sequence)])
        legacy = self.score([("a", sequence), ("b", sequence)], local_weight=0)
        self.assertEqual(current["summary"]["seed_hardness_score"], legacy["summary"]["hardness_score"])
        self.assertEqual(current["schema_version"], 2)
        self.assertEqual(current["score_model"], "local-complexity-v2")
        script = Path(__file__).resolve().parents[1] / "reference_hardness.py"
        result = subprocess.run([sys.executable, str(script), str(self.path), "--local-weight", "0"],
                                check=True, capture_output=True, text=True)
        self.assertEqual(float(result.stdout), legacy["summary"]["hardness_score"])

    @unittest.skipUnless(shutil.which("minimap2"), "minimap2 executable not installed")
    def test_exact_sketch_parity_with_native_minimap2(self):
        rng = random.Random(51)
        sequences = ["A" * 200, "ACGT" * 100, "AT" * 150, "N" * 100,
                     "A" * 100 + "N" + "C" * 100, "ACGTUacgtu" * 50,
                     random_dna(10000), random_dna(10000)[::-1]]
        sequences += ["".join(rng.choices("ACGT", k=n)) for n in range(1, 70)]
        sequences += ["".join(rng.choices("ACGTNNRY", k=500)) for _ in range(30)]
        self.path.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(sequences)))
        index = self.path.with_suffix(".mmi")
        for k, w in ((15, 10), (19, 19), (7, 1), (8, 10), (1, 3), (28, 255)):
            with self.subTest(k=k, w=w):
                subprocess.run(["minimap2", "-k", str(k), "-w", str(w), "--idx-no-seq",
                                "-t", "1", "-d", str(index), str(self.path)],
                               check=True, capture_output=True)
                native = read_mmi_seeds(index)
                python = Counter((rid, seed.key, seed.start, seed.end, seed.strand)
                                 for rid, sequence in enumerate(sequences)
                                 for seed in minimap2_minimizers(sequence, k, w))
                self.assertEqual(python, native)


if __name__ == "__main__":
    unittest.main()
