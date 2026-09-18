#!/usr/bin/env python3
from __future__ import annotations

import os
import random
import subprocess
import tempfile
import unittest


ROOT = os.path.dirname(os.path.abspath(__file__))
BINARY = os.path.join(ROOT, "kmerpartitionmerger")


def bed_row(contig: str, start: int, end: int, target: str) -> str:
    return "\t".join(map(str, (
        contig, start, end, target, 1000, "+", "Ref_h1",
        0, end - start, end - start, end - start, 1000.0,
        100.0, 1, 0, 0,
    ))) + "\n"


def write_fasta(path: str, records) -> None:
    with open(path, "w") as out:
        for name, sequence in records:
            out.write(f">{name}\n{sequence}\n")


class KmerPartitionMergerTests(unittest.TestCase):
    def test_fully_softmasked_targets_contribute_no_kmers(self):
        if not os.path.isfile(BINARY):
            self.skipTest("run make before the test")
        with tempfile.TemporaryDirectory() as directory:
            rng = random.Random(31)
            softmasked = "".join(rng.choice("acgt") for _ in range(1500))
            fasta = os.path.join(directory, "reference.fa")
            write_fasta(fasta, [("lift", "N" * 1000)])
            first = "g1_S1_h1_contigA_0_1500"
            second = "g2_S2_h1_contigB_0_1500"
            targets = os.path.join(directory, "targets.fa")
            write_fasta(targets, [(first, softmasked), (second, softmasked)])
            bed = os.path.join(directory, "filtered.bed")
            with open(bed, "w") as out:
                out.write(bed_row("lift", 0, 100, first))
                out.write(bed_row("lift", 500, 600, second))
            output = os.path.join(directory, "merged.bed")
            groups = os.path.join(directory, "groups.tsv")
            subprocess.run([
                BINARY, "-f", fasta, "-T", targets, "-b", bed, "-o", output,
                "-g", groups, "--only-kmer", "-t", "2",
            ], check=True)

            with open(groups) as handle:
                rows = dict(
                    line.rstrip("\n").split("\t")
                    for line in handle
                    if line.strip() and not line.startswith("partition\t")
                )
            self.assertEqual(rows[first], first)
            self.assertEqual(rows[second], second)

    def test_kmers_union_lifted_and_column_four_source_regions(self):
        if not os.path.isfile(BINARY):
            self.skipTest("run make before the test")
        with tempfile.TemporaryDirectory() as directory:
            rng = random.Random(19)
            shared_unmasked = "".join(rng.choice("ACGT") for _ in range(1500))
            fasta = os.path.join(directory, "reference.fa")
            # Lifted regions contain only one distinct k-mer and do not overlap.
            write_fasta(fasta, [("lift_contig", "A" * 1000)])
            first = "g1_S1_h1_source_contig_A_0_1500"
            second = "g2_S2_h2_source_contig_B_0_1500"
            unlifted = "g3_S3_h1_source_contig_C_0_1500"
            targets = os.path.join(directory, "targets.fa")
            # All three column-4 records carry the same >1000 distinct unmasked
            # canonical 31-mers. The third has no BED row at all.
            write_fasta(targets, [
                (first, shared_unmasked),
                (second, shared_unmasked),
                (unlifted, shared_unmasked),
            ])
            bed = os.path.join(directory, "filtered.bed")
            with open(bed, "w") as out:
                out.write(bed_row("lift_contig", 0, 200, first))
                out.write(bed_row("lift_contig", 500, 700, second))
            output = os.path.join(directory, "merged.bed")
            groups = os.path.join(directory, "groups.tsv")
            subprocess.run([
                BINARY, "-f", fasta, "-T", targets, "-b", bed, "-o", output,
                "-g", groups, "--only-kmer", "-t", "2",
            ], check=True)

            with open(groups) as handle:
                rows = dict(
                    line.rstrip("\n").split("\t")
                    for line in handle
                    if line.strip() and not line.startswith("partition\t")
                )
            self.assertEqual(rows[first], "merged_g1g2g3")
            self.assertEqual(rows[second], "merged_g1g2g3")
            self.assertEqual(rows[unlifted], "merged_g1g2g3")

    def test_region_then_unmasked_kmer_merges(self):
        if not os.path.isfile(BINARY):
            self.skipTest("run make before the test")
        with tempfile.TemporaryDirectory() as directory:
            rng = random.Random(7)
            shared = "".join(rng.choice("ACGT") for _ in range(1300))
            fasta = os.path.join(directory, "reference.fa")
            with open(fasta, "w") as out:
                out.write(">chrU\n" + "A" * 2000 + "\n")
                out.write(">chrM\n" + shared + "N" * 400 + shared + "\n")
            targets = os.path.join(directory, "targets.fa")
            write_fasta(targets, [
                ("g1_Ref_h1_chrU_0_1500", "A" * 1500),
                ("g2_Ref_h1_chrU_300_1800", "A" * 1500),
                ("g3_Ref_h1_chrM_0_1300", shared),
                ("g4_Ref_h1_chrM_1700_3000", shared),
            ])
            bed = os.path.join(directory, "filtered.bed")
            with open(bed, "w") as out:
                out.write(bed_row("chrU", 0, 1500, "g1_Ref_h1_chrU_0_1500"))
                out.write(bed_row("chrU", 300, 1800, "g2_Ref_h1_chrU_300_1800"))
                out.write(bed_row("chrM", 0, 1300, "g3_Ref_h1_chrM_0_1300"))
                out.write(bed_row("chrM", 1700, 3000, "g4_Ref_h1_chrM_1700_3000"))
            output = os.path.join(directory, "merged.bed")
            groups = os.path.join(directory, "groups.tsv")
            pairs = os.path.join(directory, "pairs.tsv")
            subprocess.run([
                BINARY, "-f", fasta, "-T", targets, "-b", bed, "-o", output,
                "-g", groups, "-p", pairs, "-t", "2",
            ], check=True)

            with open(groups) as handle:
                rows = dict(
                    line.rstrip("\n").split("\t")
                    for line in handle
                    if line.strip() and not line.startswith("partition\t")
                )
            self.assertEqual(rows["g1_Ref_h1_chrU_0_1500"], "merged_g1g2")
            self.assertEqual(rows["g2_Ref_h1_chrU_300_1800"], "merged_g1g2")
            self.assertEqual(rows["g3_Ref_h1_chrM_0_1300"], "merged_g3g4")
            self.assertEqual(rows["g4_Ref_h1_chrM_1700_3000"], "merged_g3g4")
            with open(pairs) as handle:
                pair_rows = handle.read()
            self.assertIn(
                "g3_Ref_h1_chrM_0_1300\tg4_Ref_h1_chrM_1700_3000\t",
                pair_rows,
            )

    def test_alignment_only_skips_kmer_merge(self):
        if not os.path.isfile(BINARY):
            self.skipTest("run make before the test")
        with tempfile.TemporaryDirectory() as directory:
            rng = random.Random(7)
            shared = "".join(rng.choice("ACGT") for _ in range(1300))
            fasta = os.path.join(directory, "reference.fa")
            with open(fasta, "w") as out:
                out.write(">chrU\n" + "A" * 2000 + "\n")
                out.write(">chrM\n" + shared + "N" * 400 + shared + "\n")
            targets = os.path.join(directory, "targets.fa")
            write_fasta(targets, [
                ("g1_Ref_h1_chrU_0_1500", "A" * 1500),
                ("g2_Ref_h1_chrU_300_1800", "A" * 1500),
                ("g3_Ref_h1_chrM_0_1300", shared),
                ("g4_Ref_h1_chrM_1700_3000", shared),
            ])
            bed = os.path.join(directory, "filtered.bed")
            with open(bed, "w") as out:
                out.write(bed_row("chrU", 0, 1500, "g1_Ref_h1_chrU_0_1500"))
                out.write(bed_row("chrU", 300, 1800, "g2_Ref_h1_chrU_300_1800"))
                out.write(bed_row("chrM", 0, 1300, "g3_Ref_h1_chrM_0_1300"))
                out.write(bed_row("chrM", 1700, 3000, "g4_Ref_h1_chrM_1700_3000"))
            output = os.path.join(directory, "merged.bed")
            groups = os.path.join(directory, "groups.tsv")
            completed = subprocess.run([
                BINARY, "-f", fasta, "-T", targets, "-b", bed, "-o", output,
                "-g", groups, "--alignment-only",
            ], check=True, capture_output=True, text=True)

            with open(groups) as handle:
                rows = dict(
                    line.rstrip("\n").split("\t")
                    for line in handle
                    if line.strip() and not line.startswith("partition\t")
                )
            self.assertEqual(rows["g1_Ref_h1_chrU_0_1500"], "merged_g1g2")
            self.assertEqual(rows["g2_Ref_h1_chrU_300_1800"], "merged_g1g2")
            self.assertEqual(rows["g3_Ref_h1_chrM_0_1300"], "g3_Ref_h1_chrM_0_1300")
            self.assertEqual(rows["g4_Ref_h1_chrM_1700_3000"], "g4_Ref_h1_chrM_1700_3000")
            self.assertIn("skipping unmasked k-mer", completed.stderr)

    def test_only_kmer_skips_aligned_region_merge(self):
        if not os.path.isfile(BINARY):
            self.skipTest("run make before the test")
        with tempfile.TemporaryDirectory() as directory:
            rng = random.Random(7)
            shared = "".join(rng.choice("ACGT") for _ in range(1300))
            fasta = os.path.join(directory, "reference.fa")
            with open(fasta, "w") as out:
                out.write(">chrU\n" + "A" * 2000 + "\n")
                out.write(">chrM\n" + shared + "N" * 400 + shared + "\n")
            targets = os.path.join(directory, "targets.fa")
            write_fasta(targets, [
                ("g1_Ref_h1_chrU_0_1500", "A" * 1500),
                ("g2_Ref_h1_chrU_300_1800", "A" * 1500),
                ("g3_Ref_h1_chrM_0_1300", shared),
                ("g4_Ref_h1_chrM_1700_3000", shared),
            ])
            bed = os.path.join(directory, "filtered.bed")
            with open(bed, "w") as out:
                # These overlap by >1000 unmasked bases, but uppercase-only
                # sequence supplies only one distinct unmasked k-mer.
                out.write(bed_row("chrU", 0, 1500, "g1_Ref_h1_chrU_0_1500"))
                out.write(bed_row("chrU", 300, 1800, "g2_Ref_h1_chrU_300_1800"))
                # These do not overlap, but their copied sequence has
                # >1000 shared canonical unmasked 31-mers.
                out.write(bed_row("chrM", 0, 1300, "g3_Ref_h1_chrM_0_1300"))
                out.write(bed_row("chrM", 1700, 3000, "g4_Ref_h1_chrM_1700_3000"))
            output = os.path.join(directory, "merged.bed")
            groups = os.path.join(directory, "groups.tsv")
            pairs = os.path.join(directory, "pairs.tsv")
            completed = subprocess.run([
                BINARY, "-f", fasta, "-T", targets, "-b", bed, "-o", output,
                "-g", groups, "-p", pairs, "--only-kmer", "-t", "2",
            ], check=True, capture_output=True, text=True)

            with open(groups) as handle:
                rows = dict(
                    line.rstrip("\n").split("\t")
                    for line in handle
                    if line.strip() and not line.startswith("partition\t")
                )
            self.assertEqual(rows["g1_Ref_h1_chrU_0_1500"], "g1_Ref_h1_chrU_0_1500")
            self.assertEqual(rows["g2_Ref_h1_chrU_300_1800"], "g2_Ref_h1_chrU_300_1800")
            self.assertEqual(rows["g3_Ref_h1_chrM_0_1300"], "merged_g3g4")
            self.assertEqual(rows["g4_Ref_h1_chrM_1700_3000"], "merged_g3g4")
            self.assertIn("skipping aligned-reference-base", completed.stderr)
            with open(pairs) as handle:
                self.assertIn(
                    "g3_Ref_h1_chrM_0_1300\t"
                    "g4_Ref_h1_chrM_1700_3000\t",
                    handle.read(),
                )

    def test_only_kmer_and_alignment_only_are_mutually_exclusive(self):
        if not os.path.isfile(BINARY):
            self.skipTest("run make before the test")
        completed = subprocess.run([
            BINARY, "-f", "unused.fa", "-T", "unused-targets.fa",
            "-b", "unused.bed",
            "-o", "unused.out", "--only-kmer", "--alignment-only",
        ], capture_output=True, text=True)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("mutually exclusive", completed.stderr)


if __name__ == "__main__":
    unittest.main()
