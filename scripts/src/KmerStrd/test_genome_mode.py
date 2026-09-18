#!/usr/bin/env python3
from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "src" / "main.cpp"


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


class KmerStrdGenomeModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.build_directory = tempfile.TemporaryDirectory()
        cls.binary = Path(cls.build_directory.name) / "kmerstrd"
        subprocess.run(
            [
                "c++", "-std=c++17", "-O2", "-Wall", "-Wextra",
                "-pthread", str(SOURCE), "-o", str(cls.binary),
            ],
            check=True,
        )

    @classmethod
    def tearDownClass(cls) -> None:
        cls.build_directory.cleanup()

    def test_genome_mode_updates_orientation_and_preserves_failures(self) -> None:
        reference_sequence = (
            "ACGTTGCATGTCAGTACGATCGTACCTGAGCATTCGATGCTAGCTACGATCG"
            "TACGGCATATCGTACGATGACCTAGTCGATCGTAGCATGCTAGTCAGTACGA"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fa"
            query = root / "query.fa"
            source = root / "input.tsv"
            output = root / "output.tsv"
            reference.write_text(f">ref\n{reference_sequence}\n")
            query.write_text(
                f">same\n{reference_sequence}\n"
                f">reverse\n{reverse_complement(reference_sequence)}\n"
            )
            successful_same = (
                f"same-row\tx\tsame:0-{len(reference_sequence)}\t"
                f"ref:0-{len(reference_sequence)}-\tpart_1"
            )
            successful_reverse = (
                f"reverse-row\tx\treverse:0-{len(reference_sequence)}+\t"
                f"ref:0-{len(reference_sequence)}+\tpart_2"
            )
            successful_union = (
                "union-row\tx\tsame:0-45-;same:60-100+\t"
                "ref:0-45+;ref:60-100+\tpart_3"
            )
            filtered = (
                f"filtered-row\tx\treverse:0-{len(reference_sequence)}+\t"
                f"ref:0-{len(reference_sequence)}+\tsingle"
            )
            failed = "failed-row\tx\tsame:0-10-\tref:0-10+\tpart_4"
            source.write_bytes((
                successful_same + "\r\n" +
                successful_reverse + "\r\n" +
                successful_union + "\r\n" +
                filtered + "\r\n" +
                failed
            ).encode())

            completed = subprocess.run(
                [
                    str(self.binary), "--genome",
                    "-i", str(source), "-o", str(output),
                    "-r", str(reference), "-q", str(query),
                    "-rc", "4", "-qc", "3", "-oc", "3", "-t", "4",
                    "--grep", r"\tpart_",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(
                output.read_bytes(),
                (
                    successful_same.replace(
                        f"same:0-{len(reference_sequence)}",
                        f"same:0-{len(reference_sequence)}-",
                    ) + "\r\n" +
                    successful_reverse.replace(
                        f"reverse:0-{len(reference_sequence)}+",
                        f"reverse:0-{len(reference_sequence)}-",
                    ) + "\r\n" +
                    successful_union.replace(
                        "same:0-45-;same:60-100+",
                        "same:0-45+;same:60-100+",
                    ) + "\r\n" +
                    filtered + "\r\n" +
                    failed
                ).encode(),
            )
            self.assertIn("updated 3 row(s)", completed.stderr)
            self.assertIn("left 1 failed/ambiguous row(s) unchanged", completed.stderr)
            self.assertIn("skipped 1 filtered/comment/blank row(s)", completed.stderr)

    def test_grep_column_restricts_matching_to_one_based_column(self) -> None:
        reference_sequence = (
            "ACGTTGCATGTCAGTACGATCGTACCTGAGCATTCGATGCTAGCTACGATCG"
            "TACGGCATATCGTACGATGACCTAGTCGATCGTAGCATGCTAGTCAGTACGA"
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fa"
            query = root / "query.fa"
            source = root / "input.tsv"
            output = root / "output.tsv"
            reference.write_text(f">ref\n{reference_sequence}\n")
            query.write_text(f">same\n{reference_sequence}\n")
            selected = (
                f"selected\tx\tsame:0-{len(reference_sequence)}\t"
                f"ref:0-{len(reference_sequence)}-\t0\t.\tpart_2"
            )
            pattern_in_wrong_column = (
                f"part_elsewhere\tx\tsame:0-{len(reference_sequence)}\t"
                f"ref:0-{len(reference_sequence)}-\t0\t.\tprimary"
            )
            missing_column = (
                f"part_short\tx\tsame:0-{len(reference_sequence)}\t"
                f"ref:0-{len(reference_sequence)}-"
            )
            source.write_text(
                "\n".join((selected, pattern_in_wrong_column, missing_column))
                + "\n"
            )

            completed = subprocess.run(
                [
                    str(self.binary), "--genome",
                    "-i", str(source), "-o", str(output),
                    "-r", str(reference), "-q", str(query),
                    "-rc", "4", "-qc", "3", "-oc", "3",
                    "--grep", "part_", "--grepcol", "7", "-t", "2",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            expected_selected = selected.replace(
                f"same:0-{len(reference_sequence)}",
                f"same:0-{len(reference_sequence)}-",
            )
            self.assertEqual(
                output.read_text().splitlines(),
                [expected_selected, pattern_in_wrong_column, missing_column],
            )
            self.assertIn("updated 1 row(s)", completed.stderr)
            self.assertIn("skipped 2 filtered/comment/blank row(s)", completed.stderr)

    def test_grep_column_requires_grep(self) -> None:
        completed = subprocess.run(
            [
                str(self.binary), "--genome", "-i", "input.tsv",
                "-o", "output.tsv", "-r", "reference.fa",
                "-q", "query.fa", "-rc", "4", "-qc", "3",
                "-oc", "3", "--grepcol", "7",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        self.assertNotEqual(completed.returncode, 0)

    def test_genome_strand_truth_table_matches_graphcigartoref_contract(
        self,
    ) -> None:
        """KmerStrd disagreement must equal graphcigartoref reverse view."""
        reference_sequence = (
            "ACGTTGCATGTCAGTACGATCGTACCTGAGCATTCGATGCTAGCTACGATCG"
            "TACGGCATATCGTACGATGACCTAGTCGATCGTAGCATGCTAGTCAGTACGA"
        )
        reverse_sequence = reverse_complement(reference_sequence)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = root / "reference.fa"
            query = root / "query.fa"
            source = root / "input.tsv"
            output = root / "output.tsv"
            reference.write_text(f">ref\n{reference_sequence}\n")
            query.write_text(
                f">same\n{reference_sequence}\n"
                f">reverse\n{reverse_sequence}\n"
            )
            rows = []
            expected_reverse = []
            for ref_strand in ("+", "-"):
                for query_name, is_reverse in (
                    ("same", False), ("reverse", True),
                ):
                    rows.append(
                        f"{query_name}-{ref_strand}\tx\t"
                        f"{query_name}:0-{len(reference_sequence)}+\t"
                        f"ref:0-{len(reference_sequence)}{ref_strand}\t"
                        "part_contract"
                    )
                    expected_reverse.append(is_reverse)
            source.write_text("\n".join(rows) + "\n")
            completed = subprocess.run(
                [
                    str(self.binary), "--genome",
                    "-i", str(source), "-o", str(output),
                    "-r", str(reference), "-q", str(query),
                    "-rc", "4", "-qc", "3", "-oc", "3",
                    "--grep", r"\tpart_", "-t", "4",
                ],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            actual_rows = output.read_text().splitlines()
            self.assertEqual(len(actual_rows), len(rows))
            for row, should_reverse in zip(
                actual_rows, expected_reverse,
            ):
                fields = row.split("\t")
                query_strand = fields[2][-1]
                reference_strand = fields[3][-1]
                # graphcigartoref.py uses precisely this inequality when it
                # selects reverse_ref_view for the pair.
                self.assertEqual(
                    query_strand != reference_strand,
                    should_reverse,
                    row,
                )

    def test_legacy_mode_remains_available(self) -> None:
        sequence = "ACGTTGCATGTCAGTACGATCGTACCTGAGCATTCGATGC"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fa"
            output = root / "output.fa"
            source.write_text(f">record contig:0-{len(sequence)}+\n{sequence}\n")
            completed = subprocess.run(
                [str(self.binary), "-i", str(source), "-o", str(output), "-t", "2"],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertEqual(output.read_text(), source.read_text())


if __name__ == "__main__":
    unittest.main()
