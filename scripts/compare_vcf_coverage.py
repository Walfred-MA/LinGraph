#!/usr/bin/env python3
"""Audit reference/query coverage represented by a graph-derived VCF.

Reference coverage comes from the VCF's 0-based, half-open
``##referenceCoverage`` records. Query coverage comes from column 3 of the
seven-column ``graphcigartoref_persample.py`` output. When the VCF declares
``##querySequenceAlphabetFilter=ACGT_only``, a graph-CIGAR row containing an
ambiguous query base is treated exactly as graphreftovcf treats it: the whole
row is excluded from query coverage and assigned to the ``N_gap`` cause.
Reference-free loci declared by ``##alternativeLocus`` use independent local
coordinates. Unless present in the reference FAI, they are excluded from the
reference assembly audit; their query alignments still count toward coverage.

FAI files alone do not contain bases. Indexed FASTA files are therefore needed
to count masking and identify ambiguous/N gaps. If --reference-fasta or
--query-fasta is omitted, the script removes the trailing ``.fai`` from the
corresponding FAI path and uses that file.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple

from minsetref_core import IndexedFasta
from minsetref_segments import complement_intervals, merge_intervals, subtract_intervals


Interval = Tuple[int, int]
REFERENCE_CATEGORIES = ("Y_chrom", "constitutive_satellite", "other")
QUERY_CATEGORIES = ("N_gap", "scaffold_edge", "other")
COORD_RE = re.compile(r"^[<>]?(.+):(\d+)-(\d+)([+-]?)$")
REFERENCE_COVERAGE_RE = re.compile(r"^##referenceCoverage=<(.+)>$")
META_FIELD_RE = re.compile(
    r'(\w+)=(?:"((?:[^"\\]|\\.)*)"|([^,>]*))(?:,|$)'
)
AMBIGUOUS_BYTES_RE = re.compile(br"[^ACGTacgt]+")


def log(message: str) -> None:
    print(f"[compare_vcf_coverage] {message}", file=sys.stderr)


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def parse_coord(text: str) -> Optional[Tuple[str, int, int, str]]:
    match = COORD_RE.match((text or "").strip())
    if match is None:
        return None
    start = int(match.group(2))
    end = int(match.group(3))
    if end < start:
        start, end = end, start
    return match.group(1), start, end, match.group(4) or "+"


def intersect_intervals(left: Sequence[Interval], right: Sequence[Interval]) -> List[Interval]:
    left_values = merge_intervals(left)
    right_values = merge_intervals(right)
    output: List[Interval] = []
    i = 0
    j = 0
    while i < len(left_values) and j < len(right_values):
        left_start, left_end = left_values[i]
        right_start, right_end = right_values[j]
        start = max(left_start, right_start)
        end = min(left_end, right_end)
        if end > start:
            output.append((start, end))
        if left_end <= right_end:
            i += 1
        else:
            j += 1
    return output


def interval_size(intervals: Iterable[Interval]) -> int:
    return sum(end - start for start, end in intervals)


def clip_interval(start: int, end: int, length: int, context: str) -> Interval:
    if end < start:
        start, end = end, start
    if end <= 0 or start >= length:
        return 0, 0
    clipped_start = max(0, start)
    clipped_end = min(length, end)
    if clipped_start != start or clipped_end != end:
        log(
            f"warning: clipped {context} {start}-{end} to "
            f"{clipped_start}-{clipped_end}"
        )
    return clipped_start, clipped_end


def build_fasta_contig_aliases(
    reader: IndexedFasta,
) -> Dict[str, str]:
    """Map FASTA names and NCBI ``chromosome N`` header aliases to FAI names.

    FAI offsets point immediately after each header. Reading a small block
    backward from those offsets recovers descriptions without scanning the
    multi-gigabase FASTA. This lets a ``chr1`` BED match an indexed record such
    as ``NC_060925.1`` when its header says ``chromosome 1``.
    """
    aliases: Dict[str, str] = {}

    def add(alias: str, name: str) -> None:
        if not alias:
            return
        previous = aliases.get(alias)
        if previous is not None and previous != name:
            raise ValueError(
                f"contig alias {alias!r} is ambiguous between "
                f"{previous!r} and {name!r}"
            )
        aliases[alias] = name

    for name in reader.names():
        add(name, name)
        add(name.lower(), name)
        _length, sequence_offset, _line_bases, _line_width = reader.index[name]
        read_start = max(0, sequence_offset - 8192)
        raw = os.pread(reader._fd, sequence_offset - read_start, read_start)
        text = raw.decode("utf-8", "replace").rstrip("\r\n")
        header = text.splitlines()[-1] if text else ""
        if header.startswith(">"):
            chromosome = re.search(
                r"\bchromosome\s+([0-9]+|X|Y|M|MT)\b",
                header,
                flags=re.IGNORECASE,
            )
            if chromosome is not None:
                token = chromosome.group(1).upper()
                token = "M" if token == "MT" else token
                for alias in (token, f"chr{token}"):
                    add(alias, name)
                    add(alias.lower(), name)
    return aliases


def resolve_contig(name: str, aliases: Mapping[str, str]) -> Optional[str]:
    return aliases.get(name) or aliases.get(name.lower())


def parse_meta_fields(payload: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    position = 0
    while position < len(payload):
        match = META_FIELD_RE.match(payload, position)
        if match is None:
            raise ValueError(f"malformed VCF metadata payload: {payload!r}")
        value = match.group(2) if match.group(2) is not None else match.group(3)
        fields[match.group(1)] = bytes(value, "utf-8").decode("unicode_escape")
        position = match.end()
    return fields


@dataclass
class VcfCoverage:
    by_sample: Dict[str, Dict[str, List[Interval]]]
    query_alphabet_filter: str
    alternative_loci: Set[str] = field(default_factory=set)


def read_vcf_reference_coverage(path: str) -> VcfCoverage:
    by_sample: Dict[str, Dict[str, List[Interval]]] = defaultdict(
        lambda: defaultdict(list)
    )
    alphabet_filter = ""
    coordinate_system = ""
    interval_count = None
    header_samples = []
    alternative_loci: Set[str] = set()
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\r\n")
            if not line.startswith("#"):
                break
            if line.startswith("#CHROM\t"):
                header_samples = line.split("\t")[9:]
            if line.startswith("##alternativeLocus=<") and line.endswith(">"):
                fields = parse_meta_fields(line[len("##alternativeLocus=<"):-1])
                if fields.get("ID"):
                    alternative_loci.add(fields["ID"])
            if line.startswith("##referenceCoverageIntervalCount="):
                interval_count = int(line.split("=", 1)[1])
            if line.startswith("##querySequenceAlphabetFilter="):
                alphabet_filter = line.split("=", 1)[1].strip()
            elif line.startswith("##referenceCoverageCoordinateSystem="):
                coordinate_system = line.split("=", 1)[1].strip()
            else:
                match = REFERENCE_COVERAGE_RE.match(line)
                if match is None:
                    continue
                fields = parse_meta_fields(match.group(1))
                try:
                    sample = fields["Sample"]
                    chrom = fields["Chrom"]
                    start = int(fields["Start"])
                    end = int(fields["End"])
                except (KeyError, ValueError) as error:
                    raise ValueError(
                        f"malformed referenceCoverage header in {path}: {line}"
                    ) from error
                if end > start:
                    by_sample[sample][chrom].append((start, end))
    if coordinate_system and coordinate_system != "0-based-half-open":
        raise ValueError(
            f"unsupported reference coverage coordinate system: "
            f"{coordinate_system!r}"
        )
    if interval_count is not None and (by_sample or interval_count == 0):
        for sample in header_samples:
            by_sample.setdefault(sample, {})
    if not by_sample:
        raise ValueError(f"no ##referenceCoverage intervals found in {path}")
    normalized: Dict[str, Dict[str, List[Interval]]] = {}
    for sample, by_contig in by_sample.items():
        normalized[sample] = {
            contig: merge_intervals(intervals)
            for contig, intervals in by_contig.items()
        }
    return VcfCoverage(normalized, alphabet_filter, alternative_loci)


def choose_sample(coverage: VcfCoverage, requested: Optional[str]) -> str:
    samples = sorted(coverage.by_sample)
    if requested:
        if requested not in coverage.by_sample:
            raise ValueError(
                f"sample {requested!r} is absent from VCF reference coverage; "
                f"available: {','.join(samples)}"
            )
        return requested
    if len(samples) != 1:
        raise ValueError(
            "VCF contains multiple reference-coverage samples; select one "
            f"with --sample: {','.join(samples)}"
        )
    return samples[0]


def read_graphcigar_query_rows(
    path: str,
    lengths: Mapping[str, int],
    aliases: Mapping[str, str],
) -> List[Tuple[str, int, int]]:
    rows: List[Tuple[str, int, int]] = []
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 7:
                raise ValueError(
                    f"{path}:{line_number}: expected the seven-column "
                    "graphcigartoref_persample format"
                )
            coordinate = parse_coord(fields[2])
            if coordinate is None:
                raise ValueError(
                    f"{path}:{line_number}: invalid query coordinate "
                    f"{fields[2]!r}"
                )
            input_contig, start, end, _strand = coordinate
            contig = resolve_contig(input_contig, aliases)
            if contig is None or contig not in lengths:
                raise ValueError(
                    f"{path}:{line_number}: query contig {input_contig!r} is absent "
                    "from query FAI"
                )
            start, end = clip_interval(
                start, end, lengths[contig],
                f"{path}:{line_number} query interval",
            )
            if end > start:
                rows.append((contig, start, end))
    return rows


def read_bed(
    path: str,
    lengths: Mapping[str, int],
    aliases: Mapping[str, str],
) -> Dict[str, List[Interval]]:
    by_contig: Dict[str, List[Interval]] = defaultdict(list)
    unknown = set()
    with open_text(path) as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected BED3 or wider")
            input_contig = fields[0]
            contig = resolve_contig(input_contig, aliases)
            if contig is None or contig not in lengths:
                unknown.add(input_contig)
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: non-integer BED coordinates"
                ) from error
            start, end = clip_interval(
                start, end, lengths[contig], f"{path}:{line_number}",
            )
            if end > start:
                by_contig[contig].append((start, end))
    if unknown:
        examples = ",".join(sorted(unknown)[:10])
        log(
            f"warning: ignored BED contigs absent from reference FAI: "
            f"{examples}{'...' if len(unknown) > 10 else ''}"
        )
    return {
        contig: merge_intervals(intervals)
        for contig, intervals in by_contig.items()
    }


def _find_ambiguous_runs_contig(
    fasta_path: str,
    contig: str,
    length: int,
    sequence_offset: int,
) -> Tuple[str, List[Interval]]:
    """Scan one FASTA contig line-by-line for maximal non-A/C/G/T runs.

    The .fai sequence offset points to the first base of ``contig``.  A worker
    therefore seeks exactly once, then streams forward through the FASTA lines
    until ``length`` sequence bases have been consumed.  This avoids repeated
    IndexedFasta.fetch() calls and avoids a Python loop over individual bases.

    Runs touching adjacent FASTA lines are merged because newlines are not part
    of sequence coordinates.
    """
    runs: List[Interval] = []
    position = 0

    # A private file handle per worker avoids seek/read races between threads.
    # Binary mode also avoids text decoding/newline-conversion overhead.
    with open(fasta_path, "rb", buffering=1024 * 1024) as handle:
        handle.seek(sequence_offset)

        while position < length:
            raw_line = handle.readline()
            if not raw_line:
                raise ValueError(
                    f"unexpected EOF while scanning FASTA contig {contig!r}: "
                    f"read {position} of {length} bases"
                )

            # FASTA sequence lines may end in LF or CRLF.  Do not use strip(),
            # because any other byte on a sequence line should remain visible
            # to the ambiguity detector.
            sequence = raw_line.rstrip(b"\r\n")
            if not sequence:
                continue
            if sequence.startswith(b">"):
                raise ValueError(
                    f"encountered FASTA header before end of contig {contig!r}: "
                    f"read {position} of {length} bases"
                )

            # The final physical FASTA line can be longer than the number of
            # sequence bases still expected only for a malformed/inconsistent
            # index.  Clip defensively to the indexed contig length.
            remaining = length - position
            if len(sequence) > remaining:
                sequence = sequence[:remaining]

            for match in AMBIGUOUS_BYTES_RE.finditer(sequence):
                run_start = position + match.start()
                run_end = position + match.end()

                # If the previous FASTA line ended in ambiguity and this line
                # starts in ambiguity, the two matches are one genomic run.
                if runs and run_start == runs[-1][1]:
                    runs[-1] = (runs[-1][0], run_end)
                else:
                    runs.append((run_start, run_end))

            position += len(sequence)

    return contig, runs


def find_ambiguous_runs(
    fasta_path: str,
    fai_index: Mapping[str, Tuple[int, int, int, int]],
    contigs: Sequence[str],
    threads: int,
) -> Dict[str, List[Interval]]:
    """Find maximal non-A/C/G/T runs, parallelized one contig per worker.

    ``fai_index`` entries are ``(length, sequence_offset, line_bases,
    line_width)`` as stored by :class:`IndexedFasta`.  Only length and offset
    are required because each worker streams FASTA lines directly.
    """
    names = list(contigs)
    if not names:
        return {}

    workers = max(1, min(threads, len(names)))
    results: Dict[str, List[Interval]] = {}

    def submit_args(contig: str) -> Tuple[str, str, int, int]:
        try:
            length, sequence_offset, _line_bases, _line_width = fai_index[contig]
        except KeyError as error:
            raise ValueError(f"contig {contig!r} is absent from query FAI") from error
        return fasta_path, contig, length, sequence_offset

    if workers == 1:
        for contig in names:
            name, runs = _find_ambiguous_runs_contig(*submit_args(contig))
            if runs:
                results[name] = runs
    else:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(_find_ambiguous_runs_contig, *submit_args(contig)): contig
                for contig in names
            }
            for future in as_completed(futures):
                name, runs = future.result()
                if runs:
                    results[name] = runs

    # Keep FASTA order deterministic even though workers finish out of order.
    return {name: results[name] for name in names if name in results}


def interval_overlaps_sorted(
    start: int,
    end: int,
    intervals: Sequence[Interval],
    starts: Sequence[int],
) -> bool:
    if not intervals or end <= start:
        return False
    index = bisect.bisect_left(starts, end) - 1
    return index >= 0 and intervals[index][1] > start


def complement_by_contig(
    covered: Mapping[str, Sequence[Interval]],
    lengths: Mapping[str, int],
) -> Dict[str, List[Interval]]:
    return {
        contig: complement_intervals(covered.get(contig, ()), 0, length)
        for contig, length in lengths.items()
        if length > 0
    }


def categorize_reference_missing(
    missing: Mapping[str, Sequence[Interval]],
    satellites: Mapping[str, Sequence[Interval]],
    y_contigs: set,
) -> Dict[str, Dict[str, List[Interval]]]:
    output = {category: defaultdict(list) for category in REFERENCE_CATEGORIES}
    for contig, intervals in missing.items():
        if contig in y_contigs:
            output["Y_chrom"][contig].extend(intervals)
            continue
        satellite = intersect_intervals(intervals, satellites.get(contig, ()))
        other = subtract_intervals(intervals, satellite)
        output["constitutive_satellite"][contig].extend(satellite)
        output["other"][contig].extend(other)
    return {
        category: {
            contig: merge_intervals(intervals)
            for contig, intervals in by_contig.items() if intervals
        }
        for category, by_contig in output.items()
    }


def categorize_query_missing(
    missing: Mapping[str, Sequence[Interval]],
    ambiguous_runs: Mapping[str, Sequence[Interval]],
    rejected_rows: Mapping[str, Sequence[Interval]],
    lengths: Mapping[str, int],
) -> Dict[str, Dict[str, List[Interval]]]:
    output = {category: defaultdict(list) for category in QUERY_CATEGORIES}
    for contig, intervals in missing.items():
        n_mask = merge_intervals(
            list(ambiguous_runs.get(contig, ()))
            + list(rejected_rows.get(contig, ()))
        )
        n_gap = intersect_intervals(intervals, n_mask)
        non_n = subtract_intervals(intervals, n_gap)
        terminal_components = [
            (start, end) for start, end in intervals
            if start == 0 or end == lengths[contig]
        ]
        scaffold_edge = intersect_intervals(non_n, terminal_components)
        other = subtract_intervals(non_n, scaffold_edge)
        output["N_gap"][contig].extend(n_gap)
        output["scaffold_edge"][contig].extend(scaffold_edge)
        output["other"][contig].extend(other)
    return {
        category: {
            contig: merge_intervals(intervals)
            for contig, intervals in by_contig.items() if intervals
        }
        for category, by_contig in output.items()
    }


@dataclass
class CategoryStats:
    intervals: int = 0
    bases: int = 0
    unmasked: int = 0
    masked: int = 0
    ambiguous: int = 0

    def add(self, other: "CategoryStats") -> None:
        self.intervals += other.intervals
        self.bases += other.bases
        self.unmasked += other.unmasked
        self.masked += other.masked
        self.ambiguous += other.ambiguous


def sequence_stats(
    reader: IndexedFasta,
    contig: str,
    start: int,
    end: int,
    chunk_size: int,
) -> CategoryStats:
    stats = CategoryStats(intervals=1, bases=end - start)
    for chunk_start in range(start, end, chunk_size):
        sequence = reader.fetch(
            contig, chunk_start, min(end, chunk_start + chunk_size),
        )
        stats.unmasked += (
            sequence.count("A") + sequence.count("C")
            + sequence.count("G") + sequence.count("T")
        )
        stats.masked += (
            sequence.count("a") + sequence.count("c")
            + sequence.count("g") + sequence.count("t")
        )
        stats.ambiguous += len(sequence) - (
            sequence.count("A") + sequence.count("C")
            + sequence.count("G") + sequence.count("T")
            + sequence.count("a") + sequence.count("c")
            + sequence.count("g") + sequence.count("t")
        )
    return stats


def write_category_bed(
    path: str,
    source: str,
    categories: Mapping[str, Mapping[str, Sequence[Interval]]],
    reader: IndexedFasta,
    chunk_size: int,
) -> Tuple[Dict[str, CategoryStats], CategoryStats]:
    summaries = {category: CategoryStats() for category in categories}
    total = CategoryStats()
    with open(path, "wt") as output:
        output.write(
            "#chrom\tstart\tend\tcategory\tscore\tstrand\t"
            "source\tlength\tunmasked_bases\tmasked_bases\tambiguous_bases\n"
        )
        for category, by_contig in categories.items():
            for contig in reader.names():
                for start, end in by_contig.get(contig, ()):
                    stats = sequence_stats(
                        reader, contig, start, end, chunk_size,
                    )
                    summaries[category].add(stats)
                    total.add(stats)
                    output.write(
                        f"{contig}\t{start}\t{end}\t{category}\t0\t.\t"
                        f"{source}\t{stats.bases}\t{stats.unmasked}\t"
                        f"{stats.masked}\t{stats.ambiguous}\n"
                    )
    return summaries, total


def write_summary(
    path: str,
    sample: str,
    reference_stats: Mapping[str, CategoryStats],
    reference_total: CategoryStats,
    query_stats: Mapping[str, CategoryStats],
    query_total: CategoryStats,
) -> None:
    with open(path, "wt") as output:
        output.write(f"#sample\t{sample}\n")
        output.write(
            "scope\tcategory\tintervals\ttotal_bases\tunmasked_bases\t"
            "masked_bases\tambiguous_bases\n"
        )
        for scope, category, stats in (
            [("reference", "total_unannotated", reference_total)]
            + [("reference", category, reference_stats[category])
               for category in REFERENCE_CATEGORIES]
            + [("query", "total_unannotated", query_total)]
            + [("query", category, query_stats[category])
               for category in QUERY_CATEGORIES]
        ):
            output.write(
                f"{scope}\t{category}\t{stats.intervals}\t{stats.bases}\t"
                f"{stats.unmasked}\t{stats.masked}\t{stats.ambiguous}\n"
            )


def resolve_fasta(fai_path: str, fasta_path: Optional[str], label: str) -> str:
    if not os.path.isfile(fai_path):
        raise FileNotFoundError(fai_path)
    if fasta_path is None:
        if not fai_path.endswith(".fai"):
            raise ValueError(
                f"--{label}-fasta is required when --{label}-fai does not "
                "end in .fai"
            )
        fasta_path = fai_path[:-4]
    if not os.path.isfile(fasta_path):
        raise FileNotFoundError(
            f"{label} FASTA not found: {fasta_path}; an FAI alone cannot "
            "identify masking or N gaps"
        )
    return fasta_path


def parse_y_contigs(
    value: str,
    reference_names: Sequence[str],
    aliases: Mapping[str, str],
) -> set:
    requested = [token.strip() for token in value.split(",") if token.strip()]
    requested.extend(
        name for name in reference_names
        if name.lower() in {"y", "chry", "chromosomey", "chromosome_y"}
    )
    matched = {
        resolved for token in requested
        for resolved in (resolve_contig(token, aliases),)
        if resolved is not None
    }
    if not matched:
        log(
            "warning: no Y contig matched; supply accession/name(s) with "
            "--y-contigs if this reference contains Y"
        )
    return matched


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare graph-derived VCF reference coverage and graph-CIGAR "
            "query coverage against indexed FASTA assemblies"
        )
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--graphcigar", required=True)
    parser.add_argument("--reference-fai", required=True)
    parser.add_argument("--query-fai", required=True)
    parser.add_argument("--reference-fasta")
    parser.add_argument("--query-fasta")
    parser.add_argument("--satellite-bed", help="optional satellite annotation BED; omit for unannotated coverage")
    parser.add_argument("--sample", help="VCF coverage sample; inferred when unique")
    parser.add_argument(
        "--y-contigs", default="chrY,Y",
        help="comma-separated Y contig names/accessions [chrY,Y]",
    )
    parser.add_argument("-o", "--output-prefix", required=True)
    parser.add_argument(
        "--chunk-size", type=int, default=4_000_000,
        help="FASTA scan/fetch chunk size [4000000]",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=min(8, os.cpu_count() or 1),
        help="threads for line-by-line query FASTA ambiguous/N scan [min(8, CPU count)]",
    )
    args = parser.parse_args()
    if args.chunk_size < 1:
        parser.error("--chunk-size must be positive")
    if args.threads < 1:
        parser.error("--threads must be positive")
    for path in (args.vcf, args.graphcigar, args.satellite_bed):
        if path is None:
            continue
        if not os.path.isfile(path):
            parser.error(f"file not found: {path}")

    try:
        reference_fasta = resolve_fasta(
            args.reference_fai, args.reference_fasta, "reference",
        )
        query_fasta = resolve_fasta(
            args.query_fai, args.query_fasta, "query",
        )
        output_directory = os.path.dirname(os.path.abspath(args.output_prefix))
        os.makedirs(output_directory, exist_ok=True)
        with IndexedFasta(reference_fasta, args.reference_fai) as reference_reader, \
             IndexedFasta(query_fasta, args.query_fai) as query_reader:
            reference_lengths = {
                name: reference_reader.length(name)
                for name in reference_reader.names()
            }
            query_lengths = {
                name: query_reader.length(name)
                for name in query_reader.names()
            }
            reference_aliases = build_fasta_contig_aliases(reference_reader)
            query_aliases = build_fasta_contig_aliases(query_reader)

            vcf_coverage = read_vcf_reference_coverage(args.vcf)
            sample = choose_sample(vcf_coverage, args.sample)
            reference_covered: Dict[str, List[Interval]] = defaultdict(list)
            excluded_alternative_loci = 0
            for input_contig, intervals in vcf_coverage.by_sample[sample].items():
                contig = resolve_contig(input_contig, reference_aliases)
                if contig is None or contig not in reference_lengths:
                    if input_contig in vcf_coverage.alternative_loci:
                        # SourceContig/SourceStart describe provenance, not a
                        # placement on the selected reference assembly.
                        excluded_alternative_loci += 1
                        continue
                    raise ValueError(
                        f"VCF reference-coverage contig {input_contig!r} is absent "
                        "from reference FAI"
                    )
                for start, end in intervals:
                    clipped = clip_interval(
                        start, end, reference_lengths[contig],
                        f"VCF reference coverage on {contig}",
                    )
                    if clipped[1] > clipped[0]:
                        reference_covered[contig].append(clipped)
            if excluded_alternative_loci:
                log(
                    f"excluded {excluded_alternative_loci} declared alternative-locus "
                    "contig(s) from reference assembly coverage; "
                    "their query alignments remain included"
                )
            reference_covered = {
                contig: merge_intervals(intervals)
                for contig, intervals in reference_covered.items()
            }
            reference_missing = complement_by_contig(
                reference_covered, reference_lengths,
            )
            satellites = read_bed(
                args.satellite_bed, reference_lengths, reference_aliases,
            ) if args.satellite_bed else {}
            y_contigs = parse_y_contigs(
                args.y_contigs, reference_reader.names(), reference_aliases,
            )
            reference_categories = categorize_reference_missing(
                reference_missing, satellites, y_contigs,
            )

            query_rows = read_graphcigar_query_rows(
                args.graphcigar, query_lengths, query_aliases,
            )
            scan_threads = min(args.threads, max(1, len(query_lengths)))
            log(
                f"scanning query FASTA for ambiguous/N runs across "
                f"{len(query_lengths)} contigs with {scan_threads} threads"
            )
            ambiguous_runs = find_ambiguous_runs(
                query_fasta,
                query_reader.index,
                query_reader.names(),
                scan_threads,
            )
            ambiguous_starts = {
                contig: [start for start, _end in intervals]
                for contig, intervals in ambiguous_runs.items()
            }
            acgt_only = (
                vcf_coverage.query_alphabet_filter.lower() == "acgt_only"
            )
            accepted_query: Dict[str, List[Interval]] = defaultdict(list)
            rejected_query: Dict[str, List[Interval]] = defaultdict(list)
            for contig, start, end in query_rows:
                rejected = (
                    acgt_only
                    and interval_overlaps_sorted(
                        start,
                        end,
                        ambiguous_runs.get(contig, ()),
                        ambiguous_starts.get(contig, ()),
                    )
                )
                target = rejected_query if rejected else accepted_query
                target[contig].append((start, end))
            accepted_query = {
                contig: merge_intervals(intervals)
                for contig, intervals in accepted_query.items()
            }
            rejected_query = {
                contig: merge_intervals(intervals)
                for contig, intervals in rejected_query.items()
            }
            query_missing = complement_by_contig(
                accepted_query, query_lengths,
            )
            query_categories = categorize_query_missing(
                query_missing,
                ambiguous_runs,
                rejected_query,
                query_lengths,
            )

            reference_bed = args.output_prefix + ".missing_reference.bed"
            query_bed = args.output_prefix + ".missing_query.bed"
            summary_path = args.output_prefix + ".summary.tsv"
            reference_stats, reference_total = write_category_bed(
                reference_bed,
                "reference",
                reference_categories,
                reference_reader,
                args.chunk_size,
            )
            query_stats, query_total = write_category_bed(
                query_bed,
                "query",
                query_categories,
                query_reader,
                args.chunk_size,
            )
            write_summary(
                summary_path,
                sample,
                reference_stats,
                reference_total,
                query_stats,
                query_total,
            )

        log(
            f"reference unannotated={reference_total.bases} bp; "
            f"query unannotated={query_total.bases} bp"
        )
        log(f"summary: {summary_path}")
        log(f"missing reference BED: {reference_bed}")
        log(f"missing query BED: {query_bed}")
        return 0
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
