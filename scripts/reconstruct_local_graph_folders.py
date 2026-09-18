#!/usr/bin/env python3
"""Reconstruct per-partition graph folders from a local-graph package.

Inputs are ``local_graphs.tsv``, ``alternatives.fasta``, and an original
reference FASTA.  Reference and alternative sequences are loaded once into
shared RAM.  Every local graph is then written independently as::

    OUTPUT/Graphs/PARTITION/PARTITION.bed
    OUTPUT/Graphs/PARTITION/PARTITION.FA
    OUTPUT/Graphs/PARTITION/PARTITION.FA_db.*  # BLAST database
    OUTPUT/Graphs/PARTITION/PARTITION.fasta  # original partition templates
    OUTPUT/Graphs/PARTITION/PARTITIONcache.json  # when packed
    OUTPUT/summary -> input package directory
    OUTPUT/references -> input package directory/references  # when present
    OUTPUT/Graphs.list  # graph alignment paths (.FA)
    OUTPUT/Graphs.template.list  # hotspot targets (.fasta)
    OUTPUT/Graphs.template.list.bin  # shared KmerSearcher template cache

The BED contains the original partition-defining intervals. ``PARTITION.FA``
contains the refined reference path followed by all retained alternative/novel
paths. Refinement therefore changes graph sequence paths without silently
redefining the partition's original reference mask. Alternative headers are
copied verbatim from ``alternatives.fasta``. Paths compacted as verified
``reference_slice`` rows are fetched from the supplied reference and recover
their exact headers from ``local_graphs.tsv``.
``PARTITION.fasta`` contains only the template path(s) represented by the
partition BED. ``--original_fasta`` can supply their original chromosome
sequence when it differs from the main reference used for graph paths.
No cohort assemblies are read. Non-reference original intervals are recovered
from packed paths on their own source haplotype/contig; an incomplete package
must be corrected by its producer, not by requiring downstream assemblies.
``partition_caches.jsonl`` is optional for compatibility with older packages.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import gc
import hashlib
import json
import logging
import os
import subprocess
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from build_local_graphs import (
    PARTITION_CACHE_PACKAGE,
    REFERENCE_SLICE_RE,
    public_path_id_from_values,
    reverse_complement,
)
from graph_list_cache import ensure_graph_list_cache
from minsetref_core import wrap_fasta


LOG = logging.getLogger("reconstruct_local_graph_folders")
RUNMAKEBLASTDB = str(Path(__file__).resolve().with_name("runmakeblastdb"))

NC_TO_CHR = {
    "NC_060925.1": "chr1",
    "NC_060926.1": "chr2",
    "NC_060927.1": "chr3",
    "NC_060928.1": "chr4",
    "NC_060929.1": "chr5",
    "NC_060930.1": "chr6",
    "NC_060931.1": "chr7",
    "NC_060932.1": "chr8",
    "NC_060933.1": "chr9",
    "NC_060934.1": "chr10",
    "NC_060935.1": "chr11",
    "NC_060936.1": "chr12",
    "NC_060937.1": "chr13",
    "NC_060938.1": "chr14",
    "NC_060939.1": "chr15",
    "NC_060940.1": "chr16",
    "NC_060941.1": "chr17",
    "NC_060942.1": "chr18",
    "NC_060943.1": "chr19",
    "NC_060944.1": "chr20",
    "NC_060945.1": "chr21",
    "NC_060946.1": "chr22",
    "NC_060947.1": "chrX",
    "NC_060948.1": "chrY",
}
CHR_TO_NC = {chromosome: accession for accession, chromosome in NC_TO_CHR.items()}


def canonical_contig(name: str) -> str:
    """Return the NC_* name used by local_graphs.tsv when an alias is known."""
    return CHR_TO_NC.get(name, name)


def load_fasta_ram(path: str) -> Dict[str, Tuple[str, str]]:
    """Load FASTA as ID -> (complete header without '>', sequence)."""
    records: Dict[str, Tuple[str, str]] = {}
    description: Optional[str] = None
    pieces: List[str] = []

    def finish() -> None:
        if description is None:
            return
        record_id = description.split(None, 1)[0]
        if not record_id:
            raise ValueError(f"{path}: empty FASTA record ID")
        if record_id in records:
            raise ValueError(f"{path}: duplicate FASTA ID {record_id!r}")
        records[record_id] = (description, "".join(pieces))

    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                finish()
                description = raw[1:].rstrip("\r\n")
                pieces = []
            elif raw.strip():
                if description is None:
                    raise ValueError(
                        f"{path}:{line_number}: sequence precedes first header"
                    )
                pieces.append(raw.strip())
    finish()
    return records


def load_reference_ram(path: str) -> Dict[str, str]:
    raw = load_fasta_ram(path)
    output: Dict[str, str] = {}
    aliases = 0
    for record_id, (_description, sequence) in raw.items():
        canonical = canonical_contig(record_id)
        if canonical in output:
            raise ValueError(
                f"{path}: both {record_id!r} and another record normalize "
                f"to {canonical!r}"
            )
        output[canonical] = sequence
        aliases += canonical != record_id
    LOG.info(
        "Loaded %d reference records (%.2f GiB) into RAM; normalized %d "
        "chr* names to NC_*",
        len(output), sum(map(len, output.values())) / (1024 ** 3), aliases,
    )
    return output


def selected_reference_rows(
    rows: Sequence[Mapping[str, str]],
) -> List[Mapping[str, str]]:
    has_refined = any(row["type"] == "reference" for row in rows)
    selected_type = "reference" if has_refined else "original"
    return [row for row in rows if row["type"] == selected_type]


def original_definition_rows(
    rows: Sequence[Mapping[str, str]],
) -> List[Mapping[str, str]]:
    """Return immutable partition-defining rows, with legacy fallback."""
    original = [row for row in rows if row["type"] == "original"]
    return original or selected_reference_rows(rows)


def parsed_reference_slice(
    row: Mapping[str, str], context: str,
) -> Optional[Tuple[str, int, int, str]]:
    """Parse an optional verified reference-backed sequence descriptor."""
    encoded = row.get("reference_slice", ".")
    if encoded in {None, "", "."}:
        return None
    match = REFERENCE_SLICE_RE.fullmatch(encoded)
    if match is None:
        raise ValueError(f"{context}: invalid reference_slice {encoded!r}")
    contig, start_text, end_text, strand = match.groups()
    start, end = int(start_text), int(end_text)
    expected = int(row["source_end"]) - int(row["source_start"])
    if start < 0 or end <= start or end - start != expected:
        raise ValueError(
            f"{context}: reference_slice {encoded!r} does not match the "
            f"{expected}-base graph path"
        )
    return contig, start, end, strand


class PackedSourceIndex:
    """Find exact source intervals using sequences already in the package.

    Original BED intervals can be inside longer refined paths or span several
    stored paths. Index by haplotype AND contig so chr1 on HG38 can never be
    substituted with chr1 on CHM13. Sequence strings are not duplicated here.
    """

    def __init__(self, by_partition, alternatives):
        grouped = defaultdict(list)
        self.local = defaultdict(list)
        for partition, rows in by_partition.items():
            # Prefer selected graph rows when an original metadata row has the
            # same public ID (possibly with a different traversal strand).
            ordered = sorted(rows, key=lambda row: row['type'] == 'original')
            seen = set()
            for row in ordered:
                record_id = public_path_id_from_values(
                    partition, row['source_contig'], row['source_start'], row['source_end'],
                )
                if record_id in seen:
                    continue
                if record_id not in alternatives and parsed_reference_slice(
                    row, f'{partition}:{row["path_order"]}',
                ) is None:
                    continue
                seen.add(record_id)
                start, end = int(row['source_start']), int(row['source_end'])
                if start < 0 or end <= start or row['strand'] not in {'+', '-'}:
                    raise ValueError(f'{partition}: invalid packed source interval')
                key = (row['source_haplotype'], canonical_contig(row['source_contig']))
                grouped[key].append((start, end, partition, row))
                self.local[(partition, *key)].append((start, end, partition, row))
        for paths in self.local.values():
            paths.sort(key=lambda item: (item[0], item[1], int(item[3]['path_order'])))
        self.groups = {}
        for key, paths in grouped.items():
            paths.sort(key=lambda item: (item[0], item[1], item[2], int(item[3]['path_order'])))
            starts, maximum_ends = [], []
            maximum = -1
            for start, end, _partition, _row in paths:
                starts.append(start)
                maximum = max(maximum, end)
                maximum_ends.append(maximum)
            self.groups[key] = (paths, starts, maximum_ends)

    def covering_paths(self, row):
        key = (row['source_haplotype'], canonical_contig(row['source_contig']))
        start, end = int(row['source_start']), int(row['source_end'])
        local = [item for item in self.local.get((row.get('partition'), *key), ())
                 if item[0] < end and item[1] > start]
        exact = [item for item in local if item[0] == start and item[1] == end]
        if exact:
            return exact[:1]
        if self._covers(local, start, end):
            return local
        group = self.groups.get(key)
        if group is None:
            return []
        paths, starts, maximum_ends = group
        left = bisect.bisect_right(maximum_ends, start)
        right = bisect.bisect_left(starts, end)
        overlaps = [item for item in paths[left:right] if item[1] > start]
        return overlaps if self._covers(overlaps, start, end) else []

    @staticmethod
    def _covers(paths, start, end):
        cursor = start
        for low, high, _partition, _path in paths:
            if low > cursor:
                return False
            cursor = max(cursor, high)
        return cursor >= end

    def sequence(self, row, fetch_path):
        paths = self.covering_paths(row)
        if not paths:
            return None
        start, end = int(row['source_start']), int(row['source_end'])
        result = bytearray(end - start)
        cursor = start
        for low, high, partition, path in paths:
            _header, sequence = fetch_path(partition, path)
            if len(sequence) != high - low:
                raise ValueError(f'{partition}: packed sequence length does not match its source interval')
            a, b = max(start, low), min(end, high)
            # Extract only the needed slice before reversing a long path.
            fragment = (sequence[a - low:b - low] if path['strand'] == '+' else
                        reverse_complement(sequence[high - b:high - a]))
            fragment = fragment.encode('ascii')
            overlap_end = min(cursor, b)
            if overlap_end > a and result[a - start:overlap_end - start] != fragment[:overlap_end - a]:
                raise ValueError(
                    f'{partition}: packed paths disagree for '
                    f'{row["source_haplotype"]}:{row["source_contig"]}:{a}-{overlap_end}'
                )
            result[a - start:b - start] = fragment
            cursor = max(cursor, b)
        sequence = result.decode('ascii')
        return reverse_complement(sequence) if row['strand'] == '-' else sequence


def read_summary(path: str) -> Dict[str, List[Dict[str, str]]]:
    required = {
        "partition", "path_order", "type", "source_haplotype",
        "source_contig", "source_start", "source_end", "strand",
        "segments",
    }
    output: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    with open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: expected columns {sorted(required)}")
        for row in reader:
            output[row["partition"]].append(row)
    for partition, rows in output.items():
        rows.sort(key=lambda row: int(row["path_order"]))
        orders = [int(row["path_order"]) for row in rows]
        if len(orders) != len(set(orders)):
            raise ValueError(f"{path}: {partition} has duplicate path_order")
    if not output:
        raise ValueError(f"{path}: no local graph rows")
    return dict(output)


def read_partition_cache_package(
    path: str, known_partitions: Iterable[str],
) -> Dict[str, str]:
    """Read optional partition-keyed exact cache text records."""
    if not os.path.isfile(path):
        return {}
    known = set(known_partitions)
    caches: Dict[str, str] = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            try:
                record = json.loads(raw)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid JSON record: {error}"
                ) from error
            if not isinstance(record, dict):
                raise ValueError(f"{path}:{line_number}: expected JSON object")
            partition = record.get("partition")
            cache_text = record.get("cache_text")
            if not isinstance(partition, str) or partition not in known:
                raise ValueError(
                    f"{path}:{line_number}: unknown partition {partition!r}"
                )
            if partition in caches:
                raise ValueError(
                    f"{path}:{line_number}: duplicate partition {partition!r}"
                )
            if not isinstance(cache_text, str) or not cache_text:
                raise ValueError(
                    f"{path}:{line_number}: cache_text must be nonempty text"
                )
            try:
                parsed_cache = json.loads(cache_text)
            except json.JSONDecodeError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid embedded cache JSON: {error}"
                ) from error
            if not isinstance(parsed_cache, dict):
                raise ValueError(
                    f"{path}:{line_number}: embedded cache must be a JSON object"
                )
            caches[partition] = cache_text
    return caches


def bounded_map(function, items: Sequence, jobs: int):
    iterator = iter(items)
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
        pending = set()
        for _index in range(min(max(1, jobs), len(items))):
            try:
                pending.add(executor.submit(function, next(iterator)))
            except StopIteration:
                break
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                yield future.result()
                try:
                    pending.add(executor.submit(function, next(iterator)))
                except StopIteration:
                    pass


def blast_database_is_current(graph_fasta: str) -> bool:
    """Return whether a complete v4 BLAST database is as new as its FASTA."""
    prefix = graph_fasta + "_db"
    graph_mtime = os.stat(graph_fasta).st_mtime_ns
    v4_files = tuple(prefix + suffix for suffix in (".nhr", ".nin", ".nsq"))
    return not os.path.exists(prefix + ".ndb") and all(
        os.path.isfile(path) for path in v4_files
    ) and min(
        os.stat(path).st_mtime_ns for path in v4_files
    ) >= graph_mtime


def build_graph_blast_database(graph_fasta: str) -> None:
    """Run the bundled database builder beside this Python script."""
    if not os.path.isfile(RUNMAKEBLASTDB):
        raise FileNotFoundError(
            "required runmakeblastdb must be beside "
            f"{Path(__file__).resolve()}: {RUNMAKEBLASTDB}"
        )
    prefix = graph_fasta + "_db"
    for path in Path(graph_fasta).parent.glob(Path(prefix).name + ".n*"):
        try:
            path.unlink()
        except FileNotFoundError:
            pass
    try:
        subprocess.run(
            [
                "bash", RUNMAKEBLASTDB,
                "-in", graph_fasta,
                "-out", prefix,
                "-blastdb_version", "4",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception as error:
        # A failed builder must not leave a marker that --resume mistakes for
        # a database matching the newly reconstructed graph.
        for path in Path(graph_fasta).parent.glob(Path(prefix).name + ".n*"):
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        if isinstance(error, subprocess.CalledProcessError):
            detail = (error.stderr or "").strip()
            message = f"runmakeblastdb failed for {graph_fasta}"
            if detail:
                message += f":\n{detail}"
            raise RuntimeError(message) from error
        raise
    if not blast_database_is_current(graph_fasta):
        raise RuntimeError(
            f"runmakeblastdb did not create a complete database for {graph_fasta}"
        )


def reference_header(
    partition: str,
    row: Mapping[str, str],
    contig: str,
    contig_length: int,
) -> str:
    start, end = int(row["source_start"]), int(row["source_end"])
    strand = row["strand"]
    path_type = row["type"]
    record_id = public_path_id_from_values(
        partition, row["source_contig"], start, end,
    )
    if strand == "+":
        cigar = f">{contig}:{start}H{end - start}={contig_length - end}H"
    else:
        cigar = f"<{contig}:{contig_length - end}H{end - start}={start}H"
    coordinate = (
        f"{row['source_haplotype']}:{contig}:{start}-{end}{strand}"
    )
    return (
        f"{record_id} {coordinate} {contig}:{start}-{end}{strand} "
        f"{cigar} . {path_type}"
    )


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--package-dir", required=True,
        help="directory containing local_graphs.tsv and alternatives.fasta",
    )
    parser.add_argument("-r", "--reference-fasta", required=True)
    parser.add_argument(
        "--original_fasta", "--original-fasta", dest="original_fasta",
        help=(
            "original chromosome FASTA used to generate PARTITION.fasta from "
            "each partition BED; chr* names are resolved as NC_* aliases"
        ),
    )
    parser.add_argument(
        "-q", "--query-paths",
        help=argparse.SUPPRESS,  # Accept old invocations, but never read assemblies.
    )
    parser.add_argument(
        "--reference-haplotype", action="append", default=[],
        help=(
            "source_haplotype served by --reference-fasta; repeat for aliases. "
            "Otherwise inferred only for packages with one template haplotype."
        ),
    )
    parser.add_argument(
        "-o", "--output-dir", required=True,
        help="output root containing Graphs/, Graphs.list, Graphs.template.list[.bin], summary/ and optional references/ link",
    )
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    return args


def run(args: argparse.Namespace) -> None:
    package_dir = os.path.abspath(args.package_dir)
    output_dir = os.path.abspath(args.output_dir)
    graph_dir = Path(output_dir) / "Graphs"
    reference_path = os.path.abspath(args.reference_fasta)
    original_fasta_path = (
        os.path.abspath(args.original_fasta)
        if args.original_fasta else None
    )
    if args.query_paths:
        LOG.warning('Ignoring legacy --query-paths: reconstruction reads only the package and supplied reference FASTA')
    summary_path = os.path.join(package_dir, "local_graphs.tsv")
    alternatives_path = os.path.join(package_dir, "alternatives.fasta")
    caches_path = os.path.join(package_dir, PARTITION_CACHE_PACKAGE)
    for path in (
        summary_path, alternatives_path, reference_path,
        original_fasta_path,
    ):
        if path is None:
            continue
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    summary_link = Path(output_dir) / "summary"
    package_target = os.path.realpath(package_dir)
    if os.path.lexists(summary_link) and os.path.realpath(summary_link) != package_target:
        raise FileExistsError(
            f"{summary_link}: already exists and does not point to {package_dir}"
        )
    references_source = Path(package_dir) / "references"
    references_link = Path(output_dir) / "references"
    references_target = references_source.resolve() if references_source.is_dir() else None
    if (references_target is not None and os.path.lexists(references_link)
            and references_link.resolve() != references_target):
        raise FileExistsError(
            f"{references_link}: already exists and does not point to {references_source}"
        )

    LOG.info("Loading local graph summary: %s", summary_path)
    by_partition = read_summary(summary_path)
    partition_caches = read_partition_cache_package(
        caches_path, by_partition,
    )
    LOG.info(
        "Loaded %d optional partition cache record(s)",
        len(partition_caches),
    )
    LOG.info(
        "Loading all reference and alternative sequences into shared RAM",
    )
    reference = load_reference_ram(reference_path)
    if original_fasta_path is None:
        original_reference = None
    elif original_fasta_path == reference_path:
        original_reference = reference
        LOG.info("Reusing the in-RAM main reference for --original_fasta")
    else:
        LOG.info("Loading --original_fasta into shared RAM")
        original_reference = load_reference_ram(original_fasta_path)
    alternatives = load_fasta_ram(alternatives_path)
    LOG.info(
        "Loaded %d alternative records (%.2f GiB) into RAM",
        len(alternatives),
        sum(len(sequence) for _header, sequence in alternatives.values())
        / (1024 ** 3),
    )

    required_alternatives = {
        public_path_id_from_values(
            partition, row["source_contig"],
            row["source_start"], row["source_end"],
        )
        for partition, rows in by_partition.items()
        for row in rows
        if (
            row["type"] not in {"original", "reference"}
            and parsed_reference_slice(
                row, f"{summary_path}:{partition}:{row['path_order']}"
            ) is None
        )
    }
    optional_templates = {
        public_path_id_from_values(
            partition, row["source_contig"],
            row["source_start"], row["source_end"],
        )
        for partition, rows in by_partition.items()
        for row in rows if row['type'] in {'original', 'reference'}
    }
    allowed_alternatives = required_alternatives | optional_templates
    if (
        not required_alternatives.issubset(alternatives)
        or not set(alternatives).issubset(allowed_alternatives)
    ):
        missing = required_alternatives - set(alternatives)
        extra = set(alternatives) - allowed_alternatives
        raise ValueError(
            "alternatives.fasta/local_graphs.tsv mismatch; missing="
            + ",".join(sorted(missing)[:20])
            + "; extra=" + ",".join(sorted(extra)[:20])
        )

    packed_sources = PackedSourceIndex(by_partition, alternatives)
    unbacked_templates = [
        (partition, row)
        for partition, rows in by_partition.items()
        for row in [*selected_reference_rows(rows), *original_definition_rows(rows)]
        if not packed_sources.covering_paths(row)
    ]
    template_haplotypes = {
        row["source_haplotype"]
        for rows in by_partition.values()
        for row in [*selected_reference_rows(rows), *original_definition_rows(rows)]
    }
    main_haplotypes = set(args.reference_haplotype)
    if not main_haplotypes and len(template_haplotypes) == 1:
        # Legacy packages omit their main-reference originals. They are the
        # only source allowed to be absent from the packed sequence catalog.
        main_haplotypes = set(template_haplotypes)
    if not main_haplotypes and unbacked_templates:
        raise ValueError(
            'Use --reference-haplotype to name the supplied -r reference '
            '(for example CHM13_h1). This package has multiple template '
            'sources; absent source intervals must not be assigned to that '
            'reference by chromosome name alone. No cohort assemblies are read.'
        )
    missing = [
        f'{partition}:{row["source_haplotype"]}:{row["source_contig"]}:'
        f'{row["source_start"]}-{row["source_end"]}{row["strand"]}'
        for partition, row in unbacked_templates
        if row['source_haplotype'] not in main_haplotypes
        or canonical_contig(row['source_contig']) not in reference
    ]
    if missing:
        raise ValueError(
            'Incomplete summary package: original/template intervals are not '
            'fully covered by packed sequences or the designated reference: '
            + ','.join(missing[:10])
            + (f' (and {len(missing) - 10} more)' if len(missing) > 10 else '')
            + '. Reconstruction does not read cohort assemblies. The package '
            'producer must include these template bases in the summary. '
            'Use --reference-haplotype only to identify the supplied reference.'
        )
    LOG.info('Template sources: packed paths plus %s; no cohort assemblies',
             ','.join(sorted(main_haplotypes)) or 'verified reference_slice records')

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    input_mtime = max(
        os.stat(summary_path).st_mtime_ns,
        os.stat(alternatives_path).st_mtime_ns,
        os.stat(reference_path).st_mtime_ns,
        *(
            [os.stat(caches_path).st_mtime_ns]
            if os.path.isfile(caches_path) else []
        ),
        *(
            [os.stat(original_fasta_path).st_mtime_ns]
            if original_fasta_path is not None else []
        ),
    )

    def fetch_template(
        partition: str, row: Mapping[str, str],
    ) -> Tuple[str, int]:
        canonical = canonical_contig(row["source_contig"])
        haplotype = row["source_haplotype"]
        start, end = int(row["source_start"]), int(row["source_end"])
        if haplotype in main_haplotypes:
            try:
                chromosome = reference[canonical]
            except KeyError as error:
                raise ValueError(
                    f"{partition}: supplied main reference lacks {canonical!r} "
                    f"(table source was {row['source_contig']!r})"
                ) from error
            contig_length = len(chromosome)
            if start < 0 or end <= start or end > contig_length:
                raise ValueError(
                    f"{partition}: invalid main-reference interval "
                    f"{canonical}:{start}-{end} for length {contig_length}"
                )
            return chromosome[start:end], contig_length
        raise ValueError(f'{partition}: template {haplotype}:{canonical}:{start}-{end} is missing from the summary package')

    def fetch_reference_backed_path(
        partition: str, row: Mapping[str, str],
    ) -> Optional[Tuple[str, str]]:
        """Recreate a compacted path from the supplied main reference."""
        parsed = parsed_reference_slice(
            row, f"{summary_path}:{partition}:{row['path_order']}",
        )
        if parsed is None:
            return None
        table_contig, start, end, strand = parsed
        contig = canonical_contig(table_contig)
        try:
            chromosome = reference[contig]
        except KeyError as error:
            raise ValueError(
                f"{partition}: supplied reference lacks reference_slice "
                f"contig {table_contig!r}"
            ) from error
        if end > len(chromosome):
            raise ValueError(
                f"{partition}: reference_slice {table_contig}:{start}-{end} "
                f"exceeds reference length {len(chromosome)}"
            )
        sequence = chromosome[start:end]
        if strand == "-":
            sequence = reverse_complement(sequence)
        expected_digest = row.get("sequence_sha256", ".")
        observed_digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if (
            expected_digest in {None, "", "."}
            or observed_digest != expected_digest
        ):
            raise ValueError(
                f"{partition}: supplied reference does not match compacted "
                f"path {row['path_order']}"
            )
        description = row.get("fasta_header", ".")
        record_id = public_path_id_from_values(
            partition, row["source_contig"],
            row["source_start"], row["source_end"],
        )
        if description in {None, "", "."}:
            raise ValueError(
                f"{partition}: compacted path {record_id} lacks fasta_header"
            )
        if description.split(None, 1)[0] != record_id:
            raise ValueError(
                f"{partition}: compacted fasta_header ID does not match "
                f"{record_id}"
            )
        return description, sequence

    def fetch_packed_path(partition, row):
        record_id = public_path_id_from_values(
            partition, row['source_contig'], row['source_start'], row['source_end'],
        )
        if record_id in alternatives:
            return alternatives[record_id]
        packed = fetch_reference_backed_path(partition, row)
        if packed is None:
            raise ValueError(f'{partition}: missing packed path {record_id}')
        return packed

    def materialize_packed_interval(partition, row):
        sequence = packed_sources.sequence(row, fetch_packed_path)
        if sequence is None:
            return None
        record_id = public_path_id_from_values(
            partition, row['source_contig'], row['source_start'], row['source_end'],
        )
        coordinate = (f'{row["source_haplotype"]}:{row["source_contig"]}:'
                      f'{row["source_start"]}-{row["source_end"]}{row["strand"]}')
        # A source interval is provenance, not a claimed mapping to the main
        # reference. Preserve original bounds and do not invent chromosome length.
        return f'{record_id} {coordinate} . . . {row["type"]}', sequence

    def materialize_bed_template(
        partition: str, row: Mapping[str, str],
    ) -> Tuple[str, str]:
        """Return one annotated template record represented by the BED row."""
        record_id = public_path_id_from_values(
            partition, row['source_contig'], row['source_start'], row['source_end'],
        )
        if (record_id in alternatives
                or parsed_reference_slice(row, f'{partition}:original-template') is not None
                or row['source_haplotype'] not in main_haplotypes):
            packed = materialize_packed_interval(partition, row)
            if packed is not None:
                return packed
        if original_reference is None or row["source_haplotype"] not in main_haplotypes:
            sequence, contig_length = fetch_template(partition, row)
            if row["strand"] == "-":
                sequence = reverse_complement(sequence)
            return reference_header(
                partition, row, canonical_contig(row["source_contig"]), contig_length,
            ), sequence
        contig = canonical_contig(row["source_contig"])
        try:
            chromosome = original_reference[contig]
        except KeyError as error:
            raise ValueError(
                f"{partition}: --original_fasta lacks {contig!r}"
            ) from error
        start, end = int(row["source_start"]), int(row["source_end"])
        if start < 0 or end <= start or end > len(chromosome):
            raise ValueError(
                f"{partition}: BED template interval {contig}:{start}-{end} "
                f"exceeds --original_fasta length {len(chromosome)}"
            )
        sequence = chromosome[start:end]
        if row["strand"] == "-":
            sequence = reverse_complement(sequence)
        description = reference_header(
            partition, row, contig, len(chromosome),
        )
        return description, sequence

    def write_partition(item) -> Tuple[str, int, int, bool]:
        partition, rows = item
        directory = os.path.join(graph_dir, partition)
        bed_path = os.path.join(directory, f"{partition}.bed")
        graph_path = os.path.join(directory, f"{partition}.FA")
        cache_path = os.path.join(directory, f"{partition}cache.json")
        cache_text = partition_caches.get(partition)
        original_graph_path = os.path.join(directory, f"{partition}.fasta")
        legacy_graph_path = os.path.join(directory, "graph.FA")
        legacy_original_graph_path = os.path.join(directory, "graph.fasta")

        def remove_legacy_names() -> None:
            for legacy in (legacy_graph_path, legacy_original_graph_path):
                try:
                    os.remove(legacy)
                except FileNotFoundError:
                    pass

        if (
            args.resume
            and os.path.isfile(bed_path) and os.path.getsize(bed_path) > 0
            and os.path.isfile(graph_path) and os.path.getsize(graph_path) > 0
            and os.stat(bed_path).st_mtime_ns >= input_mtime
            and os.stat(graph_path).st_mtime_ns >= input_mtime
            and blast_database_is_current(graph_path)
            and (
                cache_text is None
                or (
                    os.path.isfile(cache_path)
                    and os.stat(cache_path).st_mtime_ns >= input_mtime
                    and Path(cache_path).read_text(encoding="utf-8")
                    == cache_text
                )
            )
            and os.path.isfile(original_graph_path)
            and os.path.getsize(original_graph_path) > 0
            and os.stat(original_graph_path).st_mtime_ns >= input_mtime
        ):
            remove_legacy_names()
            return partition, 0, 0, True
        has_refined = any(row["type"] == "reference" for row in rows)
        graph_template_rows = selected_reference_rows(rows)
        bed_rows = original_definition_rows(rows)
        graph_rows = [
            row for row in rows
            if row["type"] != "original" or not has_refined
        ]
        if not graph_template_rows:
            raise ValueError(f"{partition}: no reference template row")
        if not bed_rows:
            raise ValueError(f"{partition}: no original partition definition row")
        Path(directory).mkdir(parents=True, exist_ok=True)
        bed_tmp = bed_path + ".tmp"
        graph_tmp = graph_path + ".tmp"
        original_graph_tmp = original_graph_path + ".tmp"
        cache_tmp = cache_path + ".tmp"
        sequence_count = 0
        base_count = 0
        try:
            with open(bed_tmp, "wt") as bed:
                for row in bed_rows:
                    contig = canonical_contig(row["source_contig"])
                    record_id = public_path_id_from_values(
                        partition, row["source_contig"],
                        row["source_start"], row["source_end"],
                    )
                    bed.write(
                        f"{contig}\t{row['source_start']}\t{row['source_end']}\t"
                        f"{partition}\t1000\t{row['strand']}\t"
                        f"{row['source_haplotype']}\t{record_id}\t"
                        f"{row['source_start']}\t{row['source_end']}\t"
                        "source_coordinates\n"
                    )
            with open(graph_tmp, "wt") as graph:
                for row in graph_rows:
                    record_id = public_path_id_from_values(
                        partition, row["source_contig"],
                        row["source_start"], row["source_end"],
                    )
                    if record_id in alternatives:
                        description, sequence = alternatives[record_id]
                    else:
                        compacted = fetch_reference_backed_path(partition, row)
                        if compacted is not None:
                            description, sequence = compacted
                        else:
                            packed = materialize_packed_interval(partition, row)
                            if packed is not None:
                                description, sequence = packed
                                graph.write(f">{description}\n{wrap_fasta(sequence)}\n")
                                sequence_count += 1
                                base_count += len(sequence)
                                continue
                            contig = canonical_contig(row["source_contig"])
                            start, end = (
                                int(row["source_start"]),
                                int(row["source_end"]),
                            )
                            sequence, contig_length = fetch_template(
                                partition, row,
                            )
                            if row["strand"] == "-":
                                sequence = reverse_complement(sequence)
                            description = reference_header(
                                partition, row, contig, contig_length,
                            )
                    expected_length = (
                        int(row["source_end"]) - int(row["source_start"])
                    )
                    if len(sequence) != expected_length:
                        raise ValueError(
                            f"{partition}: graph path {row['path_order']} has "
                            f"sequence length {len(sequence)}, but its source "
                            f"interval is {row['source_start']}-"
                            f"{row['source_end']} ({expected_length} bases)"
                        )
                    graph.write(f">{description}\n{wrap_fasta(sequence)}\n")
                    sequence_count += 1
                    base_count += len(sequence)
            with open(original_graph_tmp, "wt") as original_graph:
                for row in bed_rows:
                    description, sequence = materialize_bed_template(
                        partition, row,
                    )
                    original_graph.write(
                        f">{description}\n{wrap_fasta(sequence)}\n"
                    )
            os.replace(bed_tmp, bed_path)
            os.replace(graph_tmp, graph_path)
            os.replace(original_graph_tmp, original_graph_path)
            build_graph_blast_database(graph_path)
            if cache_text is not None:
                with open(cache_tmp, "wt", encoding="utf-8") as cache_output:
                    cache_output.write(cache_text)
                os.replace(cache_tmp, cache_path)
            remove_legacy_names()
        finally:
            for temporary in (
                bed_tmp, graph_tmp, original_graph_tmp, cache_tmp,
            ):
                try:
                    os.remove(temporary)
                except FileNotFoundError:
                    pass
        return partition, sequence_count, base_count, False

    items = sorted(by_partition.items())
    LOG.info(
        "Reconstructing %d independent local graph folders with %d workers",
        len(items), args.jobs,
    )
    graphs = bases = resumed = 0
    progress_step = max(1, len(items) // 20)
    for completed, (
        _partition, count, base_count, was_resumed,
    ) in enumerate(bounded_map(write_partition, items, args.jobs), 1):
        graphs += count
        bases += base_count
        resumed += was_resumed
        if completed == len(items) or completed % progress_step == 0:
            LOG.info(
                "Folder reconstruction progress: %d/%d partitions",
                completed, len(items),
            )
    LOG.info(
        "Reconstructed %d local graph folders and their BLAST databases with "
        "%d FASTA paths and %d bases; resumed %d folders",
        len(items), graphs, bases, resumed,
    )
    valid_graph_fastas = [
        graph_dir / partition / f"{partition}.FA"
        for partition, _rows in items
    ]
    # KmerSearcher builds a cohort-wide target hash. Release the multi-GiB
    # reconstruction sequence stores before starting that independent process.
    del alternatives
    del reference
    del original_reference
    gc.collect()
    LOG.info("Preparing graph alignment list and template-only hotspot index for %d partitions", len(items))
    listing, cache, _rebuilt = ensure_graph_list_cache(
        graph_dir,
        jobs=args.jobs,
        graph_fastas=valid_graph_fastas,
        script_dir=Path(__file__).resolve().parent,
    )
    LOG.info(
        "Generated graph alignment list and KmerSearcher template cache: %s, %s",
        listing, cache,
    )
    if not os.path.lexists(summary_link):
        summary_link.symlink_to(
            os.path.relpath(package_target, os.path.realpath(output_dir)),
            target_is_directory=True,
        )
    LOG.info("Summary package available at %s", summary_link)
    if references_target is not None:
        if not os.path.lexists(references_link):
            references_link.symlink_to(
                os.path.relpath(references_target, os.path.realpath(output_dir)),
                target_is_directory=True,
            )
        LOG.info("Reference caches available at %s", references_link)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    try:
        run(args)
    except Exception as error:
        LOG.error("%s", error)
        if args.log_level == "DEBUG":
            raise
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
