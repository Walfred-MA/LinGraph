#!/usr/bin/env python3
"""Single-sample graph-reference CIGAR to VCF conversion.

This entry point reuses graphreftovcf.py's parsing, clustering, and VCF
formatting logic while changing sequence access for a single assembly:

* the complete reference and assembly FASTAs are loaded into RAM once;
* input-row allele names become zero-copy views onto assembly contig intervals;
* input CIGAR rows are parsed in parallel by fork workers sharing those bases
  copy-on-write.

The input may use the individual-block seven-column layout emitted by
graphcigartoref_persample.py: query, reference, query coordinate, reference
coordinate, query alignment region, reference alignment region, and CIGAR.
The core parser continues to accept the older legacy order as well.

Use ``--fasta-query`` for the original assembly FASTA.  It replaces the old
per-allele records FASTA normally passed with ``-s/--seq``.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import sys
from typing import Dict, List, Optional, Tuple

import graphreftovcf as core


AlleleView = Tuple[str, int, int, str]
_FASTA_CACHE: Dict[str, Dict[str, str]] = {}
_VIEWS_BY_FASTA: Dict[str, Dict[str, AlleleView]] = {}
_ORIGINAL_READ_FASTA_METADATA = core.read_fasta_metadata


def _load_fasta(path: str) -> Dict[str, str]:
    absolute = os.path.abspath(os.path.expanduser(path))
    cached = _FASTA_CACHE.get(absolute)
    if cached is not None:
        return cached

    sequences: Dict[str, str] = {}
    name: Optional[str] = None
    chunks: List[str] = []
    with core.open_text(absolute) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith("#"):
                continue
            if raw.startswith(">"):
                if name is not None:
                    if name in sequences:
                        raise ValueError(
                            f"{absolute}:{line_number}: duplicate FASTA record {name!r}"
                        )
                    sequences[name] = "".join(chunks)
                fields = raw[1:].strip().split()
                if not fields:
                    raise ValueError(f"{absolute}:{line_number}: empty FASTA header")
                name = fields[0]
                chunks = []
            elif raw.strip():
                if name is None:
                    raise ValueError(
                        f"{absolute}:{line_number}: sequence before FASTA header"
                    )
                chunks.append(raw.strip())
    if name is not None:
        if name in sequences:
            raise ValueError(f"{absolute}: duplicate FASTA record {name!r}")
        sequences[name] = "".join(chunks)
    if not sequences:
        raise ValueError(f"no FASTA records found: {absolute}")
    _FASTA_CACHE[absolute] = sequences
    print(
        f"[graphreftovcf_persample] loaded {absolute} into RAM: "
        f"{len(sequences)} record(s), {sum(map(len, sequences.values()))} bases",
        file=sys.stderr,
    )
    return sequences


class InMemoryFastaReader:
    """IndexedFastaReader-compatible full-memory FASTA plus virtual views."""

    def __init__(self, fasta_path: str, index_path: Optional[str] = None):
        del index_path
        self.fasta_path = os.path.abspath(os.path.expanduser(fasta_path))
        self.index_path = "<in-memory>"
        self.sequences = _load_fasta(self.fasta_path)
        self.views = _VIEWS_BY_FASTA.get(self.fasta_path, {})
        self.index: Dict[str, Tuple[int, int, int, int]] = {
            name: (len(sequence), 0, max(1, len(sequence)), max(1, len(sequence)))
            for name, sequence in self.sequences.items()
        }
        for name, (_contig, start, end, _strand) in self.views.items():
            length = max(0, end - start)
            self.index[name] = (length, 0, max(1, length), max(1, length))

    def __contains__(self, name: str) -> bool:
        return name in self.index

    def length(self, name: str) -> int:
        return self.index[name][0]

    def fetch(self, name: str, start: int, end: int, strand: str = "+") -> str:
        if name in self.views:
            contig, view_start, view_end, view_strand = self.views[name]
            source = self.sequences.get(contig)
            if source is None:
                raise KeyError(
                    f"assembly contig {contig!r} for allele view {name!r} "
                    f"is absent from {self.fasta_path}"
                )
            view_length = view_end - view_start
            start = max(0, min(int(start), view_length))
            end = max(start, min(int(end), view_length))
            if view_strand == "-":
                sequence = core.revcomp(
                    source[view_end - end:view_end - start]
                )
            else:
                sequence = source[view_start + start:view_start + end]
        else:
            source = self.sequences.get(name)
            if source is None:
                raise KeyError(f"{name!r} not found in {self.fasta_path}")
            start = max(0, min(int(start), len(source)))
            end = max(start, min(int(end), len(source)))
            sequence = source[start:end]
        return core.revcomp(sequence) if strand == "-" else sequence


def _read_fasta_metadata_with_views(
    path: str,
    require_coord: bool,
    max_impute: int,
    attach_reader: bool = True,
):
    records, order, reader = _ORIGINAL_READ_FASTA_METADATA(
        path,
        require_coord=require_coord,
        max_impute=max_impute,
        attach_reader=attach_reader,
    )
    if reader is None or not isinstance(reader, InMemoryFastaReader):
        return records, order, reader
    for name, (contig, start, end, strand) in reader.views.items():
        if contig not in reader.sequences:
            raise KeyError(
                f"assembly contig {contig!r} for input query {name!r} is "
                f"absent from {reader.fasta_path}"
            )
        if start < 0 or end < start or end > len(reader.sequences[contig]):
            raise ValueError(
                f"input query {name!r} view {contig}:{start}-{end}{strand} "
                f"is outside contig length {len(reader.sequences[contig])}"
            )
        record = core.SeqRecord(
            name=name,
            allelename=name,
            contig=contig,
            start=start,
            end=end,
            strand=strand,
            seqlen=end - start,
            reader=reader,
            fasta_name=name,
            source="in-memory-assembly-view",
        )
        core.add_alias(records, record)
        if name not in order:
            order.append(name)
    return records, order, reader


def _option_value(arguments: List[str], names: Tuple[str, ...]) -> Optional[str]:
    for index, token in enumerate(arguments):
        for name in names:
            if token == name:
                if index + 1 >= len(arguments):
                    raise ValueError(f"{name} requires a value")
                return arguments[index + 1]
            if token.startswith(name + "="):
                return token.split("=", 1)[1]
    return None


def _build_allele_views(input_path: str) -> Dict[str, AlleleView]:
    views: Dict[str, AlleleView] = {}
    with core.open_text(input_path) as handle:
        for line_number, raw in enumerate(handle, 1):
            row = core.parse_graphcigartoref_row(raw, line_number)
            if row is None:
                continue
            coordinate = core.parse_coord_token(row.get("query_coord", ""))
            if coordinate is None:
                raise ValueError(
                    f"{input_path}:{line_number}: invalid query coordinate "
                    f"{row.get('query_coord', '')!r}"
                )
            name = row["query"]
            old = views.get(name)
            if old is not None and old != coordinate:
                if old[0] != coordinate[0] or old[3] != coordinate[3]:
                    raise ValueError(
                        f"input query {name!r} has conflicting assembly coordinates: "
                        f"{old!r} and {coordinate!r}"
                    )
                # Different final edge crops of one query share a covering
                # assembly view; parse_row_events selects each row's offset.
                coordinate = (old[0], min(old[1], coordinate[1]), max(old[2], coordinate[2]), old[3])
            views[name] = coordinate
    if not views:
        raise ValueError(f"no graph-reference CIGAR rows found: {input_path}")
    return views


def _fork_context():
    try:
        return mp.get_context("fork")
    except (ValueError, RuntimeError) as error:
        raise RuntimeError(
            "graphreftovcf_persample.py requires multiprocessing fork to "
            "share in-memory FASTA sequences"
        ) from error


def main() -> int:
    arguments = sys.argv[1:]
    if "-h" in arguments or "--help" in arguments:
        # graphreftovcf.py owns the complete option set.  It also exposes
        # --fasta-query as an alias for -s/--seq.
        core.main()
        return 0

    if "--per-query-columns" in arguments:
        raise ValueError(
            "--per-query-columns is incompatible with "
            "graphreftovcf_persample.py: graph/query records are partitions "
            "of one biological sample, not separate VCF samples; remove this "
            "option and use --columns SAMPLE_h1"
        )

    input_path = _option_value(arguments, ("--input", "-i"))
    query_fasta = _option_value(arguments, ("--fasta-query", "--seq", "-s"))
    reference_fasta = _option_value(arguments, ("--ref", "-r"))
    if not input_path:
        raise ValueError("-i/--input is required")
    if not query_fasta:
        raise ValueError("--fasta-query is required for the per-sample converter")
    if not reference_fasta:
        raise ValueError("-r/--ref is required")

    input_path = core.normalize_cli_path(input_path, "--input")
    query_fasta = core.normalize_cli_path(query_fasta, "--fasta-query")
    reference_fasta = core.normalize_cli_path(reference_fasta, "--ref")
    input_path = os.path.abspath(os.path.expanduser(input_path))
    query_fasta = os.path.abspath(os.path.expanduser(query_fasta))
    reference_fasta = os.path.abspath(os.path.expanduser(reference_fasta))
    for path in (input_path, query_fasta, reference_fasta):
        if not os.path.isfile(path):
            raise FileNotFoundError(path)

    views = _build_allele_views(input_path)
    _VIEWS_BY_FASTA[query_fasta] = views
    print(
        f"[graphreftovcf_persample] indexed {len(views)} input allele view(s) "
        "on the assembly",
        file=sys.stderr,
    )

    # A per-sample input should produce one biological VCF column even though
    # it contains tens of thousands of graph/query records. Infer that column
    # from graph_sample_haplotype_lineindex names unless the caller supplied it.
    requested_columns = _option_value(arguments, ("--columns", "-c"))
    inferred_samples = sorted({
        haplotype
        for name in views
        for haplotype in (core.get_haplotype(name),)
        if haplotype.endswith(("_h1", "_h2"))
    })
    if requested_columns is None:
        if len(inferred_samples) != 1:
            raise ValueError(
                "could not infer exactly one sample/haplotype from input "
                f"query names; found {inferred_samples[:20]!r}; provide "
                "--columns SAMPLE_h1"
            )
        arguments.extend(("--columns", inferred_samples[0]))
        sys.argv = [sys.argv[0], *arguments]
        print(
            f"[graphreftovcf_persample] inferred VCF sample column "
            f"{inferred_samples[0]!r}",
            file=sys.stderr,
        )

    if not any(
        option in arguments
        for option in (
            "--validate-cross-graph-svs",
            "--merge-cross-row-deletions",
        )
    ):
        arguments.append("--validate-cross-graph-svs")
        sys.argv = [sys.argv[0], *arguments]

    # Install the in-memory/shared implementations before graphreftovcf.main()
    # creates any readers or worker pools.
    core.IndexedFastaReader = InMemoryFastaReader
    core.read_fasta_metadata = _read_fasta_metadata_with_views
    core.mp_context = _fork_context
    core.main()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"graphreftovcf_persample.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
