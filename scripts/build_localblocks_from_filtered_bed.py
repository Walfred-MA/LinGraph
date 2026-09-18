#!/usr/bin/env python3
"""Build local block folders from one reference-to-block filtered BED.

The input is the single 16-column
``main_chroms_withnovel_blocks.fasta_aligned.filtered.bed`` produced after
mapping the reference assembly to all block targets. BED column 4 is the local
graph/folder name. The first entry in query_paths.txt is the reference assembly
whose query coordinates appear in BED columns 1-3.
"""
from __future__ import annotations

import argparse
import dataclasses
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Tuple

from build_local_graphs import reverse_complement
from minsetref import read_assembly_list
from minsetref_core import IndexedFasta, sanitize_id, wrap_fasta
from minsetref_segments import count_unmasked


LOG = logging.getLogger("build_localblocks_from_filtered_bed")
SOURCE_COORD_RE = re.compile(
    r"^(?:source=)?(.*):([^:]+):(\d+)-(\d+)([+-])?$"
)
GRAPH_PREFIX_RE = re.compile(r"^g\d+_(.+)$")
COORDINATE_NAME_RE = re.compile(r"^(.+)_([0-9]+)_([0-9]+)$")
MERGED_GRAPH_NAME_RE = re.compile(r"^merged_((?:g[0-9]+)+)$")
MAX_MERGED_GRAPH_IDS = 10


@dataclasses.dataclass(frozen=True)
class SourceInterval:
    sample: str
    contig: str
    start: int
    end: int
    strand: str = "+"


@dataclasses.dataclass(frozen=True)
class MappingRow:
    query_contig: str
    query_start: int
    query_end: int
    target: str
    source_target: str
    score: int
    strand: str
    query_sample_field: str
    target_start: int
    target_end: int


@dataclasses.dataclass(frozen=True)
class Template:
    sample: str
    contig: str
    # FASTA/source interval including extension anchors.
    start: int
    end: int
    # Unanchored interval represented by the reference BED.
    core_start: int
    core_end: int
    strand: str
    sequence: str
    origin: str
    fixed: bool = False


class AssemblyReaderCache:
    """Small LRU of source-assembly FASTA readers used for anchor fetching."""

    def __init__(self, paths: Dict[str, str], max_open: int = 32) -> None:
        self.paths = dict(paths)
        self.max_open = max(1, int(max_open))
        self.catalog: Dict[str, IndexedFasta] = {}
        self.readers: "OrderedDict[str, IndexedFasta]" = OrderedDict()

    def get(self, sample: str) -> Optional[IndexedFasta]:
        path = self.paths.get(sample)
        if path is None:
            return None
        reader = self.readers.pop(sample, None)
        if reader is None:
            reader = self.catalog.get(sample)
            if reader is None:
                reader = IndexedFasta(path)
                self.catalog[sample] = reader
            else:
                reader.reopen()
        self.readers[sample] = reader
        while len(self.readers) > self.max_open:
            _old_sample, old_reader = self.readers.popitem(last=False)
            old_reader.close()
        return reader

    def close(self) -> None:
        for reader in self.readers.values():
            reader.close()
        self.readers.clear()
        self.catalog.clear()

    def __enter__(self) -> "AssemblyReaderCache":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def _split_haplotype_contig(source: str, context: str) -> Tuple[str, str]:
    fields = source.split("_", 2)
    if len(fields) != 3 or not all(fields):
        raise ValueError(
            f"{context}: expected SAMPLE_HAPLOTYPE_CONTIG, got {source!r}"
        )
    return f"{fields[0]}_{fields[1]}", fields[2]


def _parse_coordinate_name(name: str) -> Optional[SourceInterval]:
    candidate = name
    prefixed = GRAPH_PREFIX_RE.fullmatch(candidate)
    if prefixed is not None:
        candidate = prefixed.group(1)
    match = COORDINATE_NAME_RE.fullmatch(candidate)
    if match is None:
        return None
    source, start_text, end_text = match.groups()
    start, end = int(start_text), int(end_text)
    if end <= start or source.count("_") < 2:
        return None
    sample, contig = _split_haplotype_contig(source, name)
    return SourceInterval(sample, contig, start, end, "+")


def parse_target_origin(target: str, description: str) -> Optional[SourceInterval]:
    """Recover a block's source interval from its ID or header annotation."""
    for token in description.split()[1:]:
        match = SOURCE_COORD_RE.fullmatch(token)
        if match is None:
            continue
        sample, contig, start_text, end_text, strand = match.groups()
        start, end = int(start_text), int(end_text)
        if end <= start:
            raise ValueError(
                f"target {target!r}: invalid source coordinate {token!r}"
            )
        return SourceInterval(sample, contig, start, end, strand or "+")
    return _parse_coordinate_name(target)


def read_fasta_descriptions(path: str) -> Dict[str, str]:
    descriptions: Dict[str, str] = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.startswith(">"):
                continue
            description = raw[1:].strip()
            if not description:
                raise ValueError(f"{path}:{line_number}: empty FASTA header")
            name = description.split()[0]
            if name in descriptions:
                raise ValueError(f"{path}:{line_number}: duplicate FASTA ID {name!r}")
            descriptions[name] = description
    if not descriptions:
        raise ValueError(f"{path}: no FASTA records")
    return descriptions


def _safe_graph_name(name: str, context: str) -> str:
    if (
        not name
        or name in {".", ".."}
        or os.path.basename(name) != name
        or "\x00" in name
    ):
        raise ValueError(f"{context}: unsafe target/folder name {name!r}")
    return name


def local_merged_graph_name(name: str) -> str:
    """Convert a merger group to a bounded local-graph folder/file name."""
    match = MERGED_GRAPH_NAME_RE.fullmatch(name)
    if match is None:
        return name
    identifiers = re.findall(r"g[0-9]+", match.group(1))
    compact = "merged" + "".join(identifiers[:MAX_MERGED_GRAPH_IDS])
    return compact


def sort_filtered_bed(input_bed: str, output_bed: str, workdir: str) -> None:
    """Externally sort one BED by target, then query coordinates."""
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    temporary = output_bed + ".tmp"
    command = [
        "sort", "-t", "\t", "-T", workdir,
        "-k4,4", "-k1,1", "-k2,2n", "-k3,3n", "-k6,6",
        input_bed,
    ]
    try:
        with open(temporary, "wb") as out:
            subprocess.run(command, check=True, stdout=out, env=environment)
        os.replace(temporary, output_bed)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def read_merged_groups(path: str) -> Dict[str, str]:
    """Read KmerPartitionMerger's original-to-final partition mapping."""
    groups: Dict[str, str] = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if fields == ["partition", "merged_partition"]:
                continue
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected PARTITION<TAB>MERGED_PARTITION"
                )
            partition = _safe_graph_name(fields[0], f"{path}:{line_number}")
            merged = local_merged_graph_name(
                _safe_graph_name(fields[1], f"{path}:{line_number}")
            )
            if partition in groups:
                raise ValueError(
                    f"{path}:{line_number}: duplicate partition {partition!r}"
                )
            groups[partition] = merged
    if not groups:
        raise ValueError(f"{path}: no partition mappings")
    return groups


def read_merged_bed(
    path: str,
) -> Dict[str, List[Tuple[str, int, int]]]:
    """Read KmerPartitionMerger BED4 intervals grouped by final partition."""
    intervals: DefaultDict[str, List[Tuple[str, int, int]]] = defaultdict(list)
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected BED4+")
            graph = local_merged_graph_name(
                _safe_graph_name(fields[3], f"{path}:{line_number}")
            )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid BED coordinates"
                ) from error
            if not fields[0] or start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid BED interval")
            intervals[graph].append((fields[0], start, end))
    if not intervals:
        raise ValueError(f"{path}: no merged BED intervals")
    return dict(intervals)


def remap_filtered_bed(
    input_bed: str,
    output_bed: str,
    groups: Dict[str, str],
) -> None:
    """Replace BED column 4 with its final group and retain original target."""
    seen = set()
    with open(input_bed, "rt") as source, open(output_bed, "wt") as output:
        for line_number, raw in enumerate(source, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 16:
                raise ValueError(
                    f"{input_bed}:{line_number}: expected 16-column filtered BED"
                )
            original = fields[3]
            try:
                fields[3] = groups[original]
            except KeyError as error:
                raise ValueError(
                    f"{input_bed}:{line_number}: target {original!r} is absent "
                    "from the KmerPartitionMerger groups TSV"
                ) from error
            fields.append(original)
            output.write("\t".join(fields) + "\n")
            seen.add(original)
    unused = set(groups) - seen
    if unused:
        final_groups = {groups[partition] for partition in unused}
        LOG.warning(
            "%d original target partitions have no surviving row in the "
            "filtered BED (across %d final merged groups); retaining their "
            "sequences from -t through the merger groups TSV",
            len(unused), len(final_groups),
        )


def parse_mapping_row(raw: str, path: str, line_number: int) -> MappingRow:
    fields = raw.rstrip("\n").split("\t")
    if len(fields) < 16:
        raise ValueError(
            f"{path}:{line_number}: expected the mapper's 16-column filtered BED"
        )
    target = _safe_graph_name(fields[3], f"{path}:{line_number}")
    source_target = _safe_graph_name(
        fields[16] if len(fields) > 16 else target,
        f"{path}:{line_number}",
    )
    strand = fields[5]
    if strand not in {"+", "-"}:
        raise ValueError(f"{path}:{line_number}: invalid strand {strand!r}")
    try:
        query_start, query_end = int(fields[1]), int(fields[2])
        score = int(fields[4])
        target_start, target_end = int(fields[7]), int(fields[8])
    except ValueError as error:
        raise ValueError(f"{path}:{line_number}: invalid numeric field") from error
    if query_start < 0 or query_end <= query_start:
        raise ValueError(f"{path}:{line_number}: invalid query interval")
    if target_start < 0 or target_end <= target_start:
        raise ValueError(f"{path}:{line_number}: invalid target interval")
    return MappingRow(
        fields[0], query_start, query_end, target, source_target, score, strand,
        fields[6], target_start, target_end,
    )


def iter_target_groups(path: str) -> Iterable[Tuple[str, List[MappingRow]]]:
    current: Optional[str] = None
    rows: List[MappingRow] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            row = parse_mapping_row(raw, path, line_number)
            if current is None:
                current = row.target
            elif row.target != current:
                yield current, rows
                current = row.target
                rows = []
            rows.append(row)
    if current is not None:
        yield current, rows


def _fetch_template(
    reader: IndexedFasta,
    sample: str,
    contig: str,
    core_start: int,
    core_end: int,
    sequence_start: int,
    sequence_end: int,
    strand: str,
    origin: str,
) -> Template:
    if contig not in reader.index:
        raise KeyError(f"{contig!r} is absent from reference FASTA {reader.path}")
    length = reader.length(contig)
    if (
        sequence_start < 0
        or sequence_end <= sequence_start
        or sequence_end > length
        or core_start < sequence_start
        or core_end > sequence_end
        or core_end <= core_start
    ):
        raise ValueError(
            f"invalid anchored/core interval for {sample}:{contig}: "
            f"sequence={sequence_start}-{sequence_end}, "
            f"core={core_start}-{core_end}, contig_length={length}"
        )
    sequence = reader.fetch(contig, sequence_start, sequence_end)
    if strand == "-":
        sequence = reverse_complement(sequence)
    return Template(
        sample, contig, sequence_start, sequence_end,
        core_start, core_end, strand, sequence, origin,
    )


def merge_core_intervals_with_anchors(
    intervals: Iterable[Tuple[int, int]],
    contig_length: int,
    anchor: int,
) -> List[Tuple[int, int, int, int]]:
    """Merge cores whose anchored spans overlap, then add clipped anchors."""
    merge_distance = max(1, 2 * anchor)
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if start < 0 or end <= start or end > contig_length:
            raise ValueError(
                f"invalid source interval {start}-{end} for contig length "
                f"{contig_length}"
            )
        if not merged or start - merged[-1][1] >= merge_distance:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [
        (
            core_start,
            core_end,
            max(0, core_start - anchor),
            min(contig_length, core_end + anchor),
        )
        for core_start, core_end in merged
    ]


def build_templates(
    graph: str,
    rows: Sequence[MappingRow],
    reference_sample: str,
    source_readers: AssemblyReaderCache,
    target_reader: IndexedFasta,
    target_descriptions: Dict[str, str],
    anchor: int,
    member_targets: Optional[Sequence[str]] = None,
    merged_query_intervals: Optional[Sequence[Tuple[str, int, int]]] = None,
    target_origins_only: bool = False,
) -> List[Template]:
    """Create one or more initial paths for a target block.

    Core intervals are grouped by source assembly/contig/strand, merged when
    their gap is less than twice ``anchor``, and then extended by ``anchor`` on
    each side for FASTA extraction. BED output retains the merged core bounds.
    """
    source_intervals: DefaultDict[
        Tuple[str, str, str], List[Tuple[int, int]]
    ] = defaultdict(list)
    reference_reader = source_readers.get(reference_sample)
    if reference_reader is None:
        raise KeyError(
            f"{graph}: reference sample {reference_sample!r} is absent from "
            "query_paths.txt"
        )

    def add_reference_interval(
        contig: str, start: int, end: int, strand: str,
    ) -> None:
        # Mapper/grouped BED query names come from the combined sequence
        # catalog, not necessarily from the primary reference.  Alternative
        # and discovered-novel records are already represented by their target
        # FASTA members below; never reinterpret their catalog IDs as CHM13
        # contigs.
        if contig in reference_reader.index:
            source_intervals[(reference_sample, contig, strand)].append(
                (start, end)
            )

    # A pre-grouped BED is extracted directly from the combined
    # reference/alternative/novel catalog.  In that mode the target FASTA
    # headers are the authoritative source mapping.  Its synthetic mapping
    # BED and grouped-merge BED retain catalog record names, which must not be
    # mistaken for contigs in the primary reference assembly.
    if target_origins_only:
        pass
    elif merged_query_intervals is None:
        for row in rows:
            add_reference_interval(
                row.query_contig, row.query_start, row.query_end, row.strand,
            )
    else:
        for contig, start, end in merged_query_intervals:
            add_reference_interval(contig, start, end, "+")

    templates: List[Template] = []
    if member_targets is None:
        member_targets = sorted({row.source_target for row in rows})
    for member in sorted(set(member_targets)):
        origin = parse_target_origin(member, target_descriptions[member])
        target_length = target_reader.length(member)
        if target_length <= 0:
            raise ValueError(f"target {member!r} has an empty sequence")
        from fixed_alternatives import is_fixed
        if is_fixed(target_descriptions[member]):
            # Imported bytes and orientation are authoritative even when the
            # source sample is in the cohort or anchors were requested.
            sequence = target_reader.sequence(member)
            if origin is not None and origin.end-origin.start != len(sequence):
                raise ValueError(f'{member}: fixed alternative source length mismatch')
            templates.append(Template(
                origin.sample if origin else 'target',
                origin.contig if origin else member,
                origin.start if origin else 0, origin.end if origin else len(sequence),
                origin.start if origin else 0, origin.end if origin else len(sequence),
                origin.strand if origin else '+', sequence, 'imported_alternative', True,
            ))
            continue
        if origin is None:
            sequence = target_reader.sequence(member)
            templates.append(Template(
                "target", member, 0, len(sequence), 0, len(sequence),
                "+", sequence, "target_block",
            ))
            continue
        if origin.end - origin.start != target_length:
            raise ValueError(
                f"target {member!r}: encoded source length "
                f"{origin.end - origin.start} differs from FASTA length "
                f"{target_length}"
            )
        if anchor == 0 and origin.sample != reference_sample:
            sequence = target_reader.sequence(member)
            if origin.strand == "-":
                sequence = reverse_complement(sequence)
            templates.append(Template(
                origin.sample, origin.contig, origin.start, origin.end,
                origin.start, origin.end, origin.strand, sequence,
                "target_block",
            ))
            continue
        source_intervals[(
            origin.sample, origin.contig, origin.strand,
        )].append((origin.start, origin.end))

    for (sample, contig, strand), intervals in sorted(source_intervals.items()):
        reader = source_readers.get(sample)
        if reader is None:
            raise KeyError(
                f"{graph}: source sample {sample!r} is absent from query_paths.txt; "
                "cannot fetch extension anchors"
            )
        if contig not in reader.index:
            raise KeyError(
                f"{graph}: source contig {contig!r} is absent from "
                f"{reader.path} for {sample}"
            )
        for core_start, core_end, sequence_start, sequence_end in (
            merge_core_intervals_with_anchors(
                intervals, reader.length(contig), anchor,
            )
        ):
            templates.append(_fetch_template(
                reader, sample, contig, core_start, core_end,
                sequence_start, sequence_end, strand,
                "lifted_reference" if sample == reference_sample
                else "target_block",
            ))
    templates.sort(key=lambda row: (
        row.sample != reference_sample, row.contig, row.start, row.end, row.strand,
    ))
    return templates


def graph_id_prefix(graph: str) -> str:
    legacy = re.match(r"^(g[0-9]+)_", graph)
    if legacy is not None:
        return legacy.group(1)
    return graph


def materialize_one_graph(
    graph: str,
    templates: Sequence[Template],
    output_root: str,
) -> Tuple[int, int]:
    directory = os.path.join(output_root, graph)
    Path(directory).mkdir(parents=True, exist_ok=True)
    fasta_path = os.path.join(directory, f"{graph}.fasta")
    bed_path = os.path.join(directory, f"{graph}.bed")
    fasta_tmp = os.path.join(directory, f".{graph}.fasta.tmp")
    bed_tmp = os.path.join(directory, f".{graph}.bed.tmp")
    prefix = graph_id_prefix(graph)
    bases = 0
    record_ids = set()
    try:
        with open(fasta_tmp, "wt") as fasta, open(bed_tmp, "wt") as bed:
            for template in templates:
                orientation = "_reverse" if template.strand == "-" else ""
                record_id = (
                    f"{prefix}_{sanitize_id(template.sample)}{orientation}_"
                    f"{sanitize_id(template.contig)}_"
                    f"{template.start}_{template.end}"
                )
                if record_id in record_ids:
                    raise ValueError(
                        f"{graph}: duplicate source-coordinate FASTA ID "
                        f"{record_id!r}"
                    )
                record_ids.add(record_id)
                fasta.write(
                    f">{record_id} {template.contig}:{template.start}-"
                    f"{template.end} block={graph} "
                    f"core={template.contig}:{template.core_start}-"
                    f"{template.core_end} "
                    f"role=reference sample={template.sample} "
                    f"strand={template.strand} origin={template.origin}"
                    + (" sequence_role=imported_alternative" if template.fixed else "")
                    + "\n"
                    +
                    f"{wrap_fasta(template.sequence)}\n"
                )
                bed.write(
                    f"{template.contig}\t{template.core_start}\t"
                    f"{template.core_end}\t{graph}\t"
                    f"1000\t{template.strand}\t{template.sample}\t"
                    f"{record_id}\t{template.start}\t{template.end}\t"
                    "source_coordinates\n"
                )
                bases += len(template.sequence)
        os.replace(fasta_tmp, fasta_path)
        os.replace(bed_tmp, bed_path)
    finally:
        for temporary in (fasta_tmp, bed_tmp):
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass
    return len(templates), bases


def remove_materialized_graph(graph: str, output_root: str) -> None:
    """Remove stale files previously generated for a now-filtered block."""
    directory = os.path.join(output_root, graph)
    for suffix in ("fasta", "bed"):
        try:
            os.remove(os.path.join(directory, f"{graph}.{suffix}"))
        except FileNotFoundError:
            pass
    try:
        os.rmdir(directory)
    except (FileNotFoundError, OSError):
        # Preserve a nonempty directory because it may contain unrelated files.
        pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build TARGET/TARGET.bed and TARGET/TARGET.fasta from one "
            "reference-to-block filtered BED."
        )
    )
    parser.add_argument(
        "-b", "--bed", required=True,
        help="single main_chroms_withnovel_blocks.fasta_aligned.filtered.bed",
    )
    parser.add_argument(
        "-q", "--query-paths", required=True,
        help="query_paths.txt; its first entry is the aligned reference",
    )
    parser.add_argument(
        "-t", "--target", required=True,
        help="main_chroms_withnovel_blocks.fasta used as the mapping target",
    )
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument(
        "--merged-groups",
        help=(
            "KmerPartitionMerger groups TSV; original filtered-BED targets "
            "are grouped into its final merged partition names"
        ),
    )
    parser.add_argument(
        "--merged-bed",
        help=(
            "KmerPartitionMerger BED4 output; automatically reads the "
            "companion PATH.groups.tsv"
        ),
    )
    parser.add_argument(
        "--unmask", type=int, default=1000, metavar="INT",
        help=(
            "minimum total uppercase A/C/G/T bases required across a block's "
            "output FASTA sequences (default: 1000; 0 disables filtering)"
        ),
    )
    parser.add_argument(
        "--anchor", type=int, default=0, metavar="INT",
        help=(
            "extend each FASTA source interval by this many bases per side "
            "after merging cores separated by less than twice this value; "
            "BED columns 2-3 remain unanchored (default: 0)"
        ),
    )
    parser.add_argument(
        "--target-origins-only", action="store_true",
        help=(
            "construct templates only from source= annotations in the target "
            "FASTA; intended for an exact mapping BED generated from a "
            "pre-grouped BED over a combined sequence catalog"
        ),
    )
    parser.add_argument("--keep-work", action="store_true")
    parser.add_argument(
        "--log-level", choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        default="INFO",
    )
    args = parser.parse_args(argv)
    if args.merged_bed and args.merged_groups:
        parser.error("use only one of --merged-bed and --merged-groups")
    if args.unmask < 0:
        parser.error("--unmask must be nonnegative")
    if args.anchor < 0:
        parser.error("--anchor must be nonnegative")
    return args


def run(args: argparse.Namespace) -> None:
    if shutil.which("sort") is None:
        raise FileNotFoundError("required executable not found on PATH: sort")
    input_bed = os.path.abspath(args.bed)
    target_fasta = os.path.abspath(args.target)
    required_paths = [input_bed, target_fasta, os.path.abspath(args.query_paths)]
    merged_groups = None
    merged_regions = None
    if args.merged_bed:
        merged_bed = os.path.abspath(args.merged_bed)
        merged_groups = merged_bed + ".groups.tsv"
        required_paths.extend((merged_bed, merged_groups))
        merged_regions = read_merged_bed(merged_bed)
    elif args.merged_groups:
        merged_groups = os.path.abspath(args.merged_groups)
        required_paths.append(merged_groups)
    for path in required_paths:
        if not os.path.isfile(path):
            raise FileNotFoundError(path)
    assemblies = read_assembly_list(os.path.abspath(args.query_paths))
    reference_sample = assemblies[0].sample
    reference_fasta = os.path.abspath(assemblies[0].path)
    assembly_paths = {
        assembly.sample: os.path.abspath(assembly.path)
        for assembly in assemblies
    }
    target_descriptions = read_fasta_descriptions(target_fasta)
    output_root = os.path.abspath(args.output_dir)
    Path(output_root).mkdir(parents=True, exist_ok=True)
    workdir = tempfile.mkdtemp(prefix=".build_localblocks.work.", dir=output_root)
    success = False
    try:
        sorted_bed = os.path.join(workdir, "filtered.by_target.bed")
        bed_to_sort = input_bed
        if merged_groups:
            groups = read_merged_groups(merged_groups)
            target_names = set(target_descriptions)
            missing_groups = target_names - set(groups)
            extra_groups = set(groups) - target_names
            if missing_groups or extra_groups:
                details = []
                if missing_groups:
                    details.append(
                        f"{len(missing_groups)} target FASTA records are missing"
                    )
                if extra_groups:
                    details.append(
                        f"{len(extra_groups)} group records are absent from target FASTA"
                    )
                raise ValueError(
                    f"{merged_groups} does not enumerate the target FASTA exactly: "
                    + "; ".join(details)
                )
            bed_to_sort = os.path.join(workdir, "filtered.remapped.bed")
            remap_filtered_bed(input_bed, bed_to_sort, groups)
        else:
            groups = {target: target for target in target_descriptions}
        graph_members: DefaultDict[str, List[str]] = defaultdict(list)
        for member, graph in groups.items():
            graph_members[graph].append(member)
        LOG.info(
            "Sorting one filtered BED by column 4; reference=%s (%s)",
            reference_sample, reference_fasta,
        )
        sort_filtered_bed(bed_to_sort, sorted_bed, workdir)
        graphs = 0
        graph_candidates = 0
        graphs_skipped_unmask = 0
        templates_total = 0
        bases_total = 0
        reference_bases_total = 0
        novel_bases_total = 0
        reference_core_bases_total = 0
        novel_core_bases_total = 0
        seen_graphs = set()
        with AssemblyReaderCache(assembly_paths) as source_readers, IndexedFasta(
            target_fasta,
        ) as target_reader:
            missing_targets = set(target_descriptions) - set(target_reader.index)
            if missing_targets:
                raise KeyError(
                    f"targets absent from {target_fasta}: "
                    f"{','.join(sorted(missing_targets))}"
                )

            def process_graph(graph: str, rows: Sequence[MappingRow]) -> None:
                nonlocal graphs, graph_candidates, graphs_skipped_unmask
                nonlocal templates_total, bases_total
                nonlocal reference_bases_total, novel_bases_total
                nonlocal reference_core_bases_total, novel_core_bases_total
                graph_candidates += 1
                members = graph_members.get(graph)
                if not members:
                    raise KeyError(
                        f"filtered BED group {graph!r} has no target FASTA members"
                    )
                if rows and merged_regions is not None and graph not in merged_regions:
                    raise KeyError(
                        f"merged group {graph!r} is absent from {args.merged_bed}"
                    )
                templates = build_templates(
                    graph, rows, reference_sample, source_readers,
                    target_reader, target_descriptions, args.anchor,
                    member_targets=members,
                    merged_query_intervals=(
                        merged_regions.get(graph)
                        if merged_regions is not None else None
                    ),
                    target_origins_only=args.target_origins_only,
                )
                unmasked = sum(
                    count_unmasked(template.sequence) for template in templates
                )
                if unmasked < args.unmask:
                    remove_materialized_graph(graph, output_root)
                    graphs_skipped_unmask += 1
                    LOG.debug(
                        "Skipped %s: %d unmasked bp is below --unmask %d",
                        graph, unmasked, args.unmask,
                    )
                    return
                count, bases = materialize_one_graph(
                    graph, templates, output_root,
                )
                graphs += 1
                templates_total += count
                bases_total += bases
                reference_bases_total += sum(
                    len(template.sequence)
                    for template in templates
                    if template.sample == reference_sample
                )
                novel_bases_total += sum(
                    len(template.sequence)
                    for template in templates
                    if template.sample != reference_sample
                )
                reference_core_bases_total += sum(
                    template.core_end - template.core_start
                    for template in templates
                    if template.sample == reference_sample
                )
                novel_core_bases_total += sum(
                    template.core_end - template.core_start
                    for template in templates
                    if template.sample != reference_sample
                )
                if graphs % 1000 == 0:
                    LOG.info("Built %d local block folders", graphs)

            for graph, rows in iter_target_groups(sorted_bed):
                if graph in seen_graphs:
                    raise ValueError(
                        f"sorted filtered BED contains non-contiguous duplicate "
                        f"group {graph!r}"
                    )
                seen_graphs.add(graph)
                process_graph(graph, rows)

            # A target can have no filtered lift at all. It still contributes
            # column-4 k-mers to KmerPartitionMerger and must remain available
            # as an initial local-graph path, either alone or in its merged group.
            for graph in graph_members:
                if graph in seen_graphs:
                    continue
                process_graph(graph, [])
        if graph_candidates == 0:
            raise ValueError(f"{target_fasta}: no target graph records")
        LOG.info(
            "Built %d/%d local blocks containing %d initial templates (%d bp) "
            "under %s; FASTA anchor=%d bp/side; skipped %d blocks below "
            "--unmask %d",
            graphs, graph_candidates, templates_total, bases_total, output_root,
            args.anchor, graphs_skipped_unmask, args.unmask,
        )
        LOG.info(
            "Total %s core size: %d bp (%.6f Gb); anchored FASTA size: "
            "%d bp (%.6f Gb)",
            reference_sample,
            reference_core_bases_total,
            reference_core_bases_total / 1e9,
            reference_bases_total,
            reference_bases_total / 1e9,
        )
        LOG.info(
            "Total novel-locus core size: %d bp (%.6f Gb); anchored FASTA "
            "size: %d bp (%.6f Gb)",
            novel_core_bases_total,
            novel_core_bases_total / 1e9,
            novel_bases_total,
            novel_bases_total / 1e9,
        )
        success = True
    finally:
        if success and not args.keep_work:
            shutil.rmtree(workdir, ignore_errors=True)
        elif not success:
            LOG.error("Retained failed-job work directory: %s", workdir)


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
    sys.exit(main())
