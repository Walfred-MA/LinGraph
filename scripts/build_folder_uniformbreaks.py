#!/usr/bin/env python3
"""Build folder-local gfixbreaks JSON caches from alignment/block TSVs.

All alignment samples participate by default. Use ``--onlyref`` for the
historical reference-only cache.
"""

from __future__ import annotations

import argparse
import collections as cl
import gc
import heapq
import json
import logging
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Dict, List, Optional, Sequence, Tuple

from align_partition_hotspots import OutputRow, read_graph_list
from summarize_partition_hotspot_segments import (
    ALIGNMENT_HEADER_FIELDS,
    query_name_matches_haplotype,
)
from uniform_graph_blocks import (
    INITIAL_BLOCK_HEADER,
    InitialBlock,
    load_graphfixbreaks,
    parse_alignment_components,
)


LOG = logging.getLogger("build_folder_uniformbreaks")
# Regenerate caches built with the old interval endpoints and chunk indexes.
# The JSON schema remains readable by existing graphDB consumers.
CACHE_FORMAT = "minsetref-gfixbreaks-folder-cache-v6"
_GFIXBREAKS: Optional[ModuleType] = None
_SHARED_ALIGNMENT_LINES: Dict[str, Tuple[str, ...]] = {}
_SHARED_BLOCK_LINES: Dict[str, Tuple[Tuple[int, str], ...]] = {}


@dataclass(frozen=True)
class FolderTask:
    hotspot_index: int
    graph_name: str
    graph_fasta: str
    row_count: int
    block_count: int
    output: str
    source_alignment: str
    source_blocks: str
    graphfixbreaks: str
    resume: bool
    alignment_scope: str = "all_samples"
    blocked_output: str = ""
    blocked_output_count: int = 0
    assemblysmall_scope: str = ""
    assemblysmall_cutoff: int = 0
    assemblysmall_links: Tuple[Tuple[int, int], ...] = ()
    assemblysmall_paths: Tuple[Tuple[int, ...], ...] = ()


def _init_worker(graphfixbreaks: str) -> None:
    global _GFIXBREAKS
    # The two large dictionaries are inherited read-only through fork. Turning
    # off cyclic GC prevents worker collections from dirtying their pages.
    gc.disable()
    if not _SHARED_ALIGNMENT_LINES or not _SHARED_BLOCK_LINES:
        raise RuntimeError("shared reference-line dictionaries are unavailable")
    _GFIXBREAKS = load_graphfixbreaks(graphfixbreaks)


def split_fields(raw: str) -> List[str]:
    fields = raw.rstrip("\r\n").split("\t")
    return fields if len(fields) > 1 else raw.split()


def parse_alignment_line(raw: str, context: str) -> OutputRow:
    fields = split_fields(raw)
    if len(fields) != len(ALIGNMENT_HEADER_FIELDS):
        raise ValueError(
            f"{context}: expected ten alignment columns, found {len(fields)}"
        )
    try:
        hotspot_index = int(fields[0])
        start = int(fields[2])
        end = int(fields[3])
    except ValueError as error:
        raise ValueError(f"{context}: invalid integer alignment field") from error
    if (
        hotspot_index < 1 or start < 0 or end <= start
        or fields[4] not in {"+", "-"}
    ):
        raise ValueError(f"{context}: invalid alignment coordinates/strand")
    return OutputRow(
        hotspot_index, fields[1], start, end, fields[4], fields[5],
        fields[6], fields[7], fields[8], fields[9],
    )


def parse_block_line(raw: str, block_id: int, context: str) -> InitialBlock:
    fields = split_fields(raw)
    if len(fields) != len(INITIAL_BLOCK_HEADER):
        raise ValueError(
            f"{context}: expected nine block columns, found {len(fields)}"
        )
    try:
        hotspot, qstart, qend, gstart, gend, count = (
            int(fields[index]) for index in (0, 3, 4, 6, 7, 8)
        )
    except ValueError as error:
        raise ValueError(f"{context}: invalid integer block field") from error
    if (
        hotspot < 1 or qstart < 0 or qend <= qstart
        or gstart < 0 or gend <= gstart or count < 1
        or fields[5] not in {"+", "-"}
    ):
        raise ValueError(f"{context}: invalid 0-based half-open block")
    return InitialBlock(
        block_id, hotspot, fields[1], fields[2], qstart, qend,
        fields[5], gstart, gend, count,
    )


def read_graph_path_lengths(path: str) -> Dict[str, int]:
    """Read every graph FASTA record and return its actual sequence length."""
    output: Dict[str, int] = {}
    name: Optional[str] = None
    length = 0
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith(">"):
                if name is not None:
                    output[name] = length
                fields = raw[1:].split()
                if not fields:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                name = fields[0]
                if name in output:
                    raise ValueError(f"{path}:{line_number}: duplicate path {name!r}")
                length = 0
            else:
                if name is None and raw.strip():
                    raise ValueError(f"{path}:{line_number}: sequence before header")
                length += len(raw.strip())
    if name is not None:
        output[name] = length
    if not output:
        raise ValueError(f"graph FASTA contains no records: {path}")
    return output


def oriented_overlap_on_row(
    row: OutputRow, block: InitialBlock,
) -> Optional[Tuple[int, int]]:
    if row.contig != block.query_contig:
        return None
    start = max(row.start, block.query_start)
    end = min(row.end, block.query_end)
    if end <= start:
        return None
    if row.strand == "+":
        return start - row.start, end - row.start
    return row.end - end, row.end - start


def graphinfo_from_rows(
    module: ModuleType,
    rows: Sequence[OutputRow],
    path_lengths: Dict[str, int],
) -> Tuple[
    Dict[str, str], Dict[str, list], Dict[int, str], Dict[str, list]
]:
    graphinfo: Dict[str, str] = {}
    contigspan: Dict[str, list] = {}
    row_loci: Dict[int, str] = {}
    locuslocations: Dict[str, list] = {}
    for row_index, row in enumerate(rows, 1):
        components = parse_alignment_components(row)
        if not components:
            continue
        for component in components:
            observed = path_lengths.get(component.path)
            if observed is None:
                raise ValueError(
                    f"hotspot {row.hotspot_index} row {row.query_name}: graph "
                    f"path {component.path!r} is absent from its .FA"
                )
            if observed != component.path_length:
                raise ValueError(
                    f"hotspot {row.hotspot_index} path {component.path!r}: "
                    f"name implies {component.path_length} bp but .FA has "
                    f"{observed} bp"
                )
        locus = f"reference_alignment_{row_index:06d}"
        cigar, spans = module.getspan_ongraph(
            locus, row.graph_path, row.graph_cigar,
            row.ref_positions, row.query_positions,
        )
        if not cigar or not spans:
            continue
        graphinfo[locus] = cigar
        contigspan[locus] = spans
        row_loci[row_index - 1] = locus
        reference_locus = "Ref_" + locus
        # This is the same structure produced by readcontigsinfo() for a
        # reference record in *_loci.txt.fasta.  The suffix keeps separate
        # alignment regions on the same chromosome as separate loci.
        scaffold_locus = f"{row.contig}_{row_index}"
        locuslocations[reference_locus] = [
            scaffold_locus, row.strand, row.start, row.end,
        ]
    return graphinfo, contigspan, row_loci, locuslocations


def best_row_overlap(
    rows: Sequence[OutputRow],
    row_loci: Dict[int, str],
    block: InitialBlock,
) -> Optional[Tuple[int, int, int]]:
    """Assign an input block to its best matching synthetic loci record."""
    candidates = []
    for row_index, row in enumerate(rows):
        if row_index not in row_loci:
            continue
        overlap = oriented_overlap_on_row(row, block)
        if overlap is None:
            continue
        absolute_start = max(row.start, block.query_start)
        absolute_end = min(row.end, block.query_end)
        contained = (
            row.start <= block.query_start and row.end >= block.query_end
        )
        candidates.append((
            int(contained), absolute_end - absolute_start,
            -abs(row.start - block.query_start)
            - abs(row.end - block.query_end),
            -row_index, row_index, overlap[0], overlap[1],
        ))
    if not candidates:
        return None
    selected = max(candidates)
    return selected[4], selected[5], selected[6]


def set_validregions_from_annotations(
    module: ModuleType,
    graphdb,
    regions: Sequence[Sequence[object]],
    original_loci: Dict[str, Tuple[str, int, int]],
    graphinfo: Dict[str, str],
    contigspan: Dict[str, list],
) -> None:
    """Record every graph segment traversed by each reference annotation.

    The legacy implementation projected only ``local_start`` and
    ``local_end - 1``.  If a block traversed ``A -> B -> A``, both endpoints
    mapped to A and the internal B path was absent from ``validregions``.
    Project the complete half-open query interval instead, which also removes
    the legacy inclusive-end ambiguity.
    """
    if not hasattr(module, "build_total_validregions"):
        raise RuntimeError(
            "gfixbreaks lacks build_total_validregions; use the bundled "
            "gfixbreaks.py"
        )
    graphdb.validregions = module.build_total_validregions(
        regions, original_loci, graphinfo, contigspan,
    )


def uniformbreaks_from_blocks(
    module: ModuleType,
    rows: Sequence[OutputRow],
    blocks: Sequence[InitialBlock],
    graph_fasta: str,
):
    """Adapt reference blocks/alignments to legacy ``uniformbreaks`` inputs."""
    path_lengths = read_graph_path_lengths(graph_fasta)
    graphinfo, contigspan, row_loci, locuslocations = graphinfo_from_rows(
        module, rows, path_lengths,
    )
    if not graphinfo:
        raise ValueError("no usable reference graph-alignment rows")

    genetoscaff = {}
    genetolocus = cl.defaultdict(list)
    for block in blocks:
        block_name = f"reference_block_{block.block_id:06d}"
        match = best_row_overlap(rows, row_loci, block)
        if match is None:
            continue
        row_index, local_start, local_end = match
        reference_locus = "Ref_" + row_loci[row_index]
        scaffold_locus = locuslocations[reference_locus][0]
        # Six fields exactly mirror an inputfile FASTA header loaded by
        # readcontigsinfo(): contig, strand, start, end, annotation, mapinfo.
        genetoscaff[block_name] = [
            scaffold_locus, block.strand,
            block.query_start, block.query_end, "", "",
        ]
        genetolocus[block_name].append([
            reference_locus, local_start, local_end,
        ])

    usable_block_count = len(genetoscaff)
    missing = [
        f"reference_block_{block.block_id:06d}"
        for block in blocks
        if f"reference_block_{block.block_id:06d}" not in genetolocus
    ]
    if missing:
        graph_name = os.path.basename(graph_fasta)[:-3]
        LOG.warning(
            "%s: skipping %d reference block%s that do not overlap a "
            "usable reference graph-alignment row; first: %s",
            graph_name, len(missing), "" if len(missing) == 1 else "s",
            missing[0],
        )

    # readcontigsinfo() also registers every loci FASTA record in
    # genetoscaff.  Those four-field records are ignored as input blocks by
    # annotate_regions(), but retain the original uniformbreaks structure.
    for reference_locus, location in locuslocations.items():
        genetoscaff[reference_locus] = list(location)

    # This is intentionally the same two-function core as uniformbreaks().
    # Reference alignments do not use a blacklist.
    module.mergesmall = 20_000
    results, graphdb = module.breaks_ongraph(
        genetoscaff, genetolocus, graphinfo, contigspan, set(),
    )
    regions, original_loci = module.annotate_regions(
        results, locuslocations, genetoscaff,
    )
    set_validregions_from_annotations(
        module, graphdb, regions, original_loci, graphinfo, contigspan,
    )
    return (
        graphdb, path_lengths, usable_block_count, len(missing), len(regions),
    )


def cache_is_reusable(task: FolderTask) -> bool:
    if not task.resume or not os.path.isfile(task.output):
        return False
    try:
        with open(task.output, "rt") as handle:
            data = json.load(handle)
        metadata = data.get("_minsetref", {})
        if (
            metadata.get("format") != CACHE_FORMAT
            or int(metadata.get("hotspot_index", -1)) != task.hotspot_index
            or metadata.get("alignment_scope", "all_samples")
            != task.alignment_scope
        ):
            return False
        output_time = os.path.getmtime(task.output)
        # CACHE_FORMAT represents output semantics. Implementation-only
        # optimizations in gfixbreaks.py must not invalidate tens of thousands
        # of equivalent completed caches; semantic changes bump CACHE_FORMAT.
        return output_time >= max(os.path.getmtime(path) for path in (
            task.graph_fasta, task.source_alignment,
            task.source_blocks,
        ))
    except (OSError, ValueError, TypeError, json.JSONDecodeError):
        return False


def write_graphdb_atomic(
    module: ModuleType,
    graphdb,
    task: FolderTask,
    path_lengths: Dict[str, int],
    usable_block_count: int,
    skipped_block_count: int,
    uniform_region_count: int,
) -> None:
    absolute = os.path.abspath(task.output)
    Path(os.path.dirname(absolute)).mkdir(parents=True, exist_ok=True)
    temporary = absolute + f".tmp.{os.getpid()}"
    try:
        # Serialize directly from graphDB. Calling graphdb.save_json() and
        # reading that file back with json.load() temporarily retained both
        # the complete graphDB and a second complete decoded JSON object.
        data = {
            "breaks": graphdb.breaks,
            "blockset": {
                f"{key[0]}|{key[1]}|{key[2]}": value
                for key, value in graphdb.blockset.items()
            },
            "blocksize": graphdb.blocksize,
            "chunkindex_togene": {
                str(key): [[list(item[0]), item[1]] for item in value]
                for key, value in graphdb.chunkindex_togene.items()
            },
            "genechunks_index": {
                ",".join(map(str, key)): value
                for key, value in graphdb.genechunks_index.items()
            },
            "validregions": graphdb.validregions,
            "_minsetref": {
                "format": CACHE_FORMAT,
                "coordinate_system": "0-based-half-open",
                "validregions_semantics": (
                    "reference-block-graph-closure-v1"
                    if task.alignment_scope == "only_reference" else
                    "cohort-source-block-graph-closure-v1"
                ),
                "hotspot_index": task.hotspot_index,
                "graph_name": os.path.basename(task.graph_fasta)[:-3],
                "graph_fasta": os.path.abspath(task.graph_fasta),
                "reference_alignment": os.path.abspath(task.source_alignment),
                "reference_blocks": os.path.abspath(task.source_blocks),
                "graphfixbreaks": os.path.abspath(task.graphfixbreaks),
                "alignment_scope": task.alignment_scope,
                "alignment_row_count": task.row_count,
                "block_count": task.block_count,
                "uniform_region_count": uniform_region_count,
                "reference_row_count": task.row_count,
                "reference_block_count": task.block_count,
                "usable_reference_block_count": usable_block_count,
                "skipped_reference_block_count": skipped_block_count,
                "uniform_reference_region_count": uniform_region_count,
                "graph_path_lengths": path_lengths,
            },
        }
        if task.blocked_output:
            data["_minsetref"]["blocked_output"] = os.path.abspath(
                task.blocked_output
            )
            data["_minsetref"]["blocked_output_count"] = int(
                task.blocked_output_count
            )
        if task.assemblysmall_scope:
            # Imported lazily to keep the standalone folder-cache builder
            # independent of the blocking CLI at module-import time.
            from block_partition_alignments import (
                COHORT_ASSEMBLYSMALL_KEY,
                COHORT_ASSEMBLYSMALL_VERSION,
            )
            data["_minsetref"][COHORT_ASSEMBLYSMALL_KEY] = {
                "version": COHORT_ASSEMBLYSMALL_VERSION,
                "scope": task.assemblysmall_scope,
                "cutoff": int(task.assemblysmall_cutoff),
                "within_block_links": [
                    list(link) for link in task.assemblysmall_links
                ],
                "within_block_paths": [
                    list(path) for path in task.assemblysmall_paths
                ],
            }
        mergesmall = getattr(module, "mergesmall", 20_000)
        if mergesmall != 20_000:
            data["mergesmall"] = str(mergesmall)
        with open(temporary, "wt") as output:
            json.dump(data, output, sort_keys=True, separators=(",", ":"))
            output.write("\n")
        os.replace(temporary, absolute)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def _process_task(task: FolderTask) -> Tuple[int, str, str]:
    if cache_is_reusable(task):
        return task.hotspot_index, "reused", task.output
    if _GFIXBREAKS is None:
        raise RuntimeError("gfixbreaks worker was not initialized")
    alignment_lines = _SHARED_ALIGNMENT_LINES.get(task.graph_name, ())
    block_lines = _SHARED_BLOCK_LINES.get(task.graph_name, ())
    if len(alignment_lines) != task.row_count:
        raise RuntimeError(
            f"{task.graph_name}: shared alignment-row count changed"
        )
    if len(block_lines) != task.block_count:
        raise RuntimeError(f"{task.graph_name}: shared block count changed")
    # Only this graph is parsed into Python field objects. The large raw-line
    # dictionaries remain shared and read-only in the forked parent heap.
    rows = [
        parse_alignment_line(
            raw, f"{task.source_alignment}:{task.graph_name}:row{index}",
        )
        for index, raw in enumerate(alignment_lines, 1)
    ]
    blocks = [
        parse_block_line(
            raw, block_id,
            f"{task.source_blocks}:{task.graph_name}:block{block_id}",
        )
        for block_id, raw in block_lines
    ]
    (
        graphdb, path_lengths, usable_block_count, skipped_block_count,
        uniform_region_count,
    ) = uniformbreaks_from_blocks(
        _GFIXBREAKS, rows, blocks, task.graph_fasta,
    )
    write_graphdb_atomic(
        _GFIXBREAKS, graphdb, task, path_lengths,
        usable_block_count, skipped_block_count, uniform_region_count,
    )
    return task.hotspot_index, "built", task.output


def _process_batch(tasks: Sequence[FolderTask]) -> List[Tuple[int, str, str]]:
    results = []
    for task in tasks:
        LOG.info(
            "Worker %d START %s (hotspot %d; %d rows; %d blocks)",
            os.getpid(), task.graph_name, task.hotspot_index,
            task.row_count, task.block_count,
        )
        try:
            result = _process_task(task)
        except Exception as error:
            raise RuntimeError(f"{task.graph_name}: {error}") from error
        results.append(result)
        LOG.info(
            "Worker %d DONE %s (%s)",
            os.getpid(), task.graph_name, result[1],
        )
    return results


def balanced_batches(
    tasks: Sequence[FolderTask], count: int,
) -> List[List[FolderTask]]:
    batches: List[List[FolderTask]] = [[] for _index in range(count)]
    heap = [(0, index) for index in range(count)]
    heapq.heapify(heap)
    for task in sorted(
        tasks, key=lambda value: -(value.row_count + value.block_count),
    ):
        load, index = heapq.heappop(heap)
        batches[index].append(task)
        heapq.heappush(
            heap, (load + task.row_count + task.block_count, index),
        )
    return [batch for batch in batches if batch]


def first_field(raw: str) -> str:
    stripped = raw.lstrip()
    if not stripped:
        return ""
    tab = stripped.find("\t")
    space = stripped.find(" ")
    ends = [value for value in (tab, space) if value >= 0]
    return stripped[:min(ends)] if ends else stripped.rstrip("\r\n")


def load_shared_line_maps(
    args: argparse.Namespace,
    graph_fastas: Sequence[Optional[str]],
) -> Tuple[Dict[str, int], Dict[str, int], int]:
    """Load raw source lines once into read-only graph-name dictionaries."""
    global _SHARED_ALIGNMENT_LINES, _SHARED_BLOCK_LINES
    graph_names = [
        os.path.basename(path)[:-3] if path is not None else None
        for path in graph_fastas
    ]
    duplicate_names = [
        name for name, count in cl.Counter(
            value for value in graph_names if value is not None
        ).items() if count > 1
    ]
    if duplicate_names:
        raise ValueError(
            "graph list has duplicate graph names; shared graph-name keys "
            f"would be ambiguous: {duplicate_names[0]}"
        )

    alignment_lists = cl.defaultdict(list)
    alignment_bytes = 0
    with open(args.alignment, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            field = first_field(raw)
            if not field or raw.startswith("#") or field == ALIGNMENT_HEADER_FIELDS[0]:
                continue
            try:
                hotspot_index = int(field)
            except ValueError as error:
                raise ValueError(
                    f"{args.alignment}:{line_number}: invalid hotspot index"
                ) from error
            if hotspot_index < 1 or hotspot_index > len(graph_names):
                raise ValueError(
                    f"{args.alignment}:{line_number}: hotspot index "
                    f"{hotspot_index} exceeds graph-list length "
                    f"{len(graph_names)}"
                )
            graph_name = graph_names[hotspot_index - 1]
            if graph_name is None:
                continue
            if args.onlyref:
                fields = split_fields(raw)
                if len(fields) != len(ALIGNMENT_HEADER_FIELDS):
                    raise ValueError(
                        f"{args.alignment}:{line_number}: expected ten "
                        "alignment columns"
                    )
                if not query_name_matches_haplotype(
                    fields[5], args.reference_haplotype,
                ):
                    continue
            alignment_lists[graph_name].append(raw)
            alignment_bytes += len(raw)

    block_lists = cl.defaultdict(list)
    block_bytes = 0
    block_id = 0
    with open(args.blocks, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            field = first_field(raw)
            if not field or raw.startswith("#") or field == INITIAL_BLOCK_HEADER[0]:
                continue
            block_id += 1
            fields = split_fields(raw)
            if len(fields) != len(INITIAL_BLOCK_HEADER):
                raise ValueError(
                    f"{args.blocks}:{line_number}: expected nine block columns, "
                    f"found {len(fields)}"
                )
            try:
                hotspot_index = int(fields[0])
            except ValueError as error:
                raise ValueError(
                    f"{args.blocks}:{line_number}: invalid hotspot index"
                ) from error
            if hotspot_index < 1 or hotspot_index > len(graph_names):
                raise ValueError(
                    f"{args.blocks}:{line_number}: hotspot index "
                    f"{hotspot_index} exceeds graph-list length "
                    f"{len(graph_names)}"
                )
            graph_name = graph_names[hotspot_index - 1]
            if graph_name is None:
                continue
            if fields[1] != graph_name:
                raise ValueError(
                    f"{args.blocks}:{line_number}: graph list resolves to "
                    f"{graph_name!r}, but block names {fields[1]!r}"
                )
            block_lists[graph_name].append((block_id, raw))
            block_bytes += len(raw)

    _SHARED_ALIGNMENT_LINES = {
        name: tuple(lines) for name, lines in alignment_lists.items()
    }
    _SHARED_BLOCK_LINES = {
        name: tuple(lines) for name, lines in block_lists.items()
    }
    return (
        {name: len(lines) for name, lines in _SHARED_ALIGNMENT_LINES.items()},
        {name: len(lines) for name, lines in _SHARED_BLOCK_LINES.items()},
        alignment_bytes + block_bytes,
    )


def build_tasks(args: argparse.Namespace) -> Tuple[List[FolderTask], int, int]:
    graph_fastas = read_graph_list(args.graph_list, args.graph_folder)
    row_counts, block_counts, shared_bytes = load_shared_line_maps(
        args, graph_fastas,
    )

    tasks = []
    skipped = 0
    for hotspot_index, graph_fasta in enumerate(graph_fastas, 1):
        if graph_fasta is None:
            continue
        graph_name = os.path.basename(graph_fasta)[:-3]
        if graph_name not in block_counts:
            continue
        if graph_name not in row_counts:
            raise ValueError(
                f"hotspot {hotspot_index} has reference blocks but no "
                "reference graph-alignment rows"
            )
        tasks.append(FolderTask(
            hotspot_index,
            graph_name,
            graph_fasta,
            row_counts[graph_name],
            block_counts[graph_name],
            graph_fasta[:-3] + "cache.json",
            os.path.abspath(args.alignment),
            os.path.abspath(args.blocks),
            os.path.abspath(args.graphfixbreaks),
            args.resume,
            "only_reference" if args.onlyref else "all_samples",
        ))
    skipped = sum(path is None for path in graph_fastas)
    return tasks, skipped, shared_bytes


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--alignment", required=True,
        help="ten-column align_partition_hotspots.py output",
    )
    parser.add_argument(
        "-b", "--blocks", required=True,
        help="summarize_partition_hotspot_segments.py nine-column output",
    )
    parser.add_argument(
        "-G", "--graph-folder", required=True,
        help="root containing GRAPH_NAME/GRAPH_NAME.FA",
    )
    parser.add_argument(
        "-L", "--graph-list", required=True,
        help="ordered graph list used by the alignment",
    )
    parser.add_argument(
        "--graphfixbreaks",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "gfixbreaks.py"),
        help="gfixbreaks.py implementation (default: beside this script)",
    )
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--onlyref", "--only-ref", action="store_true",
        help="restrict cache construction to named reference alignment rows",
    )
    parser.add_argument(
        "--reference-haplotype", action="append",
        help="reference haplotype used by --onlyref; repeat for aliases",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if not args.reference_haplotype:
        args.reference_haplotype = ["CHM13_h1"]
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )
    try:
        tasks, skipped, shared_bytes = build_tasks(args)
        if not tasks:
            raise ValueError("no graph folders have blocks to cache")
        workers = min(args.jobs, len(tasks))
        LOG.info(
            "Loaded %.2f GiB of alignment/block text into shared "
            "read-only dictionaries keyed by graph name",
            shared_bytes / 1024 ** 3,
        )
        LOG.info(
            "Building folder-local uniformbreaks caches for %d graphs with "
            "%d worker process%s; skipped %d missing graph%s",
            len(tasks), workers, "" if workers == 1 else "es",
            skipped, "" if skipped == 1 else "s",
        )
        built = reused = completed = 0
        if workers == 1:
            _init_worker(args.graphfixbreaks)
            result_batches = [_process_batch(tasks)]
        else:
            if "fork" not in multiprocessing.get_all_start_methods():
                raise RuntimeError(
                    "shared reference dictionaries require the fork "
                    "multiprocessing start method"
                )
            batches = balanced_batches(
                tasks, min(len(tasks), workers * 4),
            )
            result_batches = []
            # Freeze all currently tracked objects before fork. Workers inherit
            # the large raw-line dictionaries through copy-on-write pages and
            # receive only small graph-name tasks through the process queue.
            gc.collect()
            gc.freeze()
            try:
                with ProcessPoolExecutor(
                    max_workers=workers,
                    mp_context=multiprocessing.get_context("fork"),
                    initializer=_init_worker,
                    initargs=(args.graphfixbreaks,),
                ) as executor:
                    futures = [
                        executor.submit(_process_batch, batch)
                        for batch in batches
                    ]
                    report_every = max(1, len(futures) // 20)
                    for batch_number, future in enumerate(
                        as_completed(futures), 1,
                    ):
                        results = future.result()
                        result_batches.append(results)
                        completed += len(results)
                        if (
                            batch_number % report_every == 0
                            or batch_number == len(futures)
                        ):
                            LOG.info(
                                "Folder-cache progress: %d/%d graphs",
                                completed, len(tasks),
                            )
            finally:
                gc.unfreeze()
        for results in result_batches:
            for _index, status, _path in results:
                built += status == "built"
                reused += status == "reused"
        LOG.info(
            "Folder uniformbreaks complete: built %d; reused %d",
            built, reused,
        )
        return 0
    except Exception as error:
        LOG.error("%s", error)
        if args.verbose:
            LOG.exception("folder uniformbreaks generation failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
