#!/usr/bin/env python3
"""Block existing whole-CIGAR assembly alignments with graph cache.json files.

The input is the ten-column output from ``align_partition_hotspots.py``.  Each
hotspot index selects a graph through the same ordered graph list used during
alignment.  The corresponding folder-local ``GRAPH_NAMEcache.json`` is loaded
and the already-computed graph CIGAR is passed through the LocalGraphAlign
blocking path without running an aligner again.

Output is headerless BED6, matching ``AllGraphAlign.py`` when its output name
ends in ``.bed``::

    chrom  start  end  name  score  strand

Coordinates are 0-based, half-open assembly coordinates.  Use ``--header`` if
a header is desired.  Run the reference once with ``--update-cache-merges``;
this stores the final reference block paths and their graph-chunk adjacencies
in each cache. Later sample runs replay those reference decisions instead of
making sample-specific size/distance merge choices.
"""

from __future__ import annotations

import argparse
import collections as cl
import copy
import dataclasses
import gc
import json
import logging
import multiprocessing
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from types import ModuleType
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from align_partition_hotspots import OutputRow, graph_name_from_list_entry
from summarize_partition_hotspot_segments import (
    read_alignment_output,
    select_alignment_rows,
)
from uniform_graph_blocks import (
    absolute_to_oriented_query,
    graph_interval_to_query,
    load_graphfixbreaks,
    parse_alignment_components,
)


LOG = logging.getLogger("block_partition_alignments")
DEFAULT_MERGE_SMALL = 20_000


@dataclasses.dataclass(frozen=True)
class GraphTask:
    hotspot_index: int
    graph_name: str
    cache_path: str
    source_bed_path: str
    row_count: int


@dataclasses.dataclass(frozen=True)
class BlockedRow:
    input_index: int
    region_index: int
    hotspot_index: int
    prefix: str
    contig: str
    strand: str
    start: int
    end: int

    def line(self) -> str:
        return "\t".join(str(value) for value in (
            self.contig, self.start, self.end, self.prefix,
            self.hotspot_index, self.strand,
        )) + "\n"


_GFIXBREAKS: Optional[ModuleType] = None
_ROWS_BY_GRAPH: Dict[str, Tuple[Tuple[int, OutputRow], ...]] = {}
_MERGE_SMALL = DEFAULT_MERGE_SMALL
_UPDATE_CACHE_MERGES = False
_REFERENCE_HAPLOTYPES: Set[str] = {"CHM13_h1"}
REFERENCE_ASSEMBLYSMALL_KEY = "reference_assemblysmall"
REFERENCE_ASSEMBLYSMALL_VERSION = 4
COHORT_ASSEMBLYSMALL_KEY = "cohort_assemblysmall"
# Version 2 invalidates cohort controls whose initial in-memory blocking pass
# could interpret half-open valid-region endpoints as legacy inclusive.  The
# serialized cache declares them half-open, so replay could otherwise differ
# by one base.  The graphDB base remains reusable; only cohort learning/replay
# and its BED are regenerated.
COHORT_ASSEMBLYSMALL_VERSION = 2


def read_graph_names(path: str) -> List[str]:
    names: List[str] = []
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            content = raw.strip()
            if not content or content.startswith("#"):
                continue
            names.append(graph_name_from_list_entry(
                content.split()[0], f"{path}:{line_number}",
            ))
    if not names:
        raise ValueError(f"graph list is empty: {path}")
    duplicate = next(
        (name for name, count in cl.Counter(names).items() if count > 1),
        None,
    )
    if duplicate is not None:
        raise ValueError(f"graph list contains duplicate graph name {duplicate!r}")
    return names


def cache_path_for_graph(graph_root: str, graph_name: str) -> str:
    folder = os.path.join(graph_root, graph_name)
    candidates = (
        os.path.join(folder, graph_name + "cache.json"),
        os.path.join(folder, "cache.json"),
    )
    for candidate in candidates:
        if os.path.isfile(candidate):
            return candidate
    return candidates[0]


def graph_prefix(graph_name: str) -> str:
    """Match AllGraphAlign.py's basename-before-first-underscore prefix."""
    return graph_name.split("_", 1)[0]


def normalize_intervals(values) -> List[Tuple[int, int]]:
    if not values:
        return []
    if isinstance(values[0], (int, float, str)):
        if len(values) % 2:
            raise ValueError("cached validregions contains an odd coordinate count")
        intervals = [
            (int(values[index]), int(values[index + 1]))
            for index in range(0, len(values), 2)
        ]
    else:
        intervals = [(int(value[0]), int(value[1])) for value in values]
    normalized = sorted(
        (min(start, end), max(start, end))
        for start, end in intervals if start != end
    )
    merged: List[List[int]] = []
    for start, end in normalized:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _valid_query_projection_details(
    row: OutputRow,
    validregions: Mapping[str, object],
    coordinate_system: str = "",
    occurrence_aware: bool = False,
) -> List[Dict[str, object]]:
    """Project every valid graph interval, independent of known block genes.

    With ``occurrence_aware=True``, pieces are unioned within one graph
    traversal occurrence.  The compatibility default retains the historical
    all-occurrences hull.  A path name is not an occurrence identifier:
    graph-CIGARs may visit the same path more than once.
    """
    components = parse_alignment_components(row)
    output: List[Dict[str, object]] = []
    if occurrence_aware:
        groups = [
            (index, [component])
            for index, component in enumerate(components)
            if component.path in validregions
        ]
    else:
        # Compatibility mode for callers of this low-level helper.  The
        # production blocking path uses occurrence-aware projection below.
        grouped = cl.defaultdict(list)
        for component in components:
            grouped[component.path].append(component)
        groups = [
            (-1, path_components)
            for path, path_components in grouped.items()
            if path in validregions
        ]

    for component_index, path_components in groups:
        path = path_components[0].path
        if path not in validregions:
            continue
        for graph_start, cached_end in normalize_intervals(
            validregions[path],
        ):
            graph_end = (
                cached_end
                if coordinate_system == "0-based-half-open"
                else cached_end + 1
            )
            projections = [
                projection
                for component in path_components
                for projection in graph_interval_to_query(
                    row, component, graph_start, graph_end,
                )
            ]
            local_ranges = sorted(set(
                absolute_to_oriented_query(
                    row, projection.query_start, projection.query_end,
                )
                for projection in projections
            ))
            if not local_ranges:
                continue
            output.append({
                "path": path,
                "graph_start": graph_start,
                "graph_end": graph_end,
                "query_start": min(value[0] for value in local_ranges),
                "query_end": max(value[1] for value in local_ranges),
                "query_pieces": local_ranges,
                "orientation": path_components[0].marker,
                "component_index": component_index,
            })
    return sorted(output, key=lambda value: (
        int(value["query_start"]), int(value["query_end"]),
        str(value["path"]), int(value.get("component_index", 0)),
        int(value["graph_start"]),
    ))


def valid_query_intervals(
    row: OutputRow,
    validregions: Mapping[str, object],
    coordinate_system: str = "",
    occurrence_aware: bool = False,
) -> List[Tuple[int, int]]:
    """Project all cached graph spans to oriented query spans."""
    return [
        (int(value["query_start"]), int(value["query_end"]))
        for value in _valid_query_projection_details(
            row, validregions, coordinate_system, occurrence_aware,
        )
    ]


def lossless_projection_results(
    row: OutputRow, graphdb, module: ModuleType,
    query_valid_intervals: Optional[Sequence[Tuple[int, int]]] = None,
) -> Dict[str, List[list]]:
    """Build lossless blocks while retaining known reference representation."""
    if not hasattr(module, "lossless_validregion_results"):
        raise RuntimeError(
            "gfixbreaks lacks lossless_validregion_results; use the bundled "
            "gfixbreaks.py"
        )
    graph_cigar, spans = module.getspan_ongraph(
        row.query_name, row.graph_path, row.graph_cigar,
        row.ref_positions, row.query_positions,
    )
    try:
        return module.lossless_validregion_results(
            graphdb,
            {row.query_name: graph_cigar},
            {row.query_name: spans},
            occurrence_aware=True,
            query_valid_intervals=(
                None if query_valid_intervals is None
                else {row.query_name: tuple(query_valid_intervals)}
            ),
        )
    except TypeError as error:
        # Keep compatibility with older external gfixbreaks modules that do
        # not yet expose occurrence-aware projection.
        if (
            "occurrence_aware" not in str(error)
            and "query_valid_intervals" not in str(error)
        ):
            raise
        if query_valid_intervals is not None:
            raise RuntimeError(
                "gfixbreaks lacks reference source-window masking; use the "
                "bundled gfixbreaks.py"
            ) from error
        return module.lossless_validregion_results(
            graphdb,
            {row.query_name: graph_cigar},
            {row.query_name: spans},
        )


def _region_graph_chunks(region: Sequence[object]) -> List[int]:
    if len(region) < 2:
        return []
    value = region[1]
    if isinstance(value, (list, tuple)):
        return [int(item) for item in value]
    return [int(value)]


def _material_region_graph_chunks(region: Sequence[object]) -> List[int]:
    """Return chunks with more than a one-base query contribution.

    Legacy folder caches can place a graph-chunk boundary one base before the
    corresponding half-open valid-region boundary.  Lossless projection then
    reports the next chunk with a one-base query segment.  Treating that
    boundary sliver as a complete reference-supported chunk causes cached
    path replay to promote the whole terminal anchor in every sample.

    Keep the original chunk list when it cannot be paired unambiguously with
    projected query segments.  The filter is therefore limited to the exact
    legacy artifact and does not discard ordinary short graph chunks.
    """
    chunks = _region_graph_chunks(region)
    if len(region) <= 5 or not isinstance(region[5], (list, tuple)):
        return chunks
    segments = region[5]
    if len(segments) != len(chunks):
        return chunks
    material = []
    for chunk, segment in zip(chunks, segments):
        if not isinstance(segment, (list, tuple)) or len(segment) < 2:
            return chunks
        if abs(int(segment[1]) - int(segment[0])) > 1:
            material.append(chunk)
    return material


def _boundary_link(left: Sequence[object], right: Sequence[object]) -> Optional[Tuple[int, int]]:
    left_chunks = _region_graph_chunks(left)
    right_chunks = _region_graph_chunks(right)
    if not left_chunks or not right_chunks:
        return None
    return tuple(sorted((left_chunks[-1], right_chunks[0])))


def result_within_block_links(
    results: Mapping[str, Sequence[list]],
) -> Set[Tuple[int, int]]:
    """Return graph-chunk adjacencies contained inside final result blocks."""
    links: Set[Tuple[int, int]] = set()
    for regions in results.values():
        for region in regions:
            chunks = _material_region_graph_chunks(region)
            for left, right in zip(chunks, chunks[1:]):
                links.add(tuple(sorted((left, right))))
    return links


def _canonical_graph_path(chunks: Sequence[int]) -> Tuple[int, ...]:
    # Repeated boundary chunks are meaningful. For example, assemblysmall may
    # combine [123, 124] and [124, 130] into [123, 124, 124, 130]. Keeping both
    # copies lets sample-side substring matching recover the same boundary.
    path = tuple(int(value) for value in chunks)
    reverse = tuple(reversed(path))
    return min(path, reverse)


def result_within_block_paths(
    results: Mapping[str, Sequence[list]],
) -> Set[Tuple[int, ...]]:
    """Return orientation-independent graph paths of final reference blocks."""
    paths: Set[Tuple[int, ...]] = set()
    for regions in results.values():
        for region in regions:
            path = _canonical_graph_path(
                _material_region_graph_chunks(region),
            )
            # A singleton has no reference-supported internal boundary and
            # therefore cannot justify merging two sample regions.
            if len(path) > 1:
                paths.add(path)
    return paths


def result_within_block_control_paths(
    results: Mapping[str, Sequence[list]],
) -> Set[Tuple[int, ...]]:
    """Return exact final paths for same-cache cohort control replay.

    Unlike the historical cross-cache reference control, cohort learning and
    replay use the same graphDB. Retain one-base terminal chunks here because
    Assemblysmall may have explicitly absorbed that boundary and fixed replay
    must reproduce the learned result exactly.
    """
    paths: Set[Tuple[int, ...]] = set()
    for regions in results.values():
        for region in regions:
            path = _canonical_graph_path(_region_graph_chunks(region))
            if len(path) > 1:
                paths.add(path)
    return paths


def _merge_result_regions(left: list, right: list) -> list:
    merged = copy.deepcopy(left)
    for index in (0, 1):
        left_values = merged[index] if isinstance(merged[index], list) else [merged[index]]
        right_values = right[index] if isinstance(right[index], list) else [right[index]]
        merged[index] = list(left_values) + list(right_values)
    merged[3] = min(int(merged[3]), int(right[3]))
    merged[4] = max(int(merged[4]), int(right[4]))
    left_segments = merged[5] if isinstance(merged[5], list) else []
    right_segments = right[5] if isinstance(right[5], list) else []
    merged[5] = list(left_segments) + copy.deepcopy(list(right_segments))
    return merged


def apply_reference_merge_links(
    results: Mapping[str, Sequence[list]],
    merge_links: Set[Tuple[int, int]],
    block_paths: Optional[Sequence[Sequence[int]]] = None,
) -> Dict[str, List[list]]:
    """Apply graph-chunk boundary merges learned from the reference run.

    For version-4 caches, two adjacent sample region paths are merged when they
    occur as ordered, non-overlapping contiguous substrings of one final
    reference block path. Both orientations of every reference path are
    tested. After one merge the scan restarts and repeats until convergence.

    ``merge_links`` remains as the exact-boundary fallback for version-2
    caches, which did not store complete reference paths.
    """
    reference_paths: Optional[Tuple[Tuple[int, ...], ...]] = None
    if block_paths is not None:
        reference_paths = tuple(
            path for path in (
                _canonical_graph_path(value) for value in block_paths
            ) if len(path) > 1
        )

    def substring_starts(path: Tuple[int, ...], query: Tuple[int, ...]) -> List[int]:
        if not query or len(query) > len(path):
            return []
        width = len(query)
        return [
            start for start in range(len(path) - width + 1)
            if path[start:start + width] == query
        ]

    def supported_substrings(left: Sequence[int], right: Sequence[int]) -> bool:
        if reference_paths is None:
            return False
        left_path = tuple(int(value) for value in left)
        right_path = tuple(int(value) for value in right)
        if not left_path or not right_path:
            return False
        for canonical in reference_paths:
            for path in (canonical, tuple(reversed(canonical))):
                left_starts = substring_starts(path, left_path)
                right_starts = substring_starts(path, right_path)
                if any(
                    left_start + len(left_path) <= right_start
                    for left_start in left_starts
                    for right_start in right_starts
                ):
                    return True
        return False

    output: Dict[str, List[list]] = {}
    for locus, regions in results.items():
        ordered = sorted(
            (copy.deepcopy(region) for region in regions),
            key=lambda value: (int(value[3]), int(value[4])),
        )
        while True:
            changed = False
            for index in range(len(ordered) - 1):
                left = ordered[index]
                right = ordered[index + 1]
                left_chunks = _region_graph_chunks(left)
                right_chunks = _region_graph_chunks(right)
                if reference_paths is None:
                    link = _boundary_link(left, right)
                    supported = link is not None and link in merge_links
                else:
                    supported = supported_substrings(left_chunks, right_chunks)
                if supported:
                    ordered[index:index + 2] = [
                        _merge_result_regions(left, right),
                    ]
                    changed = True
                    # Restart from the beginning because the newly merged path
                    # can now match a different cached reference substring.
                    break
            if not changed:
                break
        output[locus] = ordered
    return output


def _cached_reference_merge_config(graphdb) -> Optional[Mapping[str, object]]:
    config = getattr(graphdb, "cache_metadata", {}).get(
        REFERENCE_ASSEMBLYSMALL_KEY,
    )
    if not isinstance(config, dict):
        return None
    # Versions 2 and 3 remain usable for exact-boundary replay, but do not have
    # duplicate-preserving paths for ordered substring matching.
    if int(config.get("version", 0)) not in {
        2, 3, REFERENCE_ASSEMBLYSMALL_VERSION,
    }:
        raise ValueError("unsupported cached reference assemblysmall version")
    return config


def cached_reference_merge_links(graphdb) -> Optional[Set[Tuple[int, int]]]:
    config = _cached_reference_merge_config(graphdb)
    if config is None:
        return None
    links: Set[Tuple[int, int]] = set()
    for value in config.get("within_block_links", []):
        if not isinstance(value, list) or len(value) != 2:
            raise ValueError("malformed cached reference within-block link")
        links.add(tuple(sorted((int(value[0]), int(value[1])))))
    return links


def cached_reference_block_paths(graphdb) -> Optional[Set[Tuple[int, ...]]]:
    config = _cached_reference_merge_config(graphdb)
    if config is None:
        return None
    if int(config.get("version", 0)) < 4:
        return None
    paths: Set[Tuple[int, ...]] = set()
    for value in config.get("within_block_paths", []):
        if not isinstance(value, list) or not value:
            raise ValueError("malformed cached reference within-block path")
        path = _canonical_graph_path([int(item) for item in value])
        if len(path) > 1:
            paths.add(path)
    return paths


def cached_cohort_assemblysmall_config(
    graphdb,
) -> Optional[Mapping[str, object]]:
    """Return the frozen cohort Assemblysmall control stored in a cache."""
    config = getattr(graphdb, "cache_metadata", {}).get(
        COHORT_ASSEMBLYSMALL_KEY,
    )
    if not isinstance(config, dict):
        return None
    if int(config.get("version", 0)) != COHORT_ASSEMBLYSMALL_VERSION:
        raise ValueError("unsupported cached cohort assemblysmall version")
    if config.get("scope") not in {"all_samples", "only_reference"}:
        raise ValueError("malformed cached cohort assemblysmall scope")
    return config


def cached_cohort_merge_links(graphdb) -> Optional[Set[Tuple[int, int]]]:
    config = cached_cohort_assemblysmall_config(graphdb)
    if config is None:
        return None
    links: Set[Tuple[int, int]] = set()
    for value in config.get("within_block_links", []):
        if not isinstance(value, list) or len(value) != 2:
            raise ValueError("malformed cached cohort within-block link")
        links.add(tuple(sorted((int(value[0]), int(value[1])))))
    return links


def cached_cohort_block_paths(graphdb) -> Optional[Set[Tuple[int, ...]]]:
    config = cached_cohort_assemblysmall_config(graphdb)
    if config is None:
        return None
    paths: Set[Tuple[int, ...]] = set()
    for value in config.get("within_block_paths", []):
        if not isinstance(value, list) or not value:
            raise ValueError("malformed cached cohort within-block path")
        path = _canonical_graph_path([int(item) for item in value])
        if len(path) > 1:
            paths.add(path)
    return paths


def expand_validregions_from_reference_paths(
    graphdb, block_paths: Optional[Sequence[Sequence[int]]],
) -> int:
    """Add every known graph chunk used by a final reference block.

    A valid graph interval can occur more than once in a reference alignment.
    Its query-coordinate hull may therefore contain other graph paths.  Those
    internal paths are part of the final reference block, but older caches
    recorded only the original interval in ``validregions``.  A sample lacking
    the repeated endpoint path then lost the internal paths entirely.

    Version-4 reference assemblysmall paths retain the exact chunk sequence of
    each final reference block.  Use those paths to close ``validregions`` over
    all nonnegative (real) chunks.  Negative chunks are query-only fallback
    regions and have no graph interval to add.

    The in-memory encoding is kept in the cache's declared coordinate system.
    The return value is the number of graph paths whose encoded intervals
    changed.
    """
    if not block_paths:
        return 0
    selected_chunks = {
        int(chunk)
        for path in block_paths
        for chunk in path
        if int(chunk) >= 0
    }
    if not selected_chunks:
        return 0

    metadata = getattr(graphdb, "cache_metadata", {})
    coordinate_system = metadata.get("coordinate_system", "")
    half_open = coordinate_system == "0-based-half-open"
    path_lengths = metadata.get("graph_path_lengths", {})
    if not isinstance(path_lengths, Mapping):
        path_lengths = {}

    original_intervals: Dict[str, List[Tuple[int, int]]] = {}
    for path, existing in graphdb.validregions.items():
        if not existing:
            continue
        if len(existing) % 2:
            raise ValueError(
                f"cached validregions for {path!r} has an odd coordinate count"
            )
        decoded = []
        for index in range(0, len(existing), 2):
            start = int(existing[index])
            end = int(existing[index + 1])
            if not half_open:
                end += 1
            if end > start:
                decoded.append((start, end))
        if decoded:
            original_intervals[str(path)] = decoded

    additions = cl.defaultdict(list)
    for (path, start, end), chunk_index in graphdb.blockset.items():
        if int(chunk_index) not in selected_chunks:
            continue
        start, end = int(start), int(end)
        # Folder-cache v3 inherited inclusive terminal breakpoints from
        # gfixbreaks while declaring validregions half-open.  Recover the last
        # base when the selected chunk ends at path_length - 1.  Native
        # half-open chunks already end at path_length and remain unchanged.
        path_length = path_lengths.get(str(path))
        if path_length is not None:
            path_length = int(path_length)
            if end == path_length - 1:
                end += 1
            # gfixbreaks historically made the terminal chunk end at
            # path_length + 1.  A graph interval must never include that
            # synthetic coordinate: doing so makes a terminal insertion look
            # as though it lies inside the graph path.
            start = max(0, min(start, path_length))
            end = max(0, min(end, path_length))
        if end <= start:
            continue

        existing_intervals = original_intervals.get(str(path), ())
        overlap = sum(
            max(0, min(end, existing_end) - max(start, existing_start))
            for existing_start, existing_end in existing_intervals
        )
        # A one-base overlap is the legacy inclusive/half-open boundary
        # artifact described in _material_region_graph_chunks.  It cannot
        # justify promoting the rest of a coarse graph chunk.  Disjoint chunks
        # are retained because they can be genuine internal paths enclosed by
        # a repeated reference alignment.
        if existing_intervals and 0 < overlap <= 1:
            continue
        additions[str(path)].append((start, end))

    changed = 0
    for path, new_intervals in additions.items():
        existing = graphdb.validregions.get(path, [])
        decoded = list(original_intervals.get(path, ()))
        merged: List[List[int]] = []
        for start, end in sorted([*decoded, *new_intervals]):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        encoded = [
            coordinate
            for start, end in merged
            for coordinate in (start, end if half_open else end - 1)
        ]
        if list(existing) != encoded:
            graphdb.validregions[path] = encoded
            changed += 1
    return changed


def write_reference_merge_links(
    cache_path: str,
    cutoff: int,
    merge_links: Set[Tuple[int, int]],
    block_paths: Optional[Set[Tuple[int, ...]]] = None,
) -> None:
    """Atomically add reference assemblysmall decisions to one graph cache."""
    with open(cache_path, "rt") as handle:
        data = json.load(handle)
    metadata = data.setdefault("_minsetref", {})
    metadata[REFERENCE_ASSEMBLYSMALL_KEY] = {
        "version": REFERENCE_ASSEMBLYSMALL_VERSION,
        "cutoff": int(cutoff),
        "within_block_links": [list(link) for link in sorted(merge_links)],
        "within_block_paths": [
            list(path) for path in sorted(block_paths or set())
        ],
    }
    temporary = cache_path + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            json.dump(data, output, sort_keys=True, separators=(",", ":"))
            output.write("\n")
        os.chmod(temporary, os.stat(cache_path).st_mode)
        os.replace(temporary, cache_path)
    finally:
        if os.path.exists(temporary):
            os.remove(temporary)


def block_one_alignment(
    module: ModuleType,
    graphdb,
    graph_name: str,
    input_index: int,
    row: OutputRow,
    merge_small: int,
    reference_merge_links: Optional[Set[Tuple[int, int]]] = None,
    reference_block_paths: Optional[Set[Tuple[int, ...]]] = None,
    learn_reference_merges: bool = False,
    learned_reference_paths: Optional[Set[Tuple[int, ...]]] = None,
    learned_control_paths: Optional[Set[Tuple[int, ...]]] = None,
    reference_query_intervals: Optional[Sequence[Tuple[int, int]]] = None,
    fixed_merge_control: bool = False,
) -> Tuple[List[BlockedRow], Set[Tuple[int, int]]]:
    if not row.graph_path or row.graph_path == "*":
        return [], set()
    # Region membership comes only from valid graph-interval projection.
    # overlap_genes is intentionally not consulted here: known graph genes
    # determine representation, never whether aligned query sequence exists.
    if reference_query_intervals is not None and not reference_query_intervals:
        return [], set()
    results = lossless_projection_results(
        row, graphdb, module, reference_query_intervals,
    )
    if fixed_merge_control and reference_merge_links is None:
        raise ValueError("fixed merge control requires cached merge links")
    if learn_reference_merges or reference_merge_links is None:
        changed = True
        while changed:
            results, changed_loci = module.assemblysmall(results, merge_small)
            changed = bool(changed_loci)
    else:
        results = apply_reference_merge_links(
            results, reference_merge_links, reference_block_paths,
        )
        if not fixed_merge_control:
            # Historical reference-cache mode replays reference paths and then
            # permits a sample-specific Assemblysmall pass. Cohort consistency
            # uses fixed_merge_control=True and deliberately skips this pass.
            changed = True
            while changed:
                results, changed_loci = module.assemblysmall(
                    results, merge_small, only_unrepresented=False,
                )
                changed = bool(changed_loci)
    learned_links = (
        result_within_block_links(results)
        if learn_reference_merges else set()
    )
    if learn_reference_merges and learned_reference_paths is not None:
        learned_reference_paths.update(result_within_block_paths(results))
    if learn_reference_merges and learned_control_paths is not None:
        learned_control_paths.update(
            result_within_block_control_paths(results)
        )

    query_length = row.end - row.start
    prefix = graph_prefix(graph_name)
    blocked = []
    seen = set()
    regions = sorted(
        results.get(row.query_name, ()),
        key=lambda value: (int(value[3]), int(value[4])),
    )
    for region_index, region in enumerate(regions, 1):
        local_start, local_end = sorted((int(region[3]), int(region[4])))
        if not (0 <= local_start < local_end <= query_length):
            raise ValueError(
                f"hotspot {row.hotspot_index} {row.query_name}: cached block "
                f"{local_start}-{local_end} exceeds query length {query_length}"
            )
        if row.strand == "+":
            start, end = row.start + local_start, row.start + local_end
        else:
            start, end = row.end - local_end, row.end - local_start
        key = (start, end)
        if key in seen:
            continue
        seen.add(key)
        blocked.append(BlockedRow(
            input_index, region_index, row.hotspot_index, prefix,
            row.contig, row.strand, start, end,
        ))
    return blocked, learned_links


def _init_worker(
    graphfixbreaks: str,
    merge_small: int,
    update_cache_merges: bool = False,
    reference_haplotypes: Sequence[str] = ("CHM13_h1",),
) -> None:
    global _GFIXBREAKS, _MERGE_SMALL, _UPDATE_CACHE_MERGES
    global _REFERENCE_HAPLOTYPES
    gc.disable()
    if not _ROWS_BY_GRAPH:
        raise RuntimeError("shared alignment-row dictionary is unavailable")
    _GFIXBREAKS = load_graphfixbreaks(graphfixbreaks)
    _MERGE_SMALL = merge_small
    _UPDATE_CACHE_MERGES = bool(update_cache_merges)
    _REFERENCE_HAPLOTYPES = set(reference_haplotypes)


def canonical_graph_target(
    graph_name: str, reference_haplotypes: Sequence[str],
) -> Optional[Tuple[str, str, int, int]]:
    """Parse ``HAPLOTYPE_CONTIG_START_END`` from a graph folder name."""
    for haplotype in sorted(reference_haplotypes, key=len, reverse=True):
        marker = "_" + haplotype + "_"
        marker_index = graph_name.rfind(marker)
        if marker_index < 0:
            continue
        suffix = graph_name[marker_index + len(marker):]
        try:
            contig, start_text, end_text = suffix.rsplit("_", 2)
            start, end = int(start_text), int(end_text)
        except (ValueError, TypeError):
            continue
        if contig and 0 <= start < end:
            return haplotype, contig, start, end
    return None


def reference_source_query_intervals(
    row: OutputRow, source_bed_path: str, graph_name: str,
) -> List[Tuple[int, int]]:
    """Return authoritative graph-source windows in oriented row coordinates."""
    if not os.path.isfile(source_bed_path):
        raise FileNotFoundError(
            f"reference source BED is missing: {source_bed_path}"
        )
    target = canonical_graph_target(graph_name, tuple(_REFERENCE_HAPLOTYPES))
    intervals: List[Tuple[int, int]] = []
    target_coverage: List[Tuple[int, int]] = []
    with open(source_bed_path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 7:
                fields = raw.split()
            if len(fields) < 7:
                raise ValueError(
                    f"{source_bed_path}:{line_number}: expected at least "
                    "seven BED columns"
                )
            if fields[6] not in _REFERENCE_HAPLOTYPES:
                continue
            if target is not None:
                target_haplotype, target_contig, target_start, target_end = target
                if fields[6] != target_haplotype or fields[0] != target_contig:
                    continue
            if fields[0] != row.contig:
                continue
            try:
                source_start, source_end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(
                    f"{source_bed_path}:{line_number}: invalid BED coordinates"
                ) from error
            if target is not None:
                source_start = max(source_start, target_start)
                source_end = min(source_end, target_end)
                if source_end > source_start:
                    target_coverage.append((source_start, source_end))
            absolute_start = max(row.start, source_start)
            absolute_end = min(row.end, source_end)
            if absolute_end <= absolute_start:
                continue
            if row.strand == "+":
                local_start = absolute_start - row.start
                local_end = absolute_end - row.start
            else:
                local_start = row.end - absolute_end
                local_end = row.end - absolute_start
            intervals.append((local_start, local_end))
    if target is not None and row.contig == target[1]:
        merged_coverage: List[List[int]] = []
        for start, end in sorted(set(target_coverage)):
            if merged_coverage and start <= merged_coverage[-1][1]:
                merged_coverage[-1][1] = max(merged_coverage[-1][1], end)
            else:
                merged_coverage.append([start, end])
        if merged_coverage != [[target[2], target[3]]]:
            raise ValueError(
                f"{source_bed_path}: canonical graph target "
                f"{target[1]}:{target[2]}-{target[3]} is not fully represented "
                "by its reference BED rows"
            )
    merged: List[List[int]] = []
    for start, end in sorted(set(intervals)):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _process_task(
    task: GraphTask,
) -> Tuple[List[BlockedRow], int, int, int, int]:
    if _GFIXBREAKS is None:
        raise RuntimeError("gfixbreaks worker was not initialized")
    rows = _ROWS_BY_GRAPH.get(task.graph_name, ())
    if len(rows) != task.row_count:
        raise RuntimeError(f"{task.graph_name}: shared row count changed")
    graphdb = _GFIXBREAKS.graphDB.load_json(task.cache_path)
    reference_merge_links = (
        None if _UPDATE_CACHE_MERGES
        else cached_reference_merge_links(graphdb)
    )
    reference_block_paths = (
        None if _UPDATE_CACHE_MERGES
        else cached_reference_block_paths(graphdb)
    )
    if reference_block_paths is not None:
        expand_validregions_from_reference_paths(
            graphdb, reference_block_paths,
        )
    output = []
    empty = 0
    learned_links: Set[Tuple[int, int]] = set()
    learned_paths: Set[Tuple[int, ...]] = set()
    source_excluded = 0
    for input_index, row in rows:
        try:
            source_intervals = (
                reference_source_query_intervals(
                    row, task.source_bed_path, task.graph_name,
                )
                if _UPDATE_CACHE_MERGES else None
            )
            blocked, row_learned_links = block_one_alignment(
                _GFIXBREAKS, graphdb, task.graph_name,
                input_index, row, _MERGE_SMALL,
                reference_merge_links=reference_merge_links,
                reference_block_paths=reference_block_paths,
                learn_reference_merges=_UPDATE_CACHE_MERGES,
                learned_reference_paths=learned_paths,
                reference_query_intervals=source_intervals,
            )
        except Exception as error:
            raise ValueError(
                f"graph {task.graph_name}, input row {input_index + 1}, "
                f"hotspot {row.hotspot_index}, query {row.query_name}, "
                f"{row.contig}:{row.start}-{row.end}{row.strand}: {error}"
            ) from error
        output.extend(blocked)
        learned_links.update(row_learned_links)
        empty += not blocked
        source_excluded += int(_UPDATE_CACHE_MERGES and not source_intervals)
    cache_updated = 0
    if _UPDATE_CACHE_MERGES:
        write_reference_merge_links(
            task.cache_path, _MERGE_SMALL, learned_links, learned_paths,
        )
        cache_updated = 1
    used_reference_cache = int(
        not _UPDATE_CACHE_MERGES and reference_merge_links is not None
    )
    return output, empty, cache_updated, used_reference_cache, source_excluded


def load_tasks(
    alignment: str,
    graph_list: str,
    graph_root: str,
    strict_missing_cache: bool,
    only_reference: bool = False,
    reference_haplotypes: Sequence[str] = ("CHM13_h1",),
) -> Tuple[List[GraphTask], int, int]:
    global _ROWS_BY_GRAPH
    graph_root = os.path.abspath(os.path.expanduser(graph_root))
    if not os.path.isdir(graph_root):
        raise NotADirectoryError(f"graph folder not found: {graph_root}")
    graph_names = read_graph_names(graph_list)
    grouped: Dict[str, List[Tuple[int, OutputRow]]] = cl.defaultdict(list)
    rows = select_alignment_rows(
        read_alignment_output(alignment),
        reference_haplotypes, only_reference,
    )
    for input_index, row in enumerate(rows):
        if row.hotspot_index > len(graph_names):
            raise IndexError(
                f"alignment hotspot index {row.hotspot_index} exceeds "
                f"graph-list length {len(graph_names)}"
            )
        graph_name = graph_names[row.hotspot_index - 1]
        grouped[graph_name].append((input_index, row))

    tasks = []
    missing = []
    for hotspot_index, graph_name in enumerate(graph_names, 1):
        graph_rows = grouped.get(graph_name)
        if not graph_rows:
            continue
        cache_path = cache_path_for_graph(graph_root, graph_name)
        if not os.path.isfile(cache_path):
            missing.append((graph_name, cache_path, len(graph_rows)))
            continue
        tasks.append(GraphTask(
            hotspot_index,
            graph_name,
            cache_path,
            os.path.join(graph_root, graph_name, graph_name + ".bed"),
            len(graph_rows),
        ))
    if missing:
        examples = ", ".join(value[0] for value in missing[:5])
        message = (
            f"{len(missing)} graph caches covering "
            f"{sum(value[2] for value in missing)} alignment rows are missing; "
            f"first: {examples}"
        )
        if strict_missing_cache:
            raise FileNotFoundError(message)
        LOG.warning("%s; those rows will be skipped", message)
    task_names = {task.graph_name for task in tasks}
    _ROWS_BY_GRAPH = {
        name: tuple(values) for name, values in grouped.items()
        if name in task_names
    }
    return tasks, len(rows), sum(value[2] for value in missing)


def write_output(
    path: str, rows: Iterable[BlockedRow], include_header: bool,
) -> None:
    absolute = os.path.abspath(path)
    Path(os.path.dirname(absolute) or ".").mkdir(parents=True, exist_ok=True)
    temporary = absolute + f".tmp.{os.getpid()}"
    try:
        with open(temporary, "wt") as output:
            if include_header:
                output.write("chrom\tstart\tend\tname\tscore\tstrand\n")
            for row in sorted(
                rows, key=lambda value: (value.input_index, value.region_index),
            ):
                output.write(row.line())
        os.replace(temporary, absolute)
    finally:
        try:
            os.remove(temporary)
        except FileNotFoundError:
            pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-i", "--input", required=True,
        help="existing ten-column align_partition_hotspots.py output",
    )
    parser.add_argument(
        "-G", "--graph-folder", required=True,
        help="root containing GRAPH_NAME/GRAPH_NAMEcache.json",
    )
    parser.add_argument(
        "-L", "--graph-list", required=True,
        help="ordered graph list used to create the input alignment",
    )
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument(
        "--graphfixbreaks",
        default=os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "gfixbreaks.py",
        ),
        help="gfixbreaks.py implementation (default: beside this script)",
    )
    parser.add_argument(
        "--merge-small", type=int, default=DEFAULT_MERGE_SMALL,
        help="assemblysmall merge cutoff (default: 20000)",
    )
    parser.add_argument(
        "--update-cache-merges",
        action="store_true",
        help=(
            "reference run: learn assemblysmall boundary merges and persist "
            "them in each graph cache for later samples"
        ),
    )
    parser.add_argument(
        "--reference-haplotype",
        action="append",
        help=(
            "reference haplotype accepted from each graph .bed when "
            "--update-cache-merges is used; repeat for aliases "
            "(default: CHM13_h1)"
        ),
    )
    parser.add_argument(
        "--onlyref", "--only-ref", action="store_true",
        help=(
            "write blocks only for named reference alignment rows; "
            "default: block all samples"
        ),
    )
    parser.add_argument("-j", "--jobs", type=int, default=16)
    parser.add_argument(
        "--strict-missing-cache", action="store_true",
        help="fail instead of warning/skipping alignment rows without a cache",
    )
    parser.add_argument(
        "--header", action="store_true",
        help="write a BED6 column header",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.merge_small < 0:
        parser.error("--merge-small must be nonnegative")
    return args


def run(args: argparse.Namespace) -> int:
    reference_haplotypes = tuple(
        args.reference_haplotype or ("CHM13_h1",)
    )
    tasks, input_rows, missing_rows = load_tasks(
        args.input, args.graph_list, args.graph_folder,
        args.strict_missing_cache,
        args.onlyref, reference_haplotypes,
    )
    if not tasks:
        if input_rows and missing_rows == input_rows:
            raise ValueError("none of the aligned graphs has a cache.json file")
        write_output(args.output, (), args.header)
        return 0
    workers = min(args.jobs, len(tasks))
    LOG.info(
        "Blocking %d alignment rows across %d graph caches with %d worker "
        "process%s",
        input_rows - missing_rows, len(tasks), workers,
        "" if workers == 1 else "es",
    )
    all_rows: List[BlockedRow] = []
    empty_rows = 0
    updated_caches = 0
    reference_merge_caches = 0
    source_excluded_rows = 0
    if workers == 1:
        _init_worker(
            args.graphfixbreaks, args.merge_small, args.update_cache_merges,
            reference_haplotypes,
        )
        results = (_process_task(task) for task in tasks)
        for (
            blocked, empty, cache_updated, used_reference_cache,
            source_excluded,
        ) in results:
            all_rows.extend(blocked)
            empty_rows += empty
            updated_caches += cache_updated
            reference_merge_caches += used_reference_cache
            source_excluded_rows += source_excluded
    else:
        if "fork" not in multiprocessing.get_all_start_methods():
            raise RuntimeError(
                "shared alignment rows require the fork multiprocessing method"
            )
        gc.collect()
        gc.freeze()
        try:
            with ProcessPoolExecutor(
                max_workers=workers,
                mp_context=multiprocessing.get_context("fork"),
                initializer=_init_worker,
                initargs=(
                    args.graphfixbreaks,
                    args.merge_small,
                    args.update_cache_merges,
                    reference_haplotypes,
                ),
            ) as executor:
                futures = [executor.submit(_process_task, task) for task in tasks]
                report_every = max(1, len(futures) // 20)
                for completed, future in enumerate(as_completed(futures), 1):
                    (
                        blocked,
                        empty,
                        cache_updated,
                        used_reference_cache,
                        source_excluded,
                    ) = future.result()
                    all_rows.extend(blocked)
                    empty_rows += empty
                    updated_caches += cache_updated
                    reference_merge_caches += used_reference_cache
                    source_excluded_rows += source_excluded
                    if completed % report_every == 0 or completed == len(futures):
                        LOG.info(
                            "Cached blocking progress: %d/%d graphs",
                            completed, len(futures),
                        )
        finally:
            gc.unfreeze()
    write_output(args.output, all_rows, args.header)
    if args.update_cache_merges:
        LOG.info(
            "Persisted reference assemblysmall decisions in %d graph cache%s",
            updated_caches,
            "" if updated_caches == 1 else "s",
        )
        if source_excluded_rows:
            LOG.info(
                "Excluded %d reference alignment row%s with no authoritative "
                "reference source interval",
                source_excluded_rows,
                "" if source_excluded_rows == 1 else "s",
            )
    else:
        LOG.info(
            "Applied persisted reference assemblysmall decisions from %d/%d "
            "graph caches",
            reference_merge_caches,
            len(tasks),
        )
        if reference_merge_caches < len(tasks):
            LOG.warning(
                "%d graph cache%s lacked reference assemblysmall metadata and "
                "used legacy sample-specific merging; rerun the reference with "
                "--update-cache-merges",
                len(tasks) - reference_merge_caches,
                "" if len(tasks) - reference_merge_caches == 1 else "s",
            )
    LOG.info(
        "Wrote %d BED6 blocks to %s; %d aligned row%s "
        "had no query projection inside a valid graph interval",
        len(all_rows), args.output, empty_rows,
        "" if empty_rows == 1 else "s",
    )
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
    )
    try:
        return run(args)
    except Exception as error:
        LOG.error("%s", error)
        if args.verbose:
            LOG.exception("cached alignment blocking failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
