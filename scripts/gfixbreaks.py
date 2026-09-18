#!/usr/bin/env python3

import copy
import os
import threading
import argparse
import concurrent.futures
import datetime
import collections as cl
import re
import glob
import json
import hashlib 
import pickle
import sqlite3
import tempfile
import bisect
from statistics import median

script_folder = os.path.dirname(os.path.abspath(__file__))
cutoffdistance = 30000
mergesmall = 20000
cache = ""
OVERLAP_PASSED_STATE_LIMIT = max(
	1, int(os.environ.get("GFIXBREAKS_PASSED_STATE_LIMIT", "10000")),
)
OVERLAP_DEDUP_STATE_LIMIT = 50000
PATH_COORDINATE_SUFFIX_RE = re.compile(r"_(\d+)_(\d+)$")
REFERENCE_ASSEMBLYSMALL_KEY = "reference_assemblysmall"
REFERENCE_ASSEMBLYSMALL_VERSION = 4


def graph_path_size(path, cigar=""):
	"""Return graph-path length from its coordinate suffix or full CIGAR.

	Reference paths normally end in ``_START_END``.  Alternative graph paths
	can instead end in an assembly contig identifier.  Whole-graph CIGARs still
	encode the complete path length: M/= /X, D, and H consume reference bases.
	"""
	match = PATH_COORDINATE_SUFFIX_RE.search(path)
	if match is not None:
		start, end = map(int, match.groups())
		if end > start:
			return end - start
	reference_size = sum(
		int(size)
		for size, operation in re.findall(r"(\d+)([M=XHDN])", cigar)
	)
	if reference_size > 0:
		return reference_size
	raise ValueError(
		f"graph path {path!r} has no _START_END suffix and its CIGAR "
		"does not encode a positive reference length"
	)


def normalize_validregions(values, coordinate_system=""):
	"""Return cached graph intervals as sorted 0-based half-open pairs.

	Legacy gfixbreaks caches stored inclusive right endpoints.  New minsetref
	caches declare ``0-based-half-open`` in ``_minsetref.coordinate_system``.
	"""
	if not values:
		return []
	if isinstance(values[0], (int, float, str)):
		if len(values) % 2:
			raise ValueError("validregions contains an odd coordinate count")
		pairs = [values[index:index + 2] for index in range(0, len(values), 2)]
	else:
		pairs = values
	intervals = []
	for value in pairs:
		if len(value) < 2:
			continue
		start, end = sorted((int(value[0]), int(value[1])))
		if coordinate_system != "0-based-half-open":
			end += 1
		if end > start:
			intervals.append((start, end))
	merged = []
	for start, end in sorted(intervals):
		if merged and start <= merged[-1][1]:
			merged[-1][1] = max(merged[-1][1], end)
		else:
			merged.append([start, end])
	return [tuple(value) for value in merged]


def _cigar_operations(cigar):
	"""Parse the sequence-free graph CIGAR representation."""
	operations = []
	cursor = 0
	for match in re.finditer(r"(\d+)([M=XHDIN])", cigar):
		if match.start() != cursor:
			raise ValueError(f"invalid graph CIGAR {cigar!r}")
		operations.append((int(match.group(1)), match.group(2)))
		cursor = match.end()
	if not operations or cursor != len(cigar):
		raise ValueError(f"invalid graph CIGAR {cigar!r}")
	return operations


def project_validregions_to_query(
    spans, cigarstr, validregions, coordinate_system="",
    occurrence_aware=False,
):
    """Project valid graph intervals onto one query alignment.

    ``occurrence_aware=True`` keeps separate visits to the same path apart.
    The legacy default still forms one hull across all visits for reference
    cache closure; blocking enables occurrence-aware mode.
    """
    cigar_tokens = re.findall(r"[><][^><]+", cigarstr)
    if len(cigar_tokens) != len(spans):
        raise ValueError(
            "graph path/CIGAR component counts differ while projecting "
            "valid regions"
        )

    components = []
    for span, token in zip(spans, cigar_tokens):
        if ":" not in token:
            raise ValueError(f"invalid graph CIGAR component {token!r}")
        components.append((span, token.split(":", 1)[1]))

    output = []
    for path, values in validregions.items():
        intervals = normalize_validregions(values, coordinate_system)
        if not intervals:
            continue
        if occurrence_aware:
            groups = [
                (index, [(span, cigar)])
                for index, (span, cigar) in enumerate(components)
                if span[0][1:] == path
            ]
        else:
            groups = [
                (-1, [item for item in components if item[0][0][1:] == path])
            ]

        for component_index, path_components in groups:
            if not path_components:
                continue
            for graph_start, graph_end in intervals:
                pieces = []
                orientations = []
                for span, cigar in path_components:
                    marker = span[0][0]
                    path_length = graph_path_size(path, cigar)
                    ref_cursor = 0 if marker == ">" else path_length
                    query_cursor = min(map(int, span[2]))
                    for size, operation in _cigar_operations(cigar):
                        consume_graph = operation in {"M", "=", "X", "H", "D", "N"}
                        consume_query = operation in {"M", "=", "X", "I"}
                        if operation == "I" and graph_start <= ref_cursor < graph_end:
                            pieces.append((query_cursor, query_cursor + size))
                            orientations.append(marker)
                        if operation in {"M", "=", "X"}:
                            if marker == ">":
                                op_start, op_end = ref_cursor, ref_cursor + size
                                left = max(graph_start, op_start)
                                right = min(graph_end, op_end)
                                if right > left:
                                    qstart = query_cursor + left - op_start
                                    qend = query_cursor + right - op_start
                            else:
                                op_start, op_end = ref_cursor - size, ref_cursor
                                left = max(graph_start, op_start)
                                right = min(graph_end, op_end)
                                if right > left:
                                    qstart = query_cursor + op_end - right
                                    qend = query_cursor + op_end - left
                            if right > left:
                                pieces.append((int(qstart), int(qend)))
                                orientations.append(marker)
                        if consume_graph:
                            ref_cursor += size if marker == ">" else -size
                        if consume_query:
                            query_cursor += size
                if not pieces:
                    continue
                pieces = sorted(set(pieces))
                output.append({
                    "path": path,
                    "graph_start": graph_start,
                    "graph_end": graph_end,
                    "query_start": min(value[0] for value in pieces),
                    "query_end": max(value[1] for value in pieces),
                    "query_pieces": pieces,
                    "orientation": orientations[0] if orientations else ">",
                    "component_index": component_index,
                })
    return sorted(output, key=lambda value: (
        value["query_start"], value["query_end"], value["path"],
        value.get("component_index", 0), value["graph_start"],
    ))


def _graph_intervals_for_query_range(cigarstr, query_start, query_end):
	"""Return every half-open graph interval traversed in a query range.

	Unlike the recovered endpoint-only cache code, this walks every graph-CIGAR
	component.  Internal paths in an ``A -> B -> A`` block are therefore kept.
	Hard clips advance the graph cursor but are not aligned graph sequence;
	deletions and reference skips inside the selected query range are retained.
	"""
	query_start, query_end = int(query_start), int(query_end)
	if query_end <= query_start:
		return {}
	output = cl.defaultdict(list)
	for token in re.findall(r"[><][^><]+", cigarstr):
		if ":" not in token:
			raise ValueError(f"invalid graph CIGAR component {token!r}")
		encoded, cigar = token.split(":", 1)
		marker = encoded[0]
		try:
			path, _aligned_graph_start, query_cursor = encoded[1:].rsplit("_", 2)
			int(_aligned_graph_start)
			query_cursor = int(query_cursor)
		except (ValueError, TypeError) as error:
			raise ValueError(
				f"invalid positioned graph CIGAR component {token!r}"
			) from error
		graph_cursor = 0 if marker == ">" else graph_path_size(path, cigar)
		for size, operation in _cigar_operations(cigar):
			consume_graph = operation in {"M", "=", "X", "H", "D", "N"}
			consume_query = operation in {"M", "=", "X", "I"}
			if operation in {"M", "=", "X"}:
				op_query_start, op_query_end = query_cursor, query_cursor + size
				left = max(query_start, op_query_start)
				right = min(query_end, op_query_end)
				if right > left:
					if marker == ">":
						graph_start = graph_cursor + left - op_query_start
						graph_end = graph_cursor + right - op_query_start
					else:
						graph_start = graph_cursor - (right - op_query_start)
						graph_end = graph_cursor - (left - op_query_start)
					output[path].append((graph_start, graph_end))
			elif (
				operation in {"D", "N"}
				and query_start <= query_cursor < query_end
			):
				if marker == ">":
					graph_start, graph_end = graph_cursor, graph_cursor + size
				else:
					graph_start, graph_end = graph_cursor - size, graph_cursor
				output[path].append((graph_start, graph_end))
			if consume_graph:
				graph_cursor += size if marker == ">" else -size
			if consume_query:
				query_cursor += size
	return output


def _merge_half_open_graph_intervals(intervals):
	merged = []
	for start, end in sorted(
		(min(int(start), int(end)), max(int(start), int(end)))
		for start, end in intervals
		if int(start) != int(end)
	):
		if merged and start <= merged[-1][1]:
			merged[-1][1] = max(merged[-1][1], end)
		else:
			merged.append([start, end])
	return merged


def build_total_validregions(regions, original_loci, graphinfo, contigspan):
	"""Build the complete half-open graph closure of reference blocks.

	Initial graph spans come from every aligned component inside each final
	reference block.  If a valid span is scattered/repeated, its complete query
	hull is then projected back to *all* graph components inside that hull.
	This repeats to convergence, ensuring that internal novel paths are recorded
	in the cache and can be recovered in later samples that lack the repeated
	endpoint path.
	"""
	graph_spans = cl.defaultdict(list)

	def resolve_locus(value):
		if value in graphinfo:
			return value
		if value.startswith("Ref_") and value[4:] in graphinfo:
			return value[4:]
		prefixed = "Ref_" + value
		if prefixed in graphinfo:
			return prefixed
		raise KeyError(value)

	for region in regions:
		name = str(region[0])
		locus, local_start, local_end = original_loci[name]
		locus = resolve_locus(str(locus))
		# annotate_regions() exposes the recovered legacy inclusive right end.
		for path, spans in _graph_intervals_for_query_range(
			graphinfo[locus], int(local_start), int(local_end) + 1,
		).items():
			graph_spans[path].extend(spans)

	graph_spans = {
		path: _merge_half_open_graph_intervals(spans)
		for path, spans in graph_spans.items()
		if spans
	}
	while graph_spans:
		before = tuple(
			(path, tuple(tuple(span) for span in spans))
			for path, spans in sorted(graph_spans.items())
		)
		encoded = {
			path: [coordinate for span in spans for coordinate in span]
			for path, spans in graph_spans.items()
		}
		additions = cl.defaultdict(list)
		for locus, cigarstr in graphinfo.items():
			spans = contigspan.get(locus)
			if spans is None and locus.startswith("Ref_"):
				spans = contigspan.get(locus[4:])
			if spans is None:
				spans = contigspan.get("Ref_" + locus)
			if not spans:
				continue
			for projection in project_validregions_to_query(
				spans, cigarstr, encoded, "0-based-half-open",
			):
				for path, projected in _graph_intervals_for_query_range(
					cigarstr,
					projection["query_start"], projection["query_end"],
				).items():
					additions[path].extend(projected)
		for path, spans in additions.items():
			graph_spans[path] = _merge_half_open_graph_intervals([
				*graph_spans.get(path, ()), *spans,
			])
		after = tuple(
			(path, tuple(tuple(span) for span in spans))
			for path, spans in sorted(graph_spans.items())
		)
		if after == before:
			break

	return {
		path: [coordinate for span in spans for coordinate in span]
		for path, spans in graph_spans.items()
	}


def graph_chunks_for_interval(graphdb, path, start, end, orientation=">"):
	"""Return cached graph chunks intersecting one valid graph interval."""
	chunks = sorted(
		(
			int(chunk_start), int(chunk_end), int(chunk_index),
		)
		for (chunk_path, chunk_start, chunk_end), chunk_index
		in graphdb.blockset.items()
		if chunk_path == path and int(chunk_start) < end and int(chunk_end) > start
	)
	values = [value[2] for value in chunks]
	return values if orientation == ">" else list(reversed(values))


def lossless_validregion_results(
	graphdb, graphinfo, contigspan, occurrence_aware=False,
	query_valid_intervals=None,
):
	"""Retain known block representation and add every uncovered query span.

	``overlap_genes`` remains responsible for the consistent reference block
	representation.  Its result is clipped to query projections of valid graph
	intervals, but neither unrecognized fallback chunks nor uncovered query
	sequence are discarded.  Any leftover portion becomes an independent block
	that ``assemblysmall`` may subsequently merge.
	"""
	coordinate_system = graphdb.cache_metadata.get("coordinate_system", "")
	results = cl.defaultdict(list)
	unknown_index = -1
	for locus in contigspan:
		projections = project_validregions_to_query(
			contigspan[locus], graphinfo[locus], graphdb.validregions,
			coordinate_system, occurrence_aware=occurrence_aware,
		)
		allowed = sorted((
			(int(value["query_start"]), int(value["query_end"]),
			 int(value.get("component_index", -1)))
			for value in projections
			if int(value["query_end"]) > int(value["query_start"])
		))
		# A reference graph interval may be visited more than once by the same
		# whole alignment.  Graph validity alone therefore cannot identify which
		# occurrence belongs to the local graph's authoritative source window.
		# Reference blocking supplies that window in query coordinates.  Apply it
		# before Assemblysmall so an edge repeat can neither become a block nor
		# teach a cached reference merge outside the source interval.  Sample
		# blocking leaves query_valid_intervals unset and retains every occurrence.
		if query_valid_intervals is not None:
			masks = sorted(
				(min(int(value[0]), int(value[1])),
				 max(int(value[0]), int(value[1])))
				for value in query_valid_intervals.get(locus, ())
				if int(value[0]) != int(value[1])
			)
			clipped_allowed = []
			for start, end, component_index in allowed:
				for mask_start, mask_end in masks:
					left = max(start, mask_start)
					right = min(end, mask_end)
					if right > left:
						clipped_allowed.append((
							left, right, component_index,
						))
			allowed = sorted(set(clipped_allowed))
		merged_allowed = []
		for start, end, component_index in allowed:
			if (
				merged_allowed
				and (
					component_index == merged_allowed[-1][2]
					# Adjacent query pieces form one continuous block even when
					# the graph traversal is split into consecutive components.
					or start == merged_allowed[-1][1]
				)
				and start <= merged_allowed[-1][1]
			):
				merged_allowed[-1][1] = max(merged_allowed[-1][1], end)
			else:
				merged_allowed.append([start, end, component_index])

		raw_regions = graphdb.overlap_genes(
			contigspan[locus], graphinfo[locus],
		)
		for region in raw_regions:
			if len(region) < 6:
				continue
			subsegments = region[5] or [[region[3], region[4]]]
			query_chunks = (
				region[0] if isinstance(region[0], list)
				else [region[0]] * len(subsegments)
			)
			graph_chunks = (
				region[1] if isinstance(region[1], list)
				else [region[1]] * len(subsegments)
			)
			for valid_start, valid_end, _component_index in merged_allowed:
				kept_ranges = []
				kept_query_chunks = []
				kept_graph_chunks = []
				for index, subsegment in enumerate(subsegments):
					if len(subsegment) < 2:
						continue
					start, end = sorted(map(int, subsegment[:2]))
					left = max(start, valid_start)
					right = min(end, valid_end)
					if right <= left:
						continue
					kept_ranges.append([left, right])
					if index < len(query_chunks):
						kept_query_chunks.append(query_chunks[index])
					if index < len(graph_chunks):
						kept_graph_chunks.append(graph_chunks[index])
				if not kept_ranges:
					continue
				kept = copy.deepcopy(region)
				kept[0] = kept_query_chunks
				kept[1] = kept_graph_chunks
				kept[3] = min(value[0] for value in kept_ranges)
				kept[4] = max(value[1] for value in kept_ranges)
				kept[5] = kept_ranges
				results[locus].append(kept)

		# Region hulls, rather than only matched subsegments, are covered.  This
		# retains novel inserted sequence already enclosed by a known block.
		covered = sorted(
			[int(region[3]), int(region[4])]
			for region in results[locus]
			if int(region[4]) > int(region[3])
		)
		merged_covered = []
		for start, end in covered:
			if merged_covered and start <= merged_covered[-1][1]:
				merged_covered[-1][1] = max(merged_covered[-1][1], end)
			else:
				merged_covered.append([start, end])
		for valid_start, valid_end, _component_index in merged_allowed:
			cursor = valid_start
			for cover_start, cover_end in merged_covered:
				if cover_end <= cursor or cover_start >= valid_end:
					continue
				if cover_start > cursor:
					results[locus].append([
						[0], [unknown_index], ("unrepresented",),
						cursor, min(cover_start, valid_end),
						[[cursor, min(cover_start, valid_end)]],
					])
					unknown_index -= 1
				cursor = max(cursor, cover_end)
				if cursor >= valid_end:
					break
			if cursor < valid_end:
				results[locus].append([
					[0], [unknown_index], ("unrepresented",),
					cursor, valid_end, [[cursor, valid_end]],
				])
				unknown_index -= 1
		results[locus].sort(key=lambda value: (value[3], value[4]))
	return results


def _result_graph_chunks(region):
	value = region[1]
	return [int(item) for item in value] if isinstance(value, list) else [int(value)]


def _canonical_graph_path(chunks):
	path = tuple(int(value) for value in chunks)
	reverse = tuple(reversed(path))
	return min(path, reverse)


def reference_assemblysmall_config(results, cutoff):
	"""Serialize the reference's final assemblysmall representation."""
	links = set()
	paths = set()
	for regions in results.values():
		for region in regions:
			chunks = _result_graph_chunks(region)
			if any(value < 0 for value in chunks):
				continue
			links.update(tuple(sorted(pair)) for pair in zip(chunks, chunks[1:]))
			path = _canonical_graph_path(chunks)
			if len(path) > 1:
				paths.add(path)
	return {
		"version": REFERENCE_ASSEMBLYSMALL_VERSION,
		"cutoff": int(cutoff),
		"within_block_links": [list(value) for value in sorted(links)],
		"within_block_paths": [list(value) for value in sorted(paths)],
	}


def _merge_result_regions(left, right):
	merged = copy.deepcopy(left)
	for index in (0, 1):
		left_values = merged[index] if isinstance(merged[index], list) else [merged[index]]
		right_values = right[index] if isinstance(right[index], list) else [right[index]]
		merged[index] = list(left_values) + list(right_values)
	merged[3] = min(int(merged[3]), int(right[3]))
	merged[4] = max(int(merged[4]), int(right[4]))
	merged[5] = list(merged[5]) + copy.deepcopy(list(right[5]))
	return merged


def apply_reference_assemblysmall(results, config):
	"""Replay cached reference assemblysmall paths on sample blocks."""
	version = int(config.get("version", 0))
	if version not in {2, 3, REFERENCE_ASSEMBLYSMALL_VERSION}:
		raise ValueError(f"unsupported reference assemblysmall version {version}")
	links = {
		tuple(sorted((int(value[0]), int(value[1]))))
		for value in config.get("within_block_links", []) if len(value) == 2
	}
	paths = None
	if version >= 4:
		paths = {
			_canonical_graph_path(value)
			for value in config.get("within_block_paths", []) if value
		}

	def starts(path, query):
		return [
			index for index in range(len(path) - len(query) + 1)
			if path[index:index + len(query)] == query
		] if query and len(query) <= len(path) else []

	def supported(left, right):
		left_path = tuple(_result_graph_chunks(left))
		right_path = tuple(_result_graph_chunks(right))
		if any(value < 0 for value in left_path + right_path):
			return False
		if paths is None:
			return tuple(sorted((left_path[-1], right_path[0]))) in links
		for canonical in paths:
			for path in (canonical, tuple(reversed(canonical))):
				if any(
					left_start + len(left_path) <= right_start
					for left_start in starts(path, left_path)
					for right_start in starts(path, right_path)
				):
					return True
		return False

	output = {}
	for locus, regions in results.items():
		ordered = sorted(copy.deepcopy(list(regions)), key=lambda value: (value[3], value[4]))
		while True:
			for index in range(len(ordered) - 1):
				if supported(ordered[index], ordered[index + 1]):
					ordered[index:index + 2] = [
						_merge_result_regions(ordered[index], ordered[index + 1]),
					]
					break
			else:
				break
		output[locus] = ordered
	return output


class _OverlapState:
	"""Compact persistent state used internally by graphDB.overlap_genes."""

	__slots__ = (
		"parent", "ichunk", "chunkindex", "involves", "qstart", "qend",
		"depth", "min_ichunk", "max_ichunk",
	)

	def __init__(self, parent, ichunk, chunkindex, involves, qstart, qend):
		self.parent = parent
		self.ichunk = ichunk
		self.chunkindex = chunkindex
		self.involves = tuple(involves)
		self.qstart = qstart
		self.qend = qend
		if parent is None:
			self.depth = 1
			self.min_ichunk = ichunk
			self.max_ichunk = ichunk
		else:
			self.depth = parent.depth + 1
			self.min_ichunk = min(parent.min_ichunk, ichunk)
			self.max_ichunk = max(parent.max_ichunk, ichunk)

	def exact_child_key(self):
		"""Key that only merges children with the identical parent/history."""
		return (
			id(self.parent), self.ichunk, self.chunkindex, self.involves,
			self.qstart, self.qend,
		)

	def max_score(self):
		return max(x[2] for x in self.involves)

	def materialize(self):
		"""Convert a selected persistent state to the legacy six-list form."""
		nodes = []
		state = self
		while state is not None:
			nodes.append(state)
			state = state.parent
		nodes.reverse()
		return [
			[x.ichunk for x in nodes],
			[x.chunkindex for x in nodes],
			list(self.involves),
			nodes[0].qstart,
			nodes[-1].qend,
			[[x.qstart, x.qend] for x in nodes],
		]


class _DPOverlapState:
	"""One polynomial-DP state for a single reference path traversal.

	Future transitions depend only on the first/last query chunks, the current
	reference path position, and the accumulated score.  ``parent`` is retained
	only to reconstruct the selected query chunks after the DP finishes.
	"""

	__slots__ = (
		"parent", "ichunk", "chunkindex", "gene", "alignindex", "score",
		"qstart", "qend", "first_ichunk", "depth",
	)

	def __init__(
		self, parent, ichunk, chunkindex, gene, alignindex, score, qstart, qend,
	):
		self.parent = parent
		self.ichunk = ichunk
		self.chunkindex = chunkindex
		self.gene = gene
		self.alignindex = alignindex
		self.score = score
		self.qstart = qstart
		self.qend = qend
		if parent is None:
			self.first_ichunk = ichunk
			self.depth = 1
		else:
			self.first_ichunk = parent.first_ichunk
			self.depth = parent.depth + 1

	def active_key(self):
		return (
			self.first_ichunk, self.ichunk, self.gene, self.alignindex,
		)

	def interval_key(self):
		return self.first_ichunk, self.ichunk

	def materialize(self):
		nodes = []
		state = self
		while state is not None:
			nodes.append(state)
			state = state.parent
		nodes.reverse()
		return [
			[x.ichunk for x in nodes],
			[x.chunkindex for x in nodes],
			[(self.gene, self.alignindex, self.score)],
			nodes[0].qstart,
			nodes[-1].qend,
			[[x.qstart, x.qend] for x in nodes],
		]


class _BestIntervalStates:
	"""Retain only the best completed state for each query interval."""

	def __init__(self):
		self.states = {}
		self.ordinal = 0

	def add(self, state):
		ordinal = self.ordinal
		self.ordinal += 1
		key = state.interval_key()
		old = self.states.get(key)
		if old is None or state.score > old[0]:
			self.states[key] = (state.score, ordinal, state)

	def extend(self, states):
		for state in states:
			self.add(state)

	def iter_sorted(self):
		for _score, _ordinal, state in sorted(
			self.states.values(), key=lambda item: (-item[0], item[1]),
		):
			yield state.materialize()


class _PassedStateStore:
	"""Bound completed overlap states in RAM, spilling excess to SQLite."""

	def __init__(self, memory_limit=None):
		if memory_limit is None:
			memory_limit = OVERLAP_PASSED_STATE_LIMIT
		self.memory_limit = max(1, int(memory_limit))
		self.states = []
		self.connection = None
		self.path = ""
		self.ordinal = 0

	def _open(self):
		if self.connection is not None:
			return
		fd, self.path = tempfile.mkstemp(
			prefix=".gfixbreaks.overlap.", suffix=".sqlite3",
		)
		os.close(fd)
		self.connection = sqlite3.connect(self.path)
		self.connection.execute("PRAGMA journal_mode=OFF")
		self.connection.execute("PRAGMA synchronous=OFF")
		self.connection.execute("PRAGMA temp_store=FILE")
		self.connection.execute("PRAGMA cache_size=-8192")
		self.connection.execute(
			"CREATE TABLE states (score INTEGER, ordinal INTEGER, state BLOB)"
		)
		for state in self.states:
			self._insert(state)
		self.states = []

	def _insert(self, state):
		payload = pickle.dumps(state.materialize(), protocol=5)
		self.connection.execute(
			"INSERT INTO states VALUES (?, ?, ?)",
			(state.max_score(), self.ordinal, sqlite3.Binary(payload)),
		)
		self.ordinal += 1

	def add(self, state):
		# determinegenes() discarded these states before sorting.
		if len(state.involves) == 0:
			return
		if self.connection is None and len(self.states) < self.memory_limit:
			self.states.append(state)
			self.ordinal += 1
			return
		if self.connection is None:
			# _insert() assigns ordinals while flushing. Rewind so their order
			# remains identical to the original passed_genes list.
			self.ordinal = 0
			self._open()
		self._insert(state)

	def extend(self, states):
		for state in states:
			self.add(state)

	def iter_sorted(self):
		if self.connection is None:
			for state in sorted(
				self.states, key=lambda x: x.max_score(), reverse=True,
			):
				yield state.materialize()
			return
		self.connection.commit()
		for payload, in self.connection.execute(
			"SELECT state FROM states ORDER BY score DESC, ordinal ASC"
		):
			yield pickle.loads(payload)

	def close(self):
		if self.connection is not None:
			self.connection.close()
			self.connection = None
		if self.path:
			try:
				os.remove(self.path)
			except FileNotFoundError:
				pass
			self.path = ""
		self.states = []


def _deduplicate_overlap_states(states):
	"""Drop only states with identical parent identity and current contents."""
	# A temporary hash table larger than this could itself become the peak-RAM
	# allocation. Skipping optional deduplication leaves results unchanged.
	if len(states) > OVERLAP_DEDUP_STATE_LIMIT:
		return states
	seen = set()
	output = []
	for state in states:
		key = state.exact_child_key()
		if key in seen:
			continue
		seen.add(key)
		output.append(state)
	return output

def process_locus2(region, haplopath, folder, outputfile):
	
	newname,  contig, start, end, ifexon, mapinfo = region
	
	if "#" not in contig: 
		haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
	else:
		haplo = "_h".join(contig.split("#")[:2])
		
	path = haplopath[haplo]
	if max(0, end) - max(0,start) <= 0:
		return
	#outputfile0 = os.path.join(folder, f"{newname}.fa")
	start = max(0,int(start))
	header = ">{}\t{}:{}-{}\t{}\t{}".format(newname, contig, start, end, ifexon, mapinfo)
	
	cmd1 = f'echo "{header}" >> {outputfile} && samtools faidx {path} {contig}:{start+1}-{end} | tail -n +2 >> {outputfile}'
	os.system(cmd1)
	

def makenewfasta2(regions, queryfile, outputfile, folder, threads):
	
	haplopath = dict()
	with open(queryfile, mode = 'r') as f:
		for line in f:
			line = line.strip().split()
			haplopath[line[0]] = line[1]
			
	os.system("> {} ".format(outputfile))
	
	alloutputs = []
	
		#with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
			#futures = []
	for region in regions:
		
			# Submit jobs to executor
		process_locus2( region, haplopath, folder, outputfile)
			#futures.append(executor.submit(process_locus2, region, haplopath, folder, outputfile))
		#alloutputs.append(os.path.join(folder, f"{region[0]}.fa"))
		# Wait for all threads to complete
		#concurrent.futures.wait(futures)
		
def makegraph(fastafile, folder, threads, reflist = []):
	
	#os.system(f"{script_folder}/kmernorm -i {fastafile} -o {fastafile}_norm.gz -w 0 -m 1")
	#os.system(f"rm {fastafile}_part*.fa || true");
	#os.system(f"python {script_folder}/normpartition.py -i {fastafile} -n {fastafile}_norm.gz -o {fastafile}_part")
	
	#files = glob.glob(f"{fastafile}_part*.fa")
	
	files = []
	if not len(files):
		files = [fastafile]
		
	for file in files:
		pass
		os.system("python {}/graphmake.py -i {}  -d {} -t {} -l 0 -p \"{}\" -u 1".format(script_folder, file, folder, threads, ",".join(reflist)))
		
	return files

def process_locus(index, line, haplopath, folder, outputfile):
	line = line.strip().split()
	name = line[0]
	line = line[1:]
	contig = line[0]
	if "#" not in contig: 
		haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
	else:
		haplo = "_h".join(contig.split("#")[:2])
		
	path = haplopath[haplo]
	strd = "-i" if line[1] == "-" else ""
	
	outputfile0 = os.path.join(folder, f"{name}.fa")
	
	line[2] = max(0,int(line[2]))
	header = ">{}\t{}:{}-{}{}\t{}".format(name, contig, line[2], line[3], line[1], line[4])
	
	if os.path.isfile(path):
		cmd1 = f'echo "{header}" > {outputfile0} && samtools faidx {path} {contig}:{line[2]+1}-{line[3]} {strd} | tail -n +2 >> {outputfile0}'
		os.system(cmd1)
		
	return outputfile0


def makenewfasta(locifile, queryfile, outputfile, folder, threads, refprefix):
	
	haplopath = dict()
	with open(queryfile, mode = 'r') as f:
		for line in f:
			line = line.strip().split()
			haplopath[line[0]] = line[1]
			
	os.system("echo > "+outputfile)
	
	refs = []
	alloutputs = []
	with open(locifile, mode='r') as f:
		with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
			futures = []
			for index, line in enumerate(f, start=1):
				name = line.split()[0]
				# Submit jobs to executor
				#process_locus(index, line, haplopath, folder, outputfile)
				futures.append(executor.submit(process_locus, name, line, haplopath, folder, outputfile))
				alloutputs.append(os.path.join(folder, f"{name}.fa"))
			# Wait for all threads to complete
			concurrent.futures.wait(futures)
			
	for result in alloutputs:
		outputfile0 = result
		cmd2 = f"cat {outputfile0} >> {outputfile}"
		os.system(cmd2)
		
		
def _merged_intervals(values, inclusive=False):
	"""Merge flat/nested intervals, optionally converting inclusive ends."""
	if not values:
		return []
	if isinstance(values[0], (int, float, str)):
		if len(values) % 2:
			raise ValueError("intervals contain an odd coordinate count")
		values = [values[i:i + 2] for i in range(0, len(values), 2)]
	intervals = sorted(
		sorted((int(value[0]), int(value[1])))
		for value in values if len(value) >= 2
	)
	merged = []
	for start, end in intervals:
		if inclusive:
			end += 1
		if end <= start:
			continue
		if merged and start <= merged[-1][1]:
			merged[-1][1] = max(merged[-1][1], end)
		else:
			merged.append([start, end])
	return merged


def Overlap2(query, ref):
	"""Intersect interval unions, counting each overlapping base once."""
	query, ref = _merged_intervals(query), _merged_intervals(ref)
	output = []
	i = j = 0
	while i < len(query) and j < len(ref):
		start = max(query[i][0], ref[j][0])
		end = min(query[i][1], ref[j][1])
		if end > start:
			output.append([start, end])
		if query[i][1] <= ref[j][1]:
			i += 1
		else:
			j += 1
	return output

def overlapvalid(regions, graphinfo, validregions, coordinate_system=""):
	
	seq_index = 1
	newregions = cl.defaultdict(list)
	
	def interval_len(a, b):
		if a is None or b is None:
			return 0
		if b < a:
			a, b = b, a
		return max(0, b - a)
	
	for locus, regs in regions.items():
		for region in regs:
			localstart, localend = region[3], region[4]
			
			graphcoordi = _query_interval_boundaries(graphinfo[locus], localstart, localend)
			
			graphspans = cl.defaultdict(list)
			for i in range(len(graphcoordi) // 2):
				path0, start0 = graphcoordi[2 * i][:2]
				path1, end1   = graphcoordi[2 * i + 1][:2]
				if path0 != path1:
					# If this happens, we can still store both, but usually it shouldn’t.
					pass
				graphspans[path0].append(sorted([start0, end1]))
				
			totalvalidoverlap = 0
			for path, pathregions in graphspans.items():
				if path not in validregions:
					continue
				overlaps = Overlap2(pathregions, normalize_validregions(validregions[path], coordinate_system))
				totalvalidoverlap += sum(interval_len(a, b) for a, b in overlaps)
				
			if totalvalidoverlap > min(1000, 0.5 * (localend - localstart)):
				newregions[locus].append([f"seq_{seq_index}", locus, localstart, localend, "", ""])
				seq_index += 1
				
	return newregions

def readscafs(assemfile):
	assemsizes = cl.defaultdict(lambda: 10000000000)
	if len(assemfile):
		with open(assemfile, mode = 'r') as f:
			for line in f:
				path = line.split()[-1]
				with open(path+".fai", mode = 'r') as f:
					for line in f:
						line = line.split()
						assemsizes[line[0]] = int(line[1])
						
	return assemsizes



def findbreaks(inputfile, assemsizes, outputfile, refprefix):
	
	ifexons = set()
	names = []
	contigs = cl.defaultdict(list)
	name_tocontig = cl.defaultdict(list)
	index = -1
	with open(inputfile, mode = 'r') as f:
		
		for line in f:
			
			if line.startswith(">"):
				
				index += 1
				name, locus = line[1:].split()[:2]
				names.append(name)
				exon = line.split()[-2]
				if exon == "Exon":
					ifexons.add(name)
					
				strd = "+"
				contig, region = locus.split(":")
				if region[-1] in ["+","-"]:
					region = region[:-1]
					strd = locus[-1]
					
				region = region.split("-")
				region = [int(region[0]), int(region[1]), name] if strd == "+" else [int(region[1]), int(region[0]), name]
				
				contigs[contig].append(region)
				
	allsize =0 
	new_contigs = cl.defaultdict(list)
	for contig, regions in contigs.items():
		
		new_regions = []
		
		strd = sorted(regions, key = lambda x:  -abs(x[0]-x[1]))[0]
		
		strd = '+' if strd[1] > strd[0] else '-'
		
		regions_sort = sorted(regions) if strd == '+' else sorted(regions, reverse = 1)
		
		lastend = regions_sort[0]
		lastregion = regions_sort[0]
		new_regions = [[regions_sort[0]]]
		
		allsize += abs(regions_sort[0][1]-regions_sort[0][0])
		for region in regions_sort[1:]:
			if max(region[0],region[1], lastend[0],lastend[1]) - min(region[0],region[1], lastend[0],lastend[1]) - abs(region[1]-region[0] ) - abs(lastend[1] - lastend[0]) < 2*cutoffdistance:
				
				new_regions[-1].append(region)
			else:
				new_regions.append([region])
				
			allsize += abs(region[0]-region[1])
			lastend = region
			
		new_contigs.update({(contig+strd,i):regions for i,regions in enumerate(new_regions)})
		
		
	haplocount = cl.defaultdict(int)
	reflist = []
	with open(outputfile, mode = 'w') as w:
		lindex = 1
		for (contig, index), regions in new_contigs.items():
			if "#" not in contig:
				haplo = "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
			else:
				haplo = "_h".join(contig.split("#")[:2])
				
			haplocount[haplo] += 1
			lname = f"{haplo}_{haplocount[haplo]}"
			
			locus = "{}".format(contig[:-1])
			start = min(regions[0][0],regions[0][1],regions[-1][0],regions[-1][1])
			end = max(regions[0][0],regions[0][1],regions[-1][0],regions[-1][1])
			names = ";".join([region[-1] for region in regions])
			
			start = max(0,start-cutoffdistance)
			end = min(end+cutoffdistance, assemsizes[locus])
			
			if locus.startswith(refprefix):
				reflist.append(lname)
				
			w.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(lname, locus,contig[-1],start,end,names))
			
			lindex += 1
			
	return reflist




def polishranges(ranges):
	
	if len(ranges) <= 1:
		return [ranges]
	
	size = abs(ranges[0][2] - ranges[0][1])
	lastend = ranges[0][2]
	breaks = set()
	for i,(index,start,end) in enumerate(ranges[1:]):
		gap = abs(start-lastend)
		if gap > size:
			breaks.add(i+1)
			size = end - start
		else:
			size += end - start - gap
		lastend = end
		
	size = abs(ranges[-1][2] - ranges[-1][1])
	lastend = ranges[-1][1]
	
	for i,(index, start,end) in enumerate(ranges[::-1][1:]):
		gap = abs(end-lastend)
		if gap > size:
			breaks.add(len(ranges) - 2 - i)
			size = end - start
		else:
			size += end - start - gap
		lastend = start
		
	eachbreak = []
	allbreaks = []
	for i,x in enumerate(ranges):
		if i in breaks:
			allbreaks.append(eachbreak)
			eachbreak = []
		eachbreak.append(x)
	allbreaks.append(eachbreak)
	
	return allbreaks

def Overlap(allaligns, allow_gap = 300):
	
	all_coordinates = [x for y in allaligns for x in y[:2]]
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	allgroups = []
	current_group = set([])
	current_group_size = 0 
	current_group_start = 0
	
	number_curr_seq = 0
	coordinate = 0
	current_group = [0,0]
	last_coordinate = -allow_gap-1
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		
		if index %2 == 0 :
			
			number_curr_seq += 1
			
			if number_curr_seq == 1:
				
				current_group_start = coordinate
				
				if coordinate - last_coordinate> allow_gap :
					
					allgroups.append(current_group)
					
					current_group = [coordinate,coordinate] 
		else:
			number_curr_seq -= 1
			current_group[1] = coordinate
			
		last_coordinate = coordinate
		
	allgroups.append(current_group)
	
	allgroups = [x for x in allgroups if x[1] - x[0] > 100]
	
	return allgroups

def qposi_togposi(cigar_str, find_qposes):
	
	l = len(find_qposes)
	if l == 0:
		
		return []
	
	allcigars = re.findall(r'[><][^><]+', cigar_str)
	
	curr_rposi = 0
	curr_qposi = 0
	new_rposi = 0
	new_qposi = 0
	
	curr_strd = 1
	curr_strd_str = ">"
	curr_size = 0
	
	find_qposes = find_qposes + [0]
	find_qpos_index = 0
	find_qpos = find_qposes[0]
	find_rposes = []
	
	for cigars in allcigars:
		
		path, cigar  = cigars.split(":")
		
		if len(find_rposes):
			
			find_rposes.append((pathname,curr_rposi,curr_qposi-1,curr_strd_str ))
			
		pathname = "_".join(path.split("_")[:-2])[1:]
		
		curr_strd_str = cigars [0]
		curr_strd = 1 if curr_strd_str  =='>' else -1
		
		curr_rposi,curr_qposi = path.split("_")[-2:]
		
		curr_rposi,curr_qposi = int(curr_rposi), int(curr_qposi)
		
		if len(find_rposes):
			find_rposes.append((pathname,curr_rposi,curr_qposi,curr_strd_str ))
			
		cigars = re.findall(r'\d+[M=XIHD]', cigars)
		
		if cigars[0][-1] == 'H':
			cigars = cigars[1:]
		if cigars[-1][-1] == 'H':
			cigars = cigars[:-1]
			
		new_rposi, new_qposi = curr_rposi,curr_qposi
		
		for cigar in cigars:
			
			curr_size = int(cigar[:-1])
			char = cigar[-1]
			
			if char == '=' or char ==  'M':
				
				new_rposi = curr_rposi + curr_strd * curr_size
				new_qposi = curr_qposi + curr_size
				
			elif char == 'I':
				
				new_qposi = curr_qposi + curr_size 
				
			elif char == 'D' or char == 'H':
				
				new_rposi = curr_rposi + curr_strd * curr_size
				
			elif char == 'X':
				
				new_rposi = curr_rposi + curr_strd * curr_size
				new_qposi = curr_qposi + curr_size
				
			elif char == 'S':
				
				if curr_qposi == 0:
					curr_qposi = curr_size
					
				curr_size = 0
				continue
			
			while find_qpos_index<l and find_qpos <= new_qposi:
				
				if char == '=' or char ==  'M' or char =='X' :
					
					find_rposes.append((pathname, curr_rposi + curr_strd * (find_qpos-curr_qposi), find_qpos, curr_strd_str ))
				else:
					
					find_rposes.append((pathname,curr_rposi, find_qpos, curr_strd_str ))
					
				find_qpos_index += 1
				find_qpos = find_qposes[find_qpos_index]
				
			if find_qpos_index == l:
				
				return find_rposes
			
			curr_size = 0
			curr_rposi = new_rposi
			curr_qposi = new_qposi
			
			
	while find_qpos_index<l:
		
		find_qpos = find_qposes[find_qpos_index]
		find_rposes.append((pathname, curr_rposi, find_qpos, curr_strd_str ))
		find_qpos_index += 1
		
	return find_rposes


def _query_interval_boundaries(cigar_str, start, end):
	"""Project a half-open query interval to paired graph/query boundaries.

	The public point projector retains its legacy return convention: internal
	component ends carry the last query base (q - 1), whereas graph coordinates
	already mark the boundary. Convert only those synthetic query endpoints;
	requested start/end points are already boundaries. Empty components at an
	exact path junction contribute no interval.
	"""
	if end <= start:
		return []
	coordinates = qposi_togposi(cigar_str, [start, end])
	output = []
	for index in range(0, len(coordinates) - 1, 2):
		left, right = coordinates[index:index + 2]
		if index + 2 < len(coordinates):
			right = (right[0], right[1], right[2] + 1, right[3])
		if right[2] > left[2]:
			output.extend((left, right))
	return output


def gposi_toqposi(cigar_str, find_rposes):
	
	l = len(find_rposes)
	if l == 0:
		
		return []
	
	find_rposes_dict = cl.defaultdict(list)
	
	cigars = re.findall(r'\d+[M=XHDI]', cigar_str)
	
	curr_qposi = 0
	new_qposi = 0
	curr_rposi = 0
	new_rposi = 0
	
	find_rposes = find_rposes+[-1]
	
	find_rpos_index = 0
	find_rpos = find_rposes[find_rpos_index] 
	curr_strd = 1
	find_qposes = []
	for cigar in cigars:
		
		curr_size = int(cigar[:-1])
		char = cigar[-1]
		
		if char == '=' or char ==  'M':
			
			new_rposi = curr_rposi + curr_size 
			new_qposi = curr_qposi + curr_size * curr_strd
		elif char == 'I':
			
			new_qposi = curr_qposi + curr_size * curr_strd
			
		elif char == 'D' or char == 'H':
			
			new_rposi = curr_rposi + curr_size 
			
		elif char == 'X':
			
			new_rposi = curr_rposi + curr_size 
			new_qposi = curr_qposi + curr_size * curr_strd
			
			
		elif char == 'S':
			
			if curr_qposi == 0:
				curr_qposi = curr_size 
				
			curr_size = 0
			continue
		
		while find_rpos_index<l and find_rpos <= new_rposi + (1-curr_strd)-1//2 :
			
			
			if char == '=' or char ==  'M' or char =='X' :
				
				find_qposes.append(curr_qposi + curr_strd * (find_rpos-curr_rposi))
			else:
				
				find_qposes.append(curr_qposi)
				
			find_rpos_index += 1
			find_rpos = find_rposes[find_rpos_index]
			
		if find_rpos_index == l:
			return find_qposes
		
		curr_size = 0
		curr_rposi = new_rposi
		curr_qposi = new_qposi
		
		
	while find_rpos_index<l:
		
		find_rpos_index += 1
		find_qposes.append(curr_qposi)
		
		
	return find_qposes


def getspan_ongraph(name, path, cigar, rposis, qposis):
	
	cigar = cigar.translate(str.maketrans('', '', 'ATCGatcgNn-'))
	
	cigars = re.findall(r'[><][^><]+', cigar)
	pathes = re.findall(r'[><][^<>]+',path)
	
	
	rposis = [ [int(x.split("_")[0]),int(x.split("_")[1])] for x in rposis.split(";")]
	
	qposis = [ [int(x.split("_")[0]),int(x.split("_")[1])] for x in qposis.split(";")]
	
	newcigars= "".join(["{}_{}_{}:{}".format(thename,rposi[0],qposi[0],cigar.split(":")[1]) if thename[0] == '>' else "{}_{}_{}:{}".format(thename,rposi[1],qposi[0],cigar.split(":")[1]) for thename,cigar,rposi,qposi in zip(pathes,cigars,rposis, qposis)])
	
	span = list(zip(pathes, [tuple([int(a) for a in x]) for x in rposis], [tuple([int(a) for a in x]) for x in qposis]))
	
	span = [x if x[0][0] =='>' else tuple([x[0],x[1],sorted(x[2],reverse = 1)])  for x in span]
	
	
	return newcigars, span


def readgraph(graphfile):
	
	contigspan = dict()
	pathsize = dict()
	graphinfo = dict()
	with open(graphfile, mode = 'r') as f:
		
		for line in f:
			
			line = line.strip().split()
			
			if len(line) < 5:
				continue
			
			name, path, cigar, rposis, qposis = line 
			graphinfo[name],contigspan[name]  = getspan_ongraph(name, path, cigar, rposis, qposis)
			
			
	return graphinfo, contigspan



def readcontigsinfo(inputfile, breakfile, assemsizes, refprefix = "NC_0609"):
	
	genetoscaff = dict()
	with open(inputfile, mode = 'r') as f:
		for line in f:
			if line.startswith(">"):
				line = line.strip().split() + ["",""]
				name = line[0][1:]
				locus = line[1]
				strd = locus[-1]
				contig, coordi = locus.split(":")
				if coordi[-1] in ['-','+']:
					coordi = coordi[:-1]
				start, end = map(int, coordi.split("-"))
				start, end = min(start,end), max(start,end)
				start, end = max(0, start),  end
				genetoscaff[name] = [contig, strd, start, end, line[2], line[3]]
				
	blacklist = set()
	locuslocations = dict()
	genetolocus = cl.defaultdict(list)
	with open(breakfile, mode = 'r') as f:
		
		scarf_index = 0
		lastend = -30000
		laststart = -30000
		lastcontig = ""
		lastgenes = set()
		index = 0
		for lindex, line in enumerate(f):
			if line.startswith(">") == False:
				if (line.count('N')+line.count('n')>100) and locusname.startswith("Ref_") == False:
					blacklist.add(locusname)
				continue
			index += 1
			
			line = line.strip().split()
			
			#locusname, contig, strd, start, end = line[0], line[1], line[2],int(line[3]), int(line[4])
			locusname, coordi, names = line[0][1:], line[1], line[-1]
			contig, (start, end), strd = coordi.split(":")[0], map(int,coordi[:-1].split(":")[1].split("-")),  coordi[-1]
			if contig != lastcontig or max(lastend,laststart,end,start) - min(lastend,laststart,end,start) > lastend - laststart + end - start:
				scarf_index += 1
				
			lastcontig = contig
			lastend = end
			laststart = start
			names = line[-1].split(";")
			
			#names = [x for x, xcoordi in genetoscaff.items() if xcoordi[0] == contig and max(xcoordi[3] ,xcoordi[2] ,end,start) - min(xcoordi[3] ,xcoordi[2],end,start) < xcoordi[3] - xcoordi[2] + end - start ]
			
			if len(names) == 1:
				pass 
				
			if refprefix in contig:
				locusname = "Ref_"+locusname
			elif start <= 1 or end >= assemsizes[contig]-1:
				blacklist.add(locusname)
				
			locuslocations[locusname ] = [contig+"_"+str(index), strd, start, end]
			
			for name in names:
				if strd == '+':
					qstart = genetoscaff[name][2] - start
					qend = genetoscaff[name][3] - start
				else:
					qstart = end - genetoscaff[name][3]
					qend = end  - genetoscaff[name][2]
					qstart,qend = min(qstart,qend), max(qstart,qend)
					
				genetolocus[name].append([locusname, qstart, qend] )
				
				genetoscaff[name][0] = contig+"_"+str(index)
				
			genetoscaff[locusname] = [contig+"_"+str(index), strd, start, end]
			
			
	return genetoscaff, genetolocus, locuslocations, blacklist

def combinebreaks_eachpath(names, breaks):
	
	all_coordinates = breaks
	
	sort_index = sorted(range(len(all_coordinates)), key = lambda x: all_coordinates[x])
	
	break_groups = []
	break_group = []
	
	break_group_index = set([])
	break_group_indexs = []
	
	last_coordinate = 0
	for index in sort_index:
		
		coordinate = all_coordinates[index]
		
		if coordinate - last_coordinate < 100:
			
			break_group.append(coordinate)
			break_group_index.add(names[index//2])
			
		else:
			
			break_groups.append(break_group)
			break_group_indexs.append(break_group_index)
			
			break_group = [coordinate]
			break_group_index = set([names[index//2]])
			
		last_coordinate = coordinate
		
	break_groups.append(break_group)
	break_group_indexs.append(break_group_index)
	
	extendgroup = cl.defaultdict(set)
	
	break_count = cl.defaultdict(int)
	break_uniformed = dict()
	break_tocontignames = cl.defaultdict(set)
	for group,names in zip(break_groups,break_group_indexs):
		
		if len(group) == 0:
			continue
		themedian = int(median(group))
		break_tocontignames[themedian] = names
		
		for coordi in group:
			break_uniformed[coordi] = themedian
			
	return break_uniformed, break_tocontignames


def combinebreaks(gbreaks):
	
	coordi_onpath = cl.defaultdict(list)
	name_onpath = cl.defaultdict(list)
	for name, breaks in gbreaks.items():
		
		breaks = [ x for y in breaks for x in y]
		
		for (path, gcoordi, pcoordi, strd,  genename) in breaks:
			
			coordi_onpath[path].append(gcoordi)
			name_onpath[path].append(name)
			
	break_uniformed_new = dict()
	for path, breaks in coordi_onpath.items():
		
		break_uniformed, break_tocontignames = combinebreaks_eachpath(name_onpath[path], breaks)
		
		for posi0, posi1 in break_uniformed.items():
			
			break_uniformed_new[(path,posi0)] = posi1
			
	return break_uniformed_new



def findvalidbreaks(break_contignames, contigspans):
	
	allcontigs_spans = [z  for y in contigspans.values() for x in y for z in [x[0]-10,x[1]+10]]
	
	allnames =  [z  for name, y in contigspans.items() for x in y for z in [name,name]] 
	
	l1 =  len(allcontigs_spans)
	
	allcoordis = allcontigs_spans+ list(break_contignames.keys())
	
	allcoordis_sort = sorted(list(range(len(allcoordis))), key = lambda x: allcoordis[x])
	
	locus_span = 0
	segment_span = 0
	
	locus_lastposi = cl.defaultdict(int)
	
	allbreaks_found = []
	current_names = set()
	lastbreak = 0
	
	breaks_allspancontigs = cl.defaultdict(list)
	for sortindex, index in enumerate(allcoordis_sort):
		
		if index < l1:
			if index % 2 == 0:
				segment_span += 1
				current_names.add(allnames[index])
			else:
				segment_span -= 1
				try:
					current_names.remove(allnames[index])
				except:
					pass
					
		else:
			coordi = allcoordis[index]
			common = break_contignames[coordi] & set(current_names)
			uncommon = set(current_names) ^ common
			breaks_allspancontigs[coordi].extend(current_names)
			
			totalbreak = len(common)
			totalcontigs = len(current_names)
			
			ifref1,ifref2 = any(x for x in common if x.startswith('Ref_')), any(x for x in uncommon if x.startswith('Ref_'))
			ifref = 10
			if ifref1 and not ifref2 :
				ifref = 1
				
			elif ifref2 and not ifref1:
				ifref = -1
			else:
				ifref = 0
				
				
			allbreaks_found.append([coordi,totalbreak, totalcontigs, ifref, uncommon , common]) 
			
			
			
	locus_addbreaks = cl.defaultdict(set)
	locus_removebreaks = cl.defaultdict(set)
	
	refbreaks = dict()
	novelbreaks = set()
	
	for abreak in allbreaks_found:
		
		coordi, totalbreak, totalcontigs, ifref, uncommon, common = abreak
		
		#mindis = min([abs(x-coordi) for x in existingbreak]+[10000000])
		ifuse = 0 
		
		if ifref == 1:
			refbreaks[coordi] = 1
			
		elif ifref == -1:
			refbreaks[coordi] = -1
			
		elif ifref == 0:
			refbreaks[coordi] = 0
			
			
		elif totalbreak > 0.5*totalcontigs:
			novelbreaks.add((totalbreak, coordi))
			
	return refbreaks, novelbreaks, breaks_allspancontigs


class graphDB:
	
	def __init__(self):
		self.breaks = cl.defaultdict(list)
		self.blocks = cl.defaultdict(list)
		self.blocksize = cl.defaultdict(list)
		self.blockset = cl.defaultdict(set)
		self.allgenes = cl.defaultdict(list)
		self.break_uniformed = cl.defaultdict(list)
		self.chunkindex_togene = cl.defaultdict(set)
		self.genechunks_index = cl.defaultdict(list)  
		self.validregions = cl.defaultdict(list)      
		self.cache_metadata = {}
		self.largechunks = set()
		self.smallchunks = set()
		
	def save_json(self, path):
		data = {
		"breaks": self.breaks,
		"blockset": {
			f"{k[0]}|{k[1]}|{k[2]}": v for k, v in self.blockset.items()
		},
		"blocksize": self.blocksize,
		"chunkindex_togene": {
			str(k): [[list(t[0]), t[1]] for t in v]
			for k, v in self.chunkindex_togene.items()
		},
		"genechunks_index": {
			",".join(map(str, k)): v for k, v in self.genechunks_index.items()
		},
		"validregions": self.validregions,
		}
		if self.cache_metadata:
			data["_minsetref"] = self.cache_metadata
		global mergesmall
		if mergesmall != 20000:
			data["mergesmall"] = str(mergesmall)
			
		with open(path, "w") as f:
			json.dump(data, f)
			
			
	@classmethod
	def load_json(cls, path):
		with open(path) as f:
			data = json.load(f)
			
		obj = cls()
		obj.breaks = cl.defaultdict(list, data["breaks"])
		
		obj.blockset = cl.defaultdict(set, {
		(parts[0], int(parts[1]), int(parts[2])): v
		for k, v in data["blockset"].items()
		for parts in [k.rsplit("|", 2)]
		})
		
		obj.blocksize = cl.defaultdict(list, {
		int(k): v for k, v in data["blocksize"].items()
		})
		obj.chunkindex_togene = cl.defaultdict(set, {
		int(k): {(tuple(x[0]), x[1]) for x in v}
		for k, v in data["chunkindex_togene"].items()
		})
		obj.genechunks_index = cl.defaultdict(list, {
		tuple(map(int, k.split(","))) if len(k) else k: v for k, v in data["genechunks_index"].items()
		})
		obj.validregions = cl.defaultdict(list, data.get("validregions", {}))
		obj.cache_metadata = data.get("_minsetref", {})
		
		# Loading a default cache must not inherit a previous cache's cutoff.
		global mergesmall
		mergesmall = int(data.get("mergesmall", 20000))
			
		return obj
	
	def load_breaks(self, break_uniformed):
		
		self.break_uniformed = break_uniformed
		
	def load_refbreaks(self, refbreaks):
		
		refbreaks = [(x[0],self.break_uniformed.get( (x[0],x[1]), x[1] ) ) for x in refbreaks]
		
		for abreak in refbreaks:
			
			self.breaks[abreak[0]].append(abreak[1])
			
		for path, breaks in self.breaks.items():
			
			self.breaks[path] = sorted(list(set(breaks)))
			
	def load_refchunks(self):
		
		index_offsite = 0
		for path, breaks in self.breaks.items():
			
			size = graph_path_size(path)
			
			breaks = [0]+breaks+[size+1]
			intervals = [(start, end) for start, end in zip(breaks, breaks[1:]) if end - start > 100]
			blocks = [(i+index_offsite, path, start, end) for i, (start, end) in enumerate(intervals)]
			
			self.blocks[path].extend(blocks)
			
			for i, block in enumerate(blocks):
				self.blockset[tuple([block[1],block[2],block[3]])] = i+index_offsite
				self.blocksize[i+index_offsite] = block[3] -  block[2]
			index_offsite += len(blocks)
			
		return self
	
	def overlap_chunks(self, chunks):
		
		chunks = sorted(chunks, key = lambda x: x[3] + x[4])
		chunks_new = []
		for chunk in chunks:
			
			
			block_index = self.blockset.get(tuple(chunk[:3]),-1)
			
			if block_index == -1:
				
				overlaps = []
				for block,index in self.blockset.items():
					
					if block[0] != chunk[0]:
						continue
					
					overlap = abs(block[2] - block[1]) + abs(chunk[2]-chunk[1]) - max(block[2], block[1], chunk[2], chunk[1]) + min(block[2], block[1], chunk[2], chunk[1])
					
					if overlap > 0.5 * max(abs(block[2] - block[1]),abs(chunk[2]-chunk[1])):
						
						overlaps.append((overlap, index))
						
				overlaps = sorted(overlaps, reverse = 1)
				
				if len(overlaps) and overlaps[0][0] > 100:
					block_index = overlaps[0][1]
					
			chunks_new.append([block_index] + chunk)
			
		return chunks_new
	
	def load_refgenes(self, loci):
		
		def cutspan_bybreaks(spans_flat, breaks):
			
			path_spans = cl.defaultdict(list)
			for (path, gposi, qposi, strd) in spans_flat:
				
				path_spans[path].append((self.break_uniformed.get((path,gposi), gposi),qposi,strd))
				
				
			allchunks = []
			for path, spans_onpath in path_spans.items():
				
				breaks_onpath = breaks[path]
				
				spans_onpath = [sorted([spans_onpath[2*i],spans_onpath[2*i+1]], key = lambda x: x[0]) for i in range(len(spans_onpath)//2)]
				spans_onpath = [(x[0],x[1],y[0][2]) for y in spans_onpath for x in y]
				
				allgcoordis = [x[0] for x in spans_onpath]
				allqcoordis = [(x[1],1 if x[2] == '>' else -1) for x in spans_onpath]
				
				l1 = len(allgcoordis)
				
				allgcoordis =  allgcoordis  + breaks_onpath
				allgcoordis_sortindex = sorted(list(range(len(allgcoordis))), key = lambda x: allgcoordis[x])
				
				current_segs = []
				
				lastqposi = cl.defaultdict(list)
				lastbreak = cl.defaultdict(int)
				for x in allgcoordis_sortindex:
					
					coordi = allgcoordis[x]
					if x < l1:
						if x %2 == 0:
							current_segs.append(x//2)
							lastbreak[x//2] = coordi
							lastqposi[x//2] = allqcoordis[x]
						else:
							allchunks.append([path, lastbreak[x//2], coordi, lastqposi[x//2][0], allqcoordis[x][0]])
							
							current_segs.remove(x//2)
					else:
						for breakindex in current_segs:
							qposi = lastqposi[breakindex][0] + lastqposi[breakindex][1]  * (coordi - lastbreak[breakindex] )
							
							allchunks.append([path, lastbreak[breakindex], coordi, lastqposi[breakindex][0], qposi])
							lastbreak[breakindex] = coordi
							lastqposi[breakindex] = [qposi, lastqposi[breakindex][1]]
							
							
			chunks = sorted([x for x in allchunks if x[2] - x[1] > 300], key = lambda x: x[3])
			return chunks
		
		self.genechunks = cl.defaultdict(list)
		self.genechunks_index = cl.defaultdict(list)
		self.chunkindex_togene = cl.defaultdict(set)
		for genes, spans, cigarstr in loci:
			
			for gene in genes:
				
				genename = gene[0][-1]
				
				gene = [tuple([x[0],x[1],x[2],x[3]]) for x in gene]
				
				chunks = cutspan_bybreaks(gene, self.breaks)
				
				chunks = [x for x in chunks if abs(x[2] - x[1]) > 300]
				
				chunks_new = self.overlap_chunks( chunks)
				
				self.genechunks[genename] = chunks_new
				
				known_chunks = tuple(x[0] for x in chunks_new if x[0] >= 0)
				self.genechunks_index[known_chunks].append(genename)
				self.genechunks_index[known_chunks[::-1]].append(genename)
				
				for i, chunkindex in enumerate(known_chunks):
					self.chunkindex_togene[chunkindex].add((known_chunks, i))
					self.chunkindex_togene[chunkindex].add((known_chunks[::-1], len(known_chunks)-1-i))
					
		self.genes = list(self.genechunks_index.keys())
		self.genenum = len(self.genechunks_index)
		
		
		return self
	
	def getchunks_fromspan(self, spans, cigarstr):
		
		
		allchunks = []
		
		cigars = re.findall(r'[><][^><]+', cigarstr)
		
		
		path_spans = []
		
		for (path, grange, qrange), cigar in zip(spans,cigars):
			path_spans.append([grange, qrange, path[1:], path[0], cigar])
			
	
		pathsizes = dict()
		for grange, qrange, path, strd, cigar in path_spans:
			
			qstart,qend = min(qrange),max(qrange)
			
			breaks_onpath = self.breaks[path]
			pathsize = graph_path_size(path, cigar)
			
			#allgcoordis = [self.break_uniformed.get((path, x),x) for y in spans_onpath for x in y[0]]
			allbreaks = [x for x in breaks_onpath if x >= grange[0] and x < grange[1]]			
			
			allgcoordis = [grange[0]] + allbreaks + [grange[1]]
			if strd == "<":
				breaks_qposi = gposi_toqposi(cigar, [pathsize - x for x in allgcoordis][::-1])
				breaks_qposi = breaks_qposi[::-1]
			else:
				breaks_qposi = gposi_toqposi(cigar, allgcoordis)
				
			breaks_qposi = [x + qstart for x in breaks_qposi]
			
			lastgposi, lastqposi = allgcoordis[0], breaks_qposi[0]
			for gposi, qposi in zip(allgcoordis[1:], breaks_qposi[1:]):
				
				if abs(qposi - lastqposi)>200:
					allchunks.append([path,lastgposi, gposi, lastqposi, qposi])
				lastgposi, lastqposi =  gposi,qposi
				
				
		allchunks = sorted([[x[0]]+sorted(x[1:3])+sorted(x[3:]) for x in allchunks if abs(x[2] - x[1])+abs(x[4]-x[3]) > 300], key = lambda x: x[3]+x[4])
		allchunks_index = self.overlap_chunks(allchunks)
		
		return allchunks_index
	
	def update_database(self, allchunks, usedchunks, results):
		
		index_tochunkindex = dict()
		next_chunk_index = max(self.blockset.values(), default=-1) + 1
		
		newchunks = []
		for i,chunk in enumerate(allchunks):
			if chunk[0] == -1:
				
				findoverlap = 0
				for pastnewchunk in newchunks:
					if pastnewchunk[1] != chunk[1]:
						continue
					overlap =  abs(chunk[3]-chunk[2]) +  abs(pastnewchunk[3]-pastnewchunk[2])  - max(pastnewchunk[3], pastnewchunk[2], chunk[3] , chunk[2]) + min(pastnewchunk[3], pastnewchunk[2], chunk[3] , chunk[2])
					if overlap > 0.5 * max(abs(chunk[3]-chunk[2]),abs(pastnewchunk[3]-pastnewchunk[2])):
						newindex = pastnewchunk[0]
						findoverlap = 1
						break
					
				if findoverlap == 0:
					newindex = next_chunk_index
					next_chunk_index += 1
					self.blocks[chunk[1]].append((newindex, chunk[1], chunk[2], chunk[3]))
					self.blockset[(chunk[1],chunk[2],chunk[3])] = newindex
					self.blocksize[newindex] = abs(chunk[3] - chunk[2])
					index_tochunkindex[i] = newindex
					chunk[0] = newindex
					newchunks.append(chunk)
				else:
					index_tochunkindex[i] = newindex
					chunk[0] = newindex
					
		novels = []
		"""
		novel = []
		lastqend = -1000000
		for i,x in enumerate(allchunks):
			if i == -1:
				qstart,qend = x[4],x[5]
				if abs(qstart - lastqend) < 30000:
					novel.append([i,qstart,qend])
				else:
					if len(novel):
						novels.append(novel)
					novel = [[i,qstart,qend]]
				lastqend = qend
				
		if len(novel):
			novels.append(novel)
		novels = [x for y in novels for x in polishranges(y) if len(x)]
		"""
		
		for result in results:
			chunks_new = tuple([index_tochunkindex.get(x,y) for x,y in zip(result[0],result[1])])
			if chunks_new in self.genechunks_index or len(chunks_new) <= 1:
				continue
			
			genename = "New_"+str(len(self.genechunks))
			
			self.genechunks[genename] = chunks_new
			self.genechunks_index[chunks_new].append(genename)
			self.genechunks_index[chunks_new[::-1]].append(genename)
			
			for i,x in enumerate(chunks_new):
				self.chunkindex_togene[x].add((chunks_new, i))
				self.chunkindex_togene[x].add((chunks_new[::-1], len(chunks_new)-1-i))
				
		for chunks_new in novels:
			
			chunks_new = [x[0] for x in chunks_new]
			chunks_new = tuple(range(min(chunks_new),max(chunks_new)+1))
			
			if chunks_new in self.genechunks_index or len(chunks_new) <= 1:
				continue
			
			genename = "Novel_"+str(len(self.genechunks))
			
			self.genechunks[genename] = chunks_new
			self.genechunks_index[chunks_new].append(genename)
			self.genechunks_index[chunks_new[::-1]].append(genename)
			
			for i,x in enumerate(chunks_new):
				self.chunkindex_togene[x].add((chunks_new, i))
				self.chunkindex_togene[x].add((chunks_new[::-1], len(chunks_new)-1-i))
				
		self.genes = set(list(self.genechunks_index.keys()))
		self.genenum = len(self.genechunks_index)
		
		
	def overlap_genes(self, spans, cigarstr, ifupdate = 0):
		"""Match query chunks to reference paths with polynomial dynamic programming.

		The legacy implementation cloned a complete candidate for both the
		"skip" and "extend" choices.  Repeated traversals of one graph path
		therefore generated exponentially many histories.  For a fixed first and
		last query chunk, reference path, and path position, only the
		highest-scoring history can affect any future transition or final greedy
		selection.  The dictionary frontier below applies that dominance rule.
		"""

		def keep_best_active(states, state):
			key = state.active_key()
			old = states.get(key)
			if old is None or state.score > old.score:
				states[key] = state

		def determinegenes(allchunks, passed_genes):
			selected_genes = []
			exclude = set()
			for gene in passed_genes.iter_sorted():
				chunk_usedindex = range(gene[0][0], gene[0][-1] + 1)
				if any(x in exclude for x in chunk_usedindex):
					continue
				selected_genes.append(gene)
				exclude.update(chunk_usedindex)

			for ichunk, chunk in enumerate(allchunks):
				chunkindex = chunk[0]
				qstart,qend = chunk[4],chunk[5]
				if chunkindex >= 0 and ichunk not in exclude:
					selected_genes.append([
						[ichunk], [chunkindex], [], qstart, qend,
						[[qstart, qend]],
					])
			return selected_genes, exclude

		# Lazily build O(1)/O(log n) transition data for each reference path.
		# This replaces tuple.index() and a repeated deletion-size summation in
		# every active candidate.
		gene_transition_cache = {}

		def transition_data(gene):
			cached = gene_transition_cache.get(gene)
			if cached is not None:
				return cached
			positions = cl.defaultdict(list)
			prefix_sizes = [0]
			for position, graph_chunk in enumerate(gene):
				positions[graph_chunk].append(position)
				prefix_sizes.append(
					prefix_sizes[-1] + self.blocksize[graph_chunk]
				)
			cached = (positions, prefix_sizes)
			gene_transition_cache[gene] = cached
			return cached

		current_genes = {}
		passed_genes = _BestIntervalStates()
		allchunks = sorted(
			self.getchunks_fromspan(spans, cigarstr),
			key=lambda x: x[4] + x[5],
		)

		for ichunk, chunk in enumerate(allchunks):
			chunkindex = chunk[0]
			qstart,qend = chunk[4],chunk[5]
			allinvolve_genes = self.chunkindex_togene.get(chunkindex, ())
			if not allinvolve_genes and not current_genes:
				continue

			continued_genes = {}
			# Longer matches were processed first by the legacy code.  Preserve
			# that order solely as a deterministic tie breaker.
			active = sorted(
				current_genes.values(), key=lambda state: state.depth, reverse=True,
			)
			for state in active:
				insersize = qstart - state.qend
				if insersize > 30000:
					passed_genes.add(state)
					continue

				positions, prefix_sizes = transition_data(state.gene)
				matching_positions = positions.get(chunkindex, ())
				position_index = bisect.bisect_right(
					matching_positions, state.alignindex,
				)
				matched = False
				if position_index < len(matching_positions):
					newindex = matching_positions[position_index]
					delsize = (
						prefix_sizes[newindex]
						- prefix_sizes[state.alignindex + 1]
					)
					newscore = (
						qend - qstart - 2 * insersize - 2 * delsize
					)
					if delsize < 30000 and newscore >= 0:
						matched = True
						child = _DPOverlapState(
							state, ichunk, chunkindex, state.gene,
							newindex, state.score + newscore, qstart, qend,
						)
						keep_best_active(continued_genes, child)

				if matched:
					if not ifupdate:
						keep_best_active(continued_genes, state)
				elif ifupdate:
					keep_best_active(continued_genes, state)
				else:
					passed_genes.add(state)

			for gene, alignindex in allinvolve_genes:
				root = _DPOverlapState(
					None, ichunk, chunkindex, gene, alignindex,
					qend - qstart, qstart, qend,
				)
				keep_best_active(continued_genes, root)
			current_genes = continued_genes

		passed_genes.extend(current_genes.values())
		results, used_chunks = determinegenes(allchunks, passed_genes)
	
	
		for result in results:
			result[2] = [(x[2],self.genechunks_index.get(x[0], [""])[0]) for x in result[2]]
			if len(result[2]):
				result[2] = result[2][0]
				
		results = sorted(results, key = lambda x: tuple(x[0]))
	
		if ifupdate:
			self.update_database(allchunks, used_chunks, results)
			
			
		return results	
	
	
def getblocks(gbreaks, break_uniformed, contigspans, genetoscaff, graphinfo, blacklists):
	
	saveorload = 0
	allloci_wlist = [x for x in gbreaks.keys() if x.startswith("Ref_")] + [x for x in gbreaks.keys() if not x.startswith("Ref_") if x not in blacklists]  
	allloci_blist = [x for x in gbreaks.keys() if not x.startswith("Ref_") if x in blacklists]
	if saveorload == 0:
		
		thegraph = graphDB()
		thegraph.load_breaks(break_uniformed)
		
		refgenes = []
		for name, gbreak in gbreaks.items():
			
			if name.startswith("Ref_"):
				
				thegraph.load_refbreaks([(x[0],x[1]) for y in gbreak for x in y])
				refgenes.append([gbreaks[name], contigspans[name], graphinfo[name]])
				
		thegraph.load_refchunks().load_refgenes(refgenes)

		counter = 0
		results = dict()
		for loci in allloci_wlist:
			results[loci] = thegraph.overlap_genes(contigspans[loci], graphinfo[loci], ifupdate=1)
		for loci in allloci_wlist:
			results[loci] = thegraph.overlap_genes(contigspans[loci], graphinfo[loci], ifupdate=0)
			
		for loci in allloci_blist:
			results[loci] = thegraph.overlap_genes(contigspans[loci], graphinfo[loci], ifupdate=1)
		for loci in allloci_blist:
			results[loci] = thegraph.overlap_genes(contigspans[loci], graphinfo[loci], ifupdate=0)
		
	"""
	if saveorload == 1:
		thegraph = graphDB()
		thegraph = thegraph.load_json(cache)

		results = dict()
		for loci in allloci:
			results[loci] = thegraph.overlap_genes(contigspans[loci], graphinfo[loci])
	"""
			
	return results,thegraph

def breaks_ongraph(genetoscaff, genetolocus, graphinfo, contigspan, blacklists):
	
	gbreaks = cl.defaultdict(list)
	segmentinfo =  cl.defaultdict(list)
	for name, infos in genetolocus.items():
		
		for info in infos:
			
			locus,start,end = info
			
			if locus.startswith("Ref_"):
				graphinfo[locus] = graphinfo[locus[4:]]
				if locus[4:] in contigspan:
					contigspan[locus] = contigspan[locus[4:]]
					del contigspan[locus[4:]]
					
			if locus not in graphinfo:
				continue
			
			coordinates = _query_interval_boundaries(graphinfo[locus], start, end)
			if not coordinates:
				continue
			
			gbreaks[locus].append([list(x)+[name] for x in coordinates])
			
			segmentinfo[locus].append([coordinates[0],coordinates[1], name])
			
	gbreaks = {locus:sorted(regions, key = lambda x: x[0][2]) for locus, regions in gbreaks.items() }
	
	breaks_andspans = cl.defaultdict(list)
	for name, data in gbreaks.items():
		breaks_andspans[name] = [x for x in data]
	for name, data in contigspan.items():
		breaks_andspans[name].append([[x[0][1:],x[1][y],x[2][y],x[0][0],""] for x in data for y in [0, 1]])
		
	break_uniformed = combinebreaks(breaks_andspans)
	results,thegraph = getblocks(gbreaks, break_uniformed, contigspan, genetoscaff, graphinfo, blacklists)
	
	return results,thegraph

def assemblysmall(allresults, cutoff=20000, only_unrepresented=False):
	"""Merge short neighboring blocks using the legacy Assemblysmall policy.

	When ``only_unrepresented`` is true, only query-only fallback blocks (the
	negative chunk indexes created by ``lossless_validregion_results``) are
	eligible to trigger a merge.  A merged block remains eligible because it
	retains that negative chunk.  This lets sample-specific novel structure be
	absorbed recursively after cached reference merges have been replayed,
	without changing boundaries between established reference blocks.
	"""
	
	allowdistance = 20000
	def is_merge_candidate(result):
		if not only_unrepresented:
			return True
		chunks = result[1] if isinstance(result[1], list) else [result[1]]
		return any(int(chunk) < 0 for chunk in chunks)

	def ifalign(thelist,x):
		
		l = len(x)
		for i in range(len(thelist)):
			if x == thelist[i:(i+l)]:
				return 1
		return 0
	
	smallchunks = set()
	bindrecords = cl.defaultdict(int)
	iffix = set()
	for loci,results in allresults.items():
		
		if len(results) == 1:
			continue
		
		for i,result in enumerate(results):
			if not is_merge_candidate(result):
				continue
			
			size = result[4] - result[3]
			
			if size < cutoff:
				
				distance1 = abs(result[3] - results[i-1][4]) if i > 0 else 2000000
				distance2 = abs(result[4] - results[i+1][3]) if i != len(results)-1 else 2000000
				
				if distance1 > allowdistance and distance2 > allowdistance:
					continue
				elif distance1 > allowdistance:
					iffix.add(loci)
					smallchunks.add(tuple(sorted(result[1])))
					bindindex2 = tuple(sorted([result[1][-1], results[i+1][1][0]]))
					bindrecords[bindindex2] += 1
					continue
				
				elif distance2 > allowdistance:
					iffix.add(loci)
					smallchunks.add(tuple(sorted(result[1])))
					bindindex1 = tuple(sorted([result[1][0], results[i-1][1][-1]]))
					bindrecords[bindindex1] += 1
					continue
				else:
					iffix.add(loci)
					smallchunks.add(tuple(sorted(result[1])))
					bindindex1 = tuple(sorted([result[1][0], results[i-1][1][-1]]))
					bindindex2 = tuple(sorted([result[1][-1], results[i+1][1][0]]))
					
					if abs(results[i+1][4] - results[i][3]) >= abs(results[i-1][3] - results[i][4]):
						bindrecords[bindindex2] += 1
					else:
						bindrecords[bindindex1] += 1
						
						
	for loci,results in allresults.items() :
		
		if loci not in iffix:
			continue
		
		connects = [0 for x in results]
		
		for i,result in enumerate(results):
			
			if (
				is_merge_candidate(result)
				and tuple(sorted(result[1])) in smallchunks
			):
				
				distance1 = abs(result[3] - results[i-1][4]) if i > 0 else 2000000
				distance2 = abs(result[4] - results[i+1][3]) if i != len(results)-1 else 2000000
				
				if distance1 > allowdistance and distance2 > allowdistance:
					continue
				elif distance1 > allowdistance:
					connects[i] = 1
					continue
				elif distance2 > allowdistance:
					connects[i-1] = 1
					continue
				else:
					
					distance1 = abs(result[3] - results[i-1][4])
					distance2 = abs(result[4] - results[i+1][3])
					
						
					if distance1 > allowdistance:
						connects[i] = 1
					elif distance2 > allowdistance:
						connects[i-1] = 1
					else:
						bindindex1 = tuple(sorted([result[1][0], results[i-1][1][-1]]))
						bindindex2 = tuple(sorted([result[1][-1], results[i+1][1][0]]))
						if bindrecords[bindindex2] > bindrecords[bindindex1]:
							connects[i] = 1
						else:
							connects[i-1] = 1
						
		lastindex = 0
		newresults = [results[0]]
		for i,result in enumerate(results[1:]):
			
			if connects[i] == 1:
				
				newresults[-1][0].extend(result[0])
				newresults[-1][1].extend(result[1])
				newresults[-1][3] = min(newresults[-1][3], result[3])
				newresults[-1][4] = max(newresults[-1][4], result[4])
				newresults[-1][5].extend(result[5])
				
			else:
				newresults.append(result)
				
		results = newresults
		
		allresults[loci] = results
		
	for loci,results in allresults.items() :
		
		if len(results) == 0:
			continue
		
		laststart, lastend  = results[0][3], results[0][4]
		for i,result in enumerate(results[1:]):
			
			start, end = result[3], result[4]
			
			distance = max(laststart, lastend,start, end ) - min(laststart, lastend,start, end ) - abs(laststart-  lastend) - abs(start - end)
			boundary_is_eligible = (
				not only_unrepresented
				or is_merge_candidate(results[i])
				or is_merge_candidate(result)
			)
			if boundary_is_eligible and distance > 1 and distance < 1000:
				if laststart < start:
					allresults[loci][i][4] = start
				else:
					allresults[loci][i-1][4] = laststart
					
			laststart, lastend = start, end
			
	return allresults, iffix


def clean_nonoverlap(breaksoncontigs, locuslocations, geneoncontigs):
	genelabs = cl.defaultdict(lambda: [0,[]])
	alloldsegments = cl.defaultdict(list)
	for loci, breaks in breaksoncontigs.items():
		
		scaflocus = locuslocations[loci]
		
		lastend = 0
		for segments in breaks:
			
			for i,segment in enumerate(segments[5]):
				
				start,end = min(segment[0],segment[1]),max(segment[0],segment[1])
				
				segment[0] = max(start, lastend)
				segment[1] = max(end, lastend)
					
				lastend = max(end,lastend)
				
		alloldsegments[loci] = geneoncontigs[scaflocus[0]]
		
		for segment in breaks:
			
			subsegments_indexs, subsegments_ranges = segment[1],  segment[5]
			
			for segmentindex, subsegments_range in zip(subsegments_indexs, subsegments_ranges):
				
				start, end = subsegments_range[0], subsegments_range[1]
				if end <= start:
					continue
				
				if scaflocus[1] == '+':
					scafstart = scaflocus[2] + start
					scafend = scaflocus[2] + end
				else:
					scafstart = scaflocus[3] - end
					scafend = scaflocus[3] - start
					
				overlap = sorted([( (x[1] - x[0]) + (scafend - scafstart) - max(x[1],scafend,x[0],scafstart) + min(x[1],scafend,x[0],scafstart), x) for i,x in enumerate(alloldsegments[loci])], reverse=1)[0]   
				
				if overlap[0] > genelabs[segmentindex][0]:
					genelabs[segmentindex] = overlap
					
					
	regions = []
	for loci, breaks in breaksoncontigs.items():
		
		scaflocus = locuslocations[loci]
		
		newsegments = []
		for i,segment in enumerate(breaks):
			
			subsegments_indexs, subsegments_ranges = segment[1],  segment[5]
			
			start = 100000000000000
			end = -1
			query_indexs_new = []
			subsegments_indexs_new = []
			subsegments_ranges_new = []
			for j, (segmentindex, subsegments_range) in enumerate(zip(subsegments_indexs, subsegments_ranges)):
				
				size = subsegments_range[1] - subsegments_range[0]
				if size <= 0:
					continue
				overlap = genelabs[segmentindex] 
				
				if len(overlap) and overlap[0] > min(1000, 0.5*size):   
					start = min(start, subsegments_range[0])
					end = max(end, subsegments_range[1])
					
					subsegments_indexs_new.append(segmentindex)
					subsegments_ranges_new.append(subsegments_range)
					query_indexs_new.append(segment[0][j])
					
			if end == -1:
				continue
			
			newsegment = segment
			
			newsegment[0] = query_indexs_new
			newsegment[1] =  subsegments_indexs_new
			newsegment[5] =  subsegments_ranges_new
			
			newsegment[3] = start
			newsegment[4] = end
			
			newsegments.append(newsegment)
			
		breaksoncontigs[loci] = newsegments
		
	return breaksoncontigs, alloldsegments, genelabs


def annotate_regions(
	breaksoncontigs, locuslocations, genetoscaff, return_block_results=False,
):
	
	geneoncontigs = cl.defaultdict(list)
	for genename, info in genetoscaff.items():
		
		if len(info) > 4:
			geneoncontigs[info[0]].append(info[2:]+[genename])
			
			groupprefix = genename.split('_')[0]
			
	novel_overlap = [x[2][1]  for loci, breaks in breaksoncontigs.items() for x in breaks if len(x[2]) ]
	novel_overlap = list(set(novel_overlap))
	breaksoncontigs, alloldsegments, genelabs = clean_nonoverlap(breaksoncontigs, locuslocations, geneoncontigs)
	
	breaksoncontigs = {k:v for k,v in breaksoncontigs.items() if len(v)}
	
	breaksoncontigs_ori = copy.deepcopy(breaksoncontigs)

	maxsize = max([max([y[4]-y[3] for y in x]+[0]) for x in breaksoncontigs.values()]+[0])
	global mergesmall
	iffix = 1
	while iffix:
		breaksoncontigs, iffix = assemblysmall(breaksoncontigs, mergesmall)
		
	maxsize = max([max([y[4]-y[3] for y in x]+[0]) for x in breaksoncontigs.values()]+[0])
	if maxsize > 100000:
		mergesmall = mergesmall//2
		breaksoncontigs = breaksoncontigs_ori
		
		maxsize = max([max([y[4]-y[3] for y in x]+[0]) for x in breaksoncontigs_ori.values()]+[0])
		iffix = 1
		while iffix:
			breaksoncontigs, iffix = assemblysmall(breaksoncontigs, mergesmall)
			
	regions = []
	for loci, breaks in breaksoncontigs.items():
		
		scaflocus = locuslocations[loci]
		
		for segment in breaks:
			
			start, end = segment[3], segment[4]
			
			if scaflocus[1] == '+':
				scafstart = scaflocus[2] + start
				scafend = scaflocus[2] + end
			else:
				scafstart = scaflocus[3] - end
				scafend = scaflocus[3] - start
				
			overlap = sorted([( (x[1] - x[0]) + (scafend - scafstart) - max(x[1],scafend,x[0],scafstart) + min(x[1],scafend,x[0],scafstart), x) for i,x in enumerate(alloldsegments[loci])], reverse=1)[0]   
			
			regions.append([scaflocus[0], int(scafstart), int(scafend)] + overlap[1][2:]+[loci,start,end])
			
	new_regions = cl.defaultdict(list)
	haplo_counter = cl.defaultdict(int)
	orginal_loci = dict()
	for region in regions:
		
		contig, start, end, ifexon, mapinfo,oldname,locus,localstart,localend = region
		
		haplo = "_h".join(contig.split("#")[:2]) if "#" in contig else "CHM13_h1" if "NC_0609" in contig else "HG38_h1"
		
		haplo_counter[haplo] += 1
		
		newname = groupprefix+"_"+ haplo +"_"+ str(haplo_counter[haplo])
		
		#header = "{}\t{}:{}-{}\t{}\t{}".format(newname, contig, start, end, ifexon, mapinfo  )
		
		new_regions[contig].append([newname,  contig, start, end, ifexon, mapinfo])
		
		# Region records and assemblysmall results use exclusive ends. Keep
		# this auxiliary map's documented legacy inclusive end so existing
		# build_total_validregions callers continue to convert it exactly once.
		orginal_loci[newname] = (locus,localstart,localend - 1)
		
	new_regions = [x[:1]+["_".join(x[1].split("_")[:-1])]+x[2:] for y in new_regions.values() for x in y]
	
	if return_block_results:
		return new_regions, orginal_loci, breaksoncontigs
	return new_regions, orginal_loci


def coordinate_uniform(cbreaks_sort):
	
	uniform_coordis = []
	
	currbreaks =  []
	lastbreak = -100
	for coordi in cbreaks_sort:
		if (coordi - lastbreak) < 100:
			currbreaks.append(coordi)
		else:
			if len(currbreaks):
				uniform = cl.Counter(currbreaks).most_common(1)[0][0]
				uniform_coordis.extend([uniform] * 1)
				
			currbreaks = [coordi]
			
		lastbreak = coordi
		
	if len(currbreaks):
		uniform = cl.Counter(currbreaks).most_common(1)[0][0]  
		uniform_coordis.extend([uniform] * 1)
		
		
	return uniform_coordis



def locatebreaksoncontigs(breaksonscaf, locuslocations, genetoscaff):
	
	oldbreaksoncontig = cl.defaultdict(list)
	for gene,obreaks in genetoscaff.items():
		
		scaf, strd, start, end = obreaks[:4]
		
		oldbreaksoncontig[scaf].extend([start, end])
		
	results = cl.defaultdict(list)
	for scafford, cbreaks in breaksonscaf.items():
		
		if len(cbreaks) == 0:
			continue
		
		cbreaks_sort_index = sorted(range(len(cbreaks)), key = lambda x: cbreaks[x])
		
		cbreaks_sort = sorted([cbreaks[x] for x in cbreaks_sort_index])
		
		cbreaks_sort_uniform = coordinate_uniform(cbreaks_sort)
		
		oldmin = min(oldbreaksoncontig[scafford]+[-100]) 
		oldmax = max(oldbreaksoncontig[scafford]+[100000000000])
		
		cbreaks_sort_uniform_rangeindex = [i for i,x in enumerate(cbreaks_sort_uniform) if x > oldmin + 100 and x< oldmax-100 ]
		
		if len(cbreaks_sort_uniform_rangeindex) == 0:
			
			if min(cbreaks_sort_uniform_rangeindex) > 0:
				cbreaks_sort_uniform_rangeindex = [min(cbreaks_sort_uniform_rangeindex)-1] + cbreaks_sort_uniform_rangeindex
			if max(cbreaks_sort_uniform_rangeindex) < len(cbreaks_sort_uniform_rangeindex) - 1:
				cbreaks_sort_uniform_rangeindex =  cbreaks_sort_uniform_rangeindex + [max(cbreaks_sort_uniform_rangeindex)+1]
				
		cbreaks_sort_uniform = [cbreaks_sort_uniform[x] for x in cbreaks_sort_uniform_rangeindex]
		
		if len(cbreaks_sort_uniform ) == 0:
			cbreaks_sort_uniform = oldbreaksoncontig[scafford]
		elif  len(cbreaks_sort_uniform ) <= 1 or oldmin < min(cbreaks_sort_uniform) - 1000 or oldmax > max(cbreaks_sort_uniform) + 1000:
			cbreaks_sort_uniform = [x for x in oldbreaksoncontig[scafford] if x < min(cbreaks_sort_uniform) -1000 ] + cbreaks_sort_uniform + [x for x in oldbreaksoncontig[scafford] if x > max(cbreaks_sort_uniform) + 1000 ]
			
			#cbreaks_sort_uniform = sorted(oldbreaksoncontig[scafford])
			
		results[scafford] = cbreaks_sort_uniform 
		
		
	return results

def MergeChunks(regions, mergedistance):
	
	regions_sort = sorted([x for x in regions], key = lambda x: (x[1], x[2]))
	
	lastcontig = ""
	lastend = -mergedistance - 1
	
	merged_regions = []
	current_region = []
	for region in regions_sort:
		
		newname,  contig, start, end, ifexon, mapinfo = region
		if contig != lastcontig:
			lastend = -mergedistance - 1
			
		if start - lastend < mergedistance:
			current_region[3] = max(current_region[3], end)
			current_region[4] = "Exon" if ifexon == "Exon" else current_region[4]
			
		else:
			merged_regions.append(current_region)
			current_region = [x for x in region]
			
		lastend = current_region[3]
		lastcontig = contig
		
		
	merged_regions = merged_regions[1:]
	merged_regions.append(current_region)
	
	
	if len([x for x in merged_regions if len(x)]) == 0:
		return 0, []
	
	maxsize = max([x[3]-x[2] for x in merged_regions])
	if maxsize < mergedistance:
		
		regions = merged_regions
		
	return maxsize, regions

def uniformbreaks(inputfile, haplomergefile, graphfile, assemsizes, refcontig):
	
	genetoscaff, genetolocus, locuslocations, blacklists = readcontigsinfo(inputfile, haplomergefile, assemsizes, refcontig)
	graphinfo, contigspan = readgraph(graphfile)
	results, thegraph = breaks_ongraph(genetoscaff, genetolocus, graphinfo, contigspan, blacklists)
	regions, orginal_loci, final_block_results = annotate_regions(
		results, locuslocations, genetoscaff, return_block_results=True,
	)
	
	#normalized = json.dumps(sorted(regions, key=lambda x: tuple(map(str, x))), sort_keys=True)
	#hashvalue = hashlib.sha256(normalized.encode()).hexdigest()
	
	if cache:
		thegraph.validregions = build_total_validregions(
			regions, orginal_loci, graphinfo, contigspan,
		)
		thegraph.cache_metadata["coordinate_system"] = "0-based-half-open"
		thegraph.cache_metadata["validregions_semantics"] = (
			"reference-block-graph-closure-v1"
		)
		thegraph.cache_metadata[REFERENCE_ASSEMBLYSMALL_KEY] = (
			reference_assemblysmall_config(final_block_results, mergesmall)
		)
		
		thegraph.save_json(cache)
		
		
	return regions

def main(args):
	
	global mergesmall 
	mergesmall = args.mergesmall
	queryfile = args.query
	
	if len(args.folder)==0: 
		folder = args.output + "_temp"+ datetime.datetime.now().strftime("%y%m%d_%H%M%S/").replace("%","_")
		os.system("rm -rf {} || true".format(folder))
		os.makedirs(folder,exist_ok=True)
		
	else:
		folder = args.folder
		os.makedirs(folder,exist_ok=True)
		
	scafsizes = readscafs(args.query)
	haplomergefile = args.output+"_loci.txt.fasta"
	ifmakegraph = 1
	if ifmakegraph == 1:
		reflist = findbreaks(args.input, scafsizes, args.output+"_loci.txt", args.ref)
		makenewfasta(args.output+"_loci.txt", args.query, haplomergefile, folder ,args.threads, args.ref)
		os.system(f"{script_folder}/kmerstrd -i {haplomergefile} -o {haplomergefile}_ -r {queryfile}  && mv {haplomergefile}_ {haplomergefile}")
		makegraph(haplomergefile, folder, args.threads, reflist)
	if args.cache:
		global cache
		cache = haplomergefile + "_graphcache.json"
		
	regions = uniformbreaks(args.input, args.output+"_loci.txt.fasta", haplomergefile+"_allgraphalign.out", scafsizes, args.ref)
	maxsize, regions = MergeChunks(regions, args.mergeall)
	if len(regions) == 0:
		print(f"gfixbreak.py Error, Missing regions: {args.input}")
	makenewfasta2(regions, args.query, args.output, folder, args.threads)
	if len(args.folder)==0:
		os.system("rm  -rf {} || true".format(folder))
		
	#os.system(f" cat {haplomergefile} | grep \"^>\" > {haplomergefile}_   &&  mv {haplomergefile}_ {haplomergefile} ")
	os.system(f"rm {haplomergefile}_norm.gz || true ")
	return

def run():
	"""
			Parse arguments and run
	"""
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
	parser.add_argument("-o", "--output", help="path to output file", dest="output",type=str, required=True)
	parser.add_argument("-t", "--threads", help="path to output file", dest="threads",type=int, default = 1)
	parser.add_argument("-q", "--query", help="path to output file", dest="query",type=str, required= True)
	parser.add_argument("-d", "--folder", help="path to output file", dest="folder",type=str, default = "")   
	parser.add_argument("-m", "--mergeall", help="path to output file", dest="mergeall",type=int, default = 80000)
	parser.add_argument("-s", "--mergesmall", help="path to output file", dest="mergesmall",type=int, default = 20000)
	parser.add_argument("-r", "--ref", help="prefix of ref contig", dest="ref",type=str, default = "NC_0609")
	parser.add_argument("--cache", help="save caches", action="store_true")
	parser.set_defaults(func=main)
	args = parser.parse_args()
	args.func(args)
	
	
if __name__ == "__main__":
	run()
