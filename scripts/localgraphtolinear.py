#!/usr/bin/env python3
"""Linearize one local graph using an MSA of all block-index paths.

The legacy ``graphtolinear.py`` graph/CIGAR projection code is retained, but
its heuristic graph walk is replaced at runtime by the deterministic
progressive token-MSA implemented here.  Tokens are graph block indexes (for
example ``12`` and the second repeated occurrence ``12_2``), not nucleotide
bases.

The default (polish=0) projection retains optional query-base payloads on I/X,
including clipped operations and bases moved into gaps between segments.
M/= denote exact matches; mismatches must be explicit X operations. Missing
payloads remain length-only operations and do not prevent linearization.
The query FASTA supplies names, order and lengths, not replacement bases.

Output uses the existing custom S/L format (despite the .gaf suffix), with
the original nine L columns followed by QS:Z:+/- for ten-column inputs.
QS records the original row strand. UC:Z preserves an unaligned row's CIGAR.
tools/reconstruct_local_linear.py can recover query sequences from this file
alone, provided all query bases are represented and every I/X payload is
complete. M/= bases (including their case) come from the graph; original
query soft-masking cannot be recovered independently. Five-column inputs
have no whole-query strand metadata. Polishing is a separate legacy path
and does not have the payload-preservation guarantee.
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from functools import partial
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


GAP = None
LEGACY_MODULE_NAME = "minsetref_legacy_graphtolinear"


def _token_base(token: str) -> str:
    return str(token).split("_", 1)[0]


def _column_score(column: Sequence[Optional[str]], token: str) -> int:
    present = [value for value in column if value is not None]
    if token in present:
        return 4
    if any(_token_base(value) == _token_base(token) for value in present):
        return 1
    return -3


def _align_profile_sequence(
    rows: Sequence[Sequence[Optional[str]]],
    sequence: Sequence[str],
    gap_score: int = -2,
) -> List[List[Optional[str]]]:
    """Needleman-Wunsch align one token sequence to an existing MSA profile."""
    if not rows:
        return [list(sequence)]
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("MSA profile rows have inconsistent widths")
    columns = [list(column) for column in zip(*rows)]
    # The score depends on membership, not multiplicity. Compute these once
    # per profile column instead of scanning every row at every DP cell.
    column_tokens = [{value for value in column if value is not None}
                     for column in columns]
    column_bases = [{_token_base(value) for value in values}
                    for values in column_tokens]
    n, m = len(columns), len(sequence)
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        scores[i][0] = scores[i - 1][0] + gap_score
        trace[i][0] = "U"
    for j in range(1, m + 1):
        scores[0][j] = scores[0][j - 1] + gap_score
        trace[0][j] = "L"
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            choices = (
                (scores[i - 1][j - 1] + (
                    4 if sequence[j - 1] in column_tokens[i - 1] else
                    1 if _token_base(sequence[j - 1]) in column_bases[i - 1] else -3
                ), "D"),
                (scores[i - 1][j] + gap_score, "U"),
                (scores[i][j - 1] + gap_score, "L"),
            )
            # D > U > L on tied scores keeps established block columns stable.
            scores[i][j], trace[i][j] = max(
                choices, key=lambda value: (value[0], "LUD".index(value[1]))
            )

    operations: List[str] = []
    i, j = n, m
    while i or j:
        operation = trace[i][j]
        if not operation:
            raise RuntimeError("failed to backtrace token MSA")
        operations.append(operation)
        if operation in {"D", "U"}:
            i -= 1
        if operation in {"D", "L"}:
            j -= 1
    operations.reverse()

    expanded = [[] for _row in rows]
    new_row: List[Optional[str]] = []
    profile_index = 0
    sequence_index = 0
    for operation in operations:
        if operation in {"D", "U"}:
            for row_index, row in enumerate(rows):
                expanded[row_index].append(row[profile_index])
            profile_index += 1
        else:
            for row in expanded:
                row.append(GAP)
        if operation in {"D", "L"}:
            new_row.append(sequence[sequence_index])
            sequence_index += 1
        else:
            new_row.append(GAP)
    return [*expanded, new_row]


def token_msa(sequences: Sequence[Sequence[str]]) -> List[List[Optional[str]]]:
    """Return a deterministic progressive multiple alignment of block paths."""
    nonempty = [list(map(str, sequence)) for sequence in sequences if sequence]
    if not nonempty:
        return []
    # Seed with the longest path, then add the most similar remaining path.
    seed_index = max(range(len(nonempty)), key=lambda index: (
        len(nonempty[index]), tuple(nonempty[index]), -index,
    ))
    seed = nonempty.pop(seed_index)
    profile: List[List[Optional[str]]] = [seed]
    seed_tokens = set(seed)
    nonempty.sort(key=lambda sequence: (
        -len(seed_tokens.intersection(sequence)), -len(sequence), tuple(sequence),
    ))
    for sequence in nonempty:
        profile = _align_profile_sequence(profile, sequence)
    return profile


def msa_linear_path(
    sequences: Sequence[Sequence[str]],
) -> Dict[str, int]:
    """Convert aligned token columns into one total block-index order."""
    msa = token_msa(sequences)
    if not msa:
        return {}
    positions: Dict[str, List[int]] = defaultdict(list)
    support = Counter()
    first_seen: Dict[str, Tuple[int, int]] = {}
    for row_index, row in enumerate(msa):
        for column_index, token in enumerate(row):
            if token is None:
                continue
            positions[token].append(column_index)
            support[token] += 1
            first_seen.setdefault(token, (row_index, column_index))

    def median_column(token: str) -> int:
        values = sorted(positions[token])
        return values[(len(values) - 1) // 2]

    ordered = sorted(positions, key=lambda token: (
        median_column(token),
        -support[token],
        first_seen[token],
        int(_token_base(token)) if _token_base(token).isdigit() else 10**18,
        token,
    ))
    return {token: index for index, token in enumerate(ordered)}


def getlinearpath(allchunks_onquery: Mapping[object, Sequence[Sequence[object]]]):
    """Legacy-compatible entry point using all sequences' block indexes."""
    sequences = [
        [str(item[0]) for item in chunks]
        for _query, chunks in sorted(
            allchunks_onquery.items(), key=lambda value: str(value[0])
        )
    ]
    return msa_linear_path(sequences)


def _legacy_candidates() -> Iterable[str]:
    configured = os.environ.get("MINSETREF_GRAPHTOLINEAR")
    if configured:
        yield configured
    here = Path(__file__).resolve().parent
    yield str(here / "graphtolinear.py")
    yield str(here.parent / "Graph" / "graphtolinear.py")
    yield str(here.parent / "Ctyper2" / "graphtolinear.py")


_CIGAR_OP = re.compile(r"(\d+)([HSMIDX=N])([A-Za-z]*)")
_COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhdbn",
)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def cigar_operations(cigar: str):
    """Read optional query payloads; missing bases remain missing.

    M and = mean exact matches in this format. X carries query bases (or
    legacy reference+query bases); D/H/N carry no query bases. A short
    payload represents a known prefix, never invented/padded sequence.
    """
    body = cigar.partition(":")[2] if ":" in cigar else cigar
    body = body.lstrip("><")
    matches = list(_CIGAR_OP.finditer(body))
    if "".join(match[0] for match in matches) != body:
        raise ValueError(f"invalid graph CIGAR: {cigar!r}")
    for match in matches:
        size, op, payload = int(match[1]), match[2], match[3]
        if op in "DHN":
            payload = ""
        elif op == "X" and size and len(payload) == 2 * size:
            payload = payload[size:]
        elif len(payload) > size:
            raise ValueError(f"payload exceeds operation length: {match[0]!r}")
        yield size, op, payload


def cigar_cutrrange(cigar: str, rstart: int, rend: int, strd: int = 1):
    """Legacy graph-coordinate clipping, retaining and slicing base payloads.

    Insertions at a cut belong to the segment on the left. Initial insertions
    are recovered by the query-gap code in findglobalcigar. H (and legacy S)
    follow the graph-CIGAR convention, not SAM clipping semantics.
    """
    operations = list(cigar_operations(cigar))
    if strd == -1:
        # A partial prefix becomes a suffix on reversal. Split it into an
        # unknown operation and a known operation so its position is retained.
        reversed_ops = []
        for size, op, payload in reversed(operations):
            if payload and len(payload) < size:
                reversed_ops.append((size - len(payload), op, ""))
                size = len(payload)
            reversed_ops.append((size, op, reverse_complement(payload)))
        operations = reversed_ops
    qstart = qend = -1
    rpos = qpos = 0
    output = []
    for size, op, payload in operations:
        old_r, old_q = rpos, qpos
        if op in "M=XS":
            rpos += size
            qpos += size
        elif op in "DHN":
            rpos += size
        elif op == "I":
            qpos += size
        if rpos <= rstart:
            continue
        if old_r > rend:
            break
        take = min(rpos, rend) - max(old_r, rstart)
        if take:
            offset = max(old_r, rstart) - old_r
            output.append(f"{take}{op}{payload[offset:offset + take]}")
        elif op == "I":
            output.append(f"{size}I{payload}")
        if qstart == -1 and rstart >= old_r:
            qstart = old_q + rstart - old_r if op in "M=XSI" else old_q
        if qend == -1 and rend < rpos:
            qend = old_q + rend - old_r if op in "M=XSI" else old_q
    if rend >= rpos:
        qend = qpos
        if rend > rpos:
            output.append(f"{rend - rpos}H")
    return qstart, qend, "".join(output)


def _gap_cigar(start, end, payload_intervals):
    """Retain known query bases in gaps while leaving unknown intervals bare."""
    output = []
    cursor = start
    for left, right, bases in payload_intervals:
        left, right = max(left, start), min(right, end)
        if right <= left:
            continue
        if left > cursor:
            output.append(f"{left - cursor}I")
        # The interval's original start is stored alongside its payload.
        offset = left - bases[0]
        output.append(f"{right - left}I{bases[1][offset:offset + right - left]}")
        cursor = right
    if cursor < end:
        output.append(f"{end - cursor}I")
    return "".join(output)


class _QueryCigar(str):
    """A CIGAR plus its query offset and oriented graph reference.

    A str subclass keeps the legacy projection interface and is picklable by
    multiprocessing. Reference bases are used only to preserve query bases
    that the graph segmentation moves into insertion gaps.
    """
    def __new__(cls, text, query_start=0, reference=None):
        value = super().__new__(cls, text)
        value.query_start = query_start
        value.reference = reference
        return value


def _fasta_records(path):
    name, pieces = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(pieces)
                name, pieces = line[1:].split()[0], []
            elif name is not None:
                pieces.append(line.strip())
    if name is not None:
        yield name, "".join(pieces)


def _query_cigar_inputs(queryfile, inputfile, reffile):
    references = dict(_fasta_records(reffile))
    reverse_references = {}
    lengths = {name: len(sequence) for name, sequence in _fasta_records(queryfile)}
    reads = {}
    with open(inputfile) as handle:
        for raw in handle:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 5 or fields[1] in ("", "*"):
                continue
            name, _path, cigar, _ranges, positions = fields
            parts = re.findall(r"([<>])([^:<>]+):([^<>]+)", cigar)
            if "".join(f"{d}{n}:{c}" for d, n, c in parts) != cigar:
                raise ValueError(f"{name}: invalid named graph CIGAR")
            if "_" in positions:
                spans = [tuple(map(int, span.split("_"))) for span in positions.split(";")]
            else:
                # Historical five-column files may have flattened intervals.
                values = list(map(int, positions.split(";")))
                spans = list(zip(values[::2], values[1::2]))
            if len(spans) != len(parts) or any(len(span) != 2 for span in spans):
                raise ValueError(f"{name}: query spans do not match CIGAR components")
            reads[name] = []
            previous = 0
            for (direction, node, body), (start, end) in zip(parts, spans):
                consumed = sum(size for size, op, _ in cigar_operations(body) if op in "M=XSI")
                if start < previous or end - start != consumed or end > lengths[name]:
                    raise ValueError(f"{name}: inconsistent query span {start}_{end}")
                previous = end
                reference = references[node]
                if direction == "<":
                    if node not in reverse_references:
                        reverse_references[node] = reverse_complement(reference)
                    reference = reverse_references[node]
                reads[name].append(_QueryCigar(
                    f"{direction}{node}:{body}", start, reference,
                ))
    return reads, lengths, list(lengths)


def findglobalcigar(queryaligns, chunks, querylength, queryread, chunkread):
    """Project existing CIGARs, preserving payloads even in boundary gaps.

    This importable worker also works with multiprocessing's spawn start
    method; workers no longer need to import the dynamically loaded module.
    Query FASTA bases are deliberately unused in this non-polishing mode.
    """
    query_starts, payload_intervals = [], []
    query_pos = 0
    for alignment in queryaligns:
        query_pos = getattr(alignment, "query_start", query_pos)
        reference = getattr(alignment, "reference", None)
        reference_pos = 0
        query_starts.append(query_pos)
        for size, op, payload in cigar_operations(alignment):
            if op in "M=XSI":
                if op in "M=" and reference is not None:
                    payload = reference[reference_pos:reference_pos + size]
                    if len(payload) != size:
                        raise ValueError("CIGAR exceeds graph reference")
                if payload:
                    payload_intervals.append(
                        (query_pos, query_pos + len(payload), (query_pos, payload))
                    )
                query_pos += size
            if op in "M=XSDHN":
                reference_pos += size
    reference_ends = [0] * len(queryaligns)
    cigars = []
    last_query_end = 0
    for name, start, end, strand, align_index, _size, _index in chunks:
        if strand == -1:
            reference = getattr(queryaligns[align_index], "reference", None)
            pathsize = (len(reference) if reference is not None else
                        int(name.split("_")[-1]) - int(name.split("_")[-2]))
            start, end = pathsize - end, pathsize - start
        previous = reference_ends[align_index]
        qstart, qend, clipped = cigar_cutrrange(
            queryaligns[align_index], max(start, previous), end,
        )
        expected = query_starts[align_index] + qstart
        prefix = f"{previous - start}H" if previous > start else ""
        if expected > last_query_end:
            gap = _gap_cigar(last_query_end, expected, payload_intervals)
            if cigars:
                cigars[-1] += gap
            else:
                prefix += gap
        cigars.append(queryaligns[align_index][0] + prefix + clipped)
        reference_ends[align_index] = end
        last_query_end = query_starts[align_index] + qend
    if cigars and querylength > last_query_end:
        cigars[-1] += _gap_cigar(last_query_end, querylength, payload_intervals)
    return 0, 0, cigars


def _load_legacy():
    source = next(
        (os.path.abspath(path) for path in _legacy_candidates() if os.path.isfile(path)),
        None,
    )
    if source is None:
        raise FileNotFoundError(
            "legacy graphtolinear.py was not found; set MINSETREF_GRAPHTOLINEAR"
        )
    existing = sys.modules.get(LEGACY_MODULE_NAME)
    if (existing is not None and os.path.abspath(
            getattr(existing, "__file__", "")) == source):
        existing.getlinearpath = getlinearpath
        existing.findglobalcigar = findglobalcigar
        return existing

    spec = importlib.util.spec_from_file_location(LEGACY_MODULE_NAME, source)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load legacy graphtolinear module: {source}")
    module = importlib.util.module_from_spec(spec)
    # multiprocessing pickles worker functions by ``module.function``.  The
    # module therefore has to be registered before exec_module defines
    # findglobalalign/findglobalcigar and before the legacy code creates a
    # worker Pool.  Without this entry, ForkingPickler reports that importing
    # minsetref_legacy_graphtolinear failed.
    sys.modules[LEGACY_MODULE_NAME] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        if sys.modules.get(LEGACY_MODULE_NAME) is module:
            del sys.modules[LEGACY_MODULE_NAME]
        raise
    module.getlinearpath = getlinearpath
    module.findglobalcigar = findglobalcigar
    return module


def graphtolinear(
    inputfile: str,
    queryfile: str,
    reffile: str,
    output: str,
    thread: int = 8,
    tempfolder: str = "",
    polish: int = 0,
):
    """Run legacy projection/output with the MSA block-order callback."""
    normalized = inputfile
    temporary: Optional[str] = None
    query_metadata = {}
    with open(inputfile, "rt") as handle:
        first = next((line for line in handle if line.strip() and not line.startswith("#")), "")
    if len(first.rstrip("\r\n").split("\t")) == 10:
        descriptor, temporary = tempfile.mkstemp(
            prefix=".localgraphtolinear.", suffix=".five-column.tsv",
            dir=os.path.dirname(os.path.abspath(output)) or ".",
        )
        os.close(descriptor)
        best: Dict[str, Tuple[int, Tuple[str, str, str, str, str]]] = {}
        with open(inputfile, "rt") as source:
            for raw in source:
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) != 10 or fields[0] == "hotspot_index":
                    continue
                span = int(fields[3]) - int(fields[2])
                row = (fields[5], fields[6], fields[7], fields[8], fields[9])
                if row[1] == "*":
                    row = (row[0], "", "", "", "")
                if fields[5] not in best or span > best[fields[5]][0]:
                    best[fields[5]] = (span, row)
                    query_metadata[fields[5]] = (fields[4], fields[7] if fields[6] == "*" else "")
        with open(temporary, "wt") as converted:
            for _name, (_span, row) in best.items():
                converted.write("\t".join(row) + "\n")
        normalized = temporary
    try:
        legacy = _load_legacy()
        original_reader = legacy.getquerycigar
        legacy.getquerycigar = partial(_query_cigar_inputs, reffile=reffile)
        try:
            result = legacy.graphtolinear(
                normalized, queryfile, reffile, output, thread, tempfolder, polish,
            )
        finally:
            legacy.getquerycigar = original_reader
        if query_metadata:
            # Keep the original nine columns used by visualization scripts.
            # QS records whole-query orientation; individual < / > markers
            # still describe each segment relative to its stored S sequence.
            with open(output, "rt") as handle:
                lines = handle.readlines()
            with open(output, "wt") as handle:
                for line in lines:
                    fields = line.rstrip("\r\n").split("\t")
                    if fields[0] == "L" and fields[1] in query_metadata:
                        strand, unaligned_cigar = query_metadata[fields[1]]
                        fields.append(f"QS:Z:{strand}")
                        if unaligned_cigar and not polish:
                            fields.append(f"UC:Z:{unaligned_cigar}")
                        line = "\t".join(fields) + "\n"
                    handle.write(line)
        return result
    finally:
        if temporary:
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="five- or ten-column graph alignment; optional I/X payloads")
    parser.add_argument("-q", "--query", required=True, help="query FASTA")
    parser.add_argument("-r", "--ref", required=True, help="graph .FA")
    parser.add_argument("-o", "--output", required=True, help="linear GAF output")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-f", "--folder", dest="tempfolder", default="")
    parser.add_argument("-p", "--polish", type=int, default=0,
                        help="0: preserve input CIGARs/payloads (default); 1: legacy realignment")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.threads < 1:
        raise ValueError("--threads must be positive")
    graphtolinear(
        args.input, args.query, args.ref, args.output,
        args.threads, args.tempfolder, args.polish,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
