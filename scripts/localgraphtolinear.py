#!/usr/bin/env python3
"""Linearize one local graph using an MSA of all block-index paths.

The legacy ``graphtolinear.py`` graph/CIGAR projection code is retained, but
its heuristic graph walk is replaced at runtime by the deterministic
progressive token-MSA implemented here.  Tokens are graph block indexes (for
example ``12`` and the second repeated occurrence ``12_2``), not nucleotide
bases.
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import tempfile
from collections import Counter, defaultdict
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
                (scores[i - 1][j - 1] + _column_score(columns[i - 1], sequence[j - 1]), "D"),
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
        with open(temporary, "wt") as converted:
            for _name, (_span, row) in best.items():
                converted.write("\t".join(row) + "\n")
        normalized = temporary
    try:
        return _load_legacy().graphtolinear(
            normalized, queryfile, reffile, output, thread, tempfolder, polish,
        )
    finally:
        if temporary:
            try:
                os.remove(temporary)
            except FileNotFoundError:
                pass


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="five-column graph alignment")
    parser.add_argument("-q", "--query", required=True, help="query FASTA")
    parser.add_argument("-r", "--ref", required=True, help="graph .FA")
    parser.add_argument("-o", "--output", required=True, help="linear GAF output")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("-f", "--folder", dest="tempfolder", default="")
    parser.add_argument("-p", "--polish", type=int, default=0)
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
