"""Nucleotide alignment scores used to rank and accept candidates.

Native aligners still generate candidates with their own search objectives.
Legacy scoring and trimming use +4 match, -8 mismatch, and a cost of
12 + 100*log2(n) for every contiguous insertion/deletion run of length n.
The one-base gap has no extension penalty and costs only the gap open, 12.
D/I sequence repair instead uses +1 match, -4 mismatch, and an affine gap
cost of 4 + (n - 1), exposed separately by di_repair_score.
Whole-pair tandem swaps use whole_pair_score (+4/-8, gaps 24+8*n).
SV-window swaps use MaskedSVScorer (+4/-8, gaps 24+8*n+160*log2(masked+1)).
M is treated as match-like only when the caller has no mismatch information;
use explicit =/X or supply a known match count for SAM/PAF CIGARs.
"""
from dataclasses import dataclass
from array import array
import math
import re
from typing import Iterable, Tuple

MATCH_SCORE = 4
MISMATCH_SCORE = -8
GAP_OPEN = 12
GAP_LOG_EXTENSION = 100
DI_MATCH_SCORE = 1
DI_MISMATCH_SCORE = -4
DI_GAP_OPEN = 4
DI_GAP_EXTEND = 1
_CIGAR_RE = re.compile(r'(\d+)([=MXIDHSNP])([A-Za-z]*)')


def gap_extension_penalty(size: int) -> float:
    if size <= 0:
        raise ValueError('gap length must be positive')
    return GAP_LOG_EXTENSION * math.log2(size)


def gap_penalty(size: int) -> float:
    return GAP_OPEN + gap_extension_penalty(size)


def operation_parts(operation):
    if hasattr(operation, 'n'):
        return int(operation.n), operation.op
    return int(operation[0]), operation[1]


def cigar_operations(body: str):
    """Parse an unnamed pairwise CIGAR, including optional I/X base payloads."""
    if ':' in body or '>' in body or '<' in body:
        raise ValueError('remove graph segment names before scoring a CIGAR')
    return [(int(m[1]), m[2]) for m in _CIGAR_RE.finditer(body)]


def scoring_runs(operations: Iterable):
    """Yield (first index, exclusive end index, length, op) for scored runs.

    Adjacent I tokens or D tokens count as one gap, including across H padding.
    Original indices are retained so trimming never cuts a scored gap in half.
    Non-scored clips/skips separate runs except for coordinate-only H padding.
    """
    pending = None
    for index, operation in enumerate(operations):
        size, op = operation_parts(operation)
        if size < 0:
            raise ValueError('CIGAR operation length must be non-negative')
        if not size or op == 'H':
            continue
        if op not in '=MXIDSNP':
            raise ValueError(f'unsupported alignment operation: {op}')
        op = '=' if op == 'M' else op
        if pending is not None and op == pending[3]:
            pending = (pending[0], index + 1, pending[2] + size, op)
        else:
            if pending is not None:
                yield pending
            pending = (index, index + 1, size, op)
    if pending is not None:
        yield pending


def run_score(size: int, op: str) -> float:
    if op in '=M':
        return MATCH_SCORE * size
    if op == 'X':
        return MISMATCH_SCORE * size
    if op in 'ID':
        return -gap_penalty(size)
    return 0.0


@dataclass(frozen=True)
class AlignmentScoreComponents:
    matches: int
    mismatches: int
    gap_lengths: Tuple[int, ...]

    @property
    def gap_opens(self) -> int:
        return len(self.gap_lengths)

    @property
    def gap_extensions(self) -> int:
        """Raw n-1 count retained for diagnostics; NOT the scoring penalty."""
        return sum(n - 1 for n in self.gap_lengths)

    @property
    def extension_penalty(self) -> float:
        return math.fsum(gap_extension_penalty(n) for n in self.gap_lengths)

    @property
    def score(self) -> float:
        # Summing complete costs avoids order/tokenization-dependent ties.
        return (MATCH_SCORE * self.matches + MISMATCH_SCORE * self.mismatches
                - math.fsum(gap_penalty(n) for n in self.gap_lengths))


def alignment_score_components(operations: Iterable) -> AlignmentScoreComponents:
    matches = mismatches = 0
    gaps = []
    for _first, _end, size, op in scoring_runs(operations):
        if op == '=':
            matches += size
        elif op == 'X':
            mismatches += size
        elif op in 'ID':
            gaps.append(size)
    return AlignmentScoreComponents(matches, mismatches, tuple(gaps))


def alignment_score(operations: Iterable) -> float:
    return alignment_score_components(operations).score


def whole_pair_score(operations: Iterable) -> int:
    """Whole-pair tandem swap: +4/-8 bases, gap cost 24 + 8*n per run."""
    components = alignment_score_components(operations)
    return (4 * components.matches - 8 * components.mismatches
            - sum(24 + 8 * n for n in components.gap_lengths))


class MaskedSVScorer:
    """Score verified SV-window CIGARs against oriented, case-preserved axes.

    Each gap costs 24 + 8*n + 160*log2(masked+1), counting lowercase bases
    in the gap on the query for I and reference for D.
    Coordinates passed to __call__ are relative to these full context strings.
    """
    def __init__(self, reference: str, query: str):
        def prefix(sequence):
            result = array('Q', [0])
            for base in sequence:
                result.append(result[-1] + base.islower())
            return result
        self.reference_masks = prefix(reference)
        self.query_masks = prefix(query)

    @staticmethod
    def gap_cost(prefix, start, size):
        end = start + size
        length = len(prefix) - 1
        if not 0 <= start < end <= length:
            raise ValueError('gap is outside its scoring sequence')
        masked = prefix[end] - prefix[start]
        return 24 + 8 * size + 160 * math.log2(masked + 1)

    def __call__(self, operations: Iterable, reference_start=0, query_start=0):
        r, q = reference_start, query_start
        if r < 0 or q < 0:
            raise ValueError('negative scoring sequence offset')
        matches = mismatches = 0
        costs = []
        for _first, _end, size, op in scoring_runs(operations):
            if op == '=':
                matches += size
            elif op == 'X':
                mismatches += size
            elif op == 'I':
                costs.append(self.gap_cost(self.query_masks, q, size))
            elif op == 'D':
                costs.append(self.gap_cost(self.reference_masks, r, size))
            else:
                raise ValueError(f'unsupported SV scoring operation: {op}')
            if op in '=XD':
                r += size
            if op in '=XI':
                q += size
        if r >= len(self.reference_masks) or q >= len(self.query_masks):
            raise ValueError('CIGAR exceeds its scoring sequences')
        return 4 * matches - 8 * mismatches - math.fsum(costs)


def di_repair_score(operations: Iterable) -> int:
    """Score D/I sequence repair; a one-base gap costs 4, each extra base 1.

    As with the SV score, join adjacent same-axis gap tokens before charging
    openings. Callers must resolve M and verify =/X against the actual bases.
    """
    components = alignment_score_components(operations)
    return (DI_MATCH_SCORE * components.matches
            + DI_MISMATCH_SCORE * components.mismatches
            - DI_GAP_OPEN * components.gap_opens
            - DI_GAP_EXTEND * components.gap_extensions)


def cigar_score(body: str, *, matches=None, mismatches=None) -> float:
    """Score a pairwise CIGAR, optionally resolving M using SAM/PAF counts."""
    components = alignment_score_components(cigar_operations(body))
    if matches is None and mismatches is None:
        return components.score
    aligned = components.matches + components.mismatches
    if matches is None:
        matches = aligned - int(mismatches)
    matches = int(matches)
    if not 0 <= matches <= aligned:
        raise ValueError('match count is outside the aligned CIGAR span')
    return AlignmentScoreComponents(matches, aligned - matches,
                                    components.gap_lengths).score
