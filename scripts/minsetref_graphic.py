#!/usr/bin/env python3
"""Graphic-path annotations for query sequence placed on original FASTAs."""
from __future__ import annotations

from typing import List, Mapping, Optional, Sequence, Tuple

from minsetref_core import IndexedFasta


PairedBlock = Tuple[int, int, int, int]
_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _append_cigar_op(
    operations: List[Tuple[int, str]],
    length: int,
    operation: str,
) -> None:
    if length <= 0:
        return
    if operations and operations[-1][1] == operation:
        prior_length, _prior_operation = operations[-1]
        operations[-1] = (prior_length + length, operation)
    else:
        operations.append((length, operation))


def _exact_base_operations(
    query_sequence: str,
    reference_sequence: str,
) -> List[Tuple[int, str]]:
    if len(query_sequence) != len(reference_sequence):
        raise ValueError(
            "paired alignment block has unequal query/reference sequence "
            f"lengths: {len(query_sequence)} != {len(reference_sequence)}"
        )
    operations: List[Tuple[int, str]] = []
    for query_base, reference_base in zip(
        query_sequence, reference_sequence,
    ):
        _append_cigar_op(
            operations,
            1,
            "="
            if query_base.upper() == reference_base.upper()
            else "X",
        )
    return operations


def _coalesce_collinear_blocks(
    paired_blocks: Sequence[PairedBlock],
    strand: str,
    reference_length: int,
) -> List[PairedBlock]:
    """Losslessly merge duplicate/overlapping blocks on one diagonal."""
    oriented: List[Tuple[int, int, int, int]] = []
    for q0, q1, ref0, ref1 in paired_blocks:
        if strand == "+":
            oriented_ref0, oriented_ref1 = ref0, ref1
        else:
            oriented_ref0 = reference_length - ref1
            oriented_ref1 = reference_length - ref0
        oriented.append((q0, q1, oriented_ref0, oriented_ref1))
    oriented.sort(key=lambda row: (row[0], row[1], row[2], row[3]))

    merged: List[List[int]] = []
    for q0, q1, ref0, ref1 in oriented:
        if merged:
            prior = merged[-1]
            same_diagonal = ref0 - q0 == prior[2] - prior[0]
            if same_diagonal and q0 <= prior[1]:
                new_end = max(prior[1], q1)
                prior[1] = new_end
                prior[3] = prior[2] + new_end - prior[0]
                continue
        merged.append([q0, q1, ref0, ref1])

    output: List[PairedBlock] = []
    for q0, q1, oriented_ref0, oriented_ref1 in merged:
        if strand == "+":
            ref0, ref1 = oriented_ref0, oriented_ref1
        else:
            ref0 = reference_length - oriented_ref1
            ref1 = reference_length - oriented_ref0
        output.append((q0, q1, ref0, ref1))
    return output


def build_graphic_path(
    *,
    query_id: str,
    query_start: int,
    query_end: int,
    query_sequence: str,
    ref_contig: str,
    strand: str,
    paired_blocks: Sequence[PairedBlock],
    reference: IndexedFasta,
) -> str:
    """Encode accepted query/reference blocks on the original reference."""
    if not paired_blocks:
        return "."
    if ref_contig not in set(reference.names()):
        raise ValueError(
            f"{query_id}: original reference lacks contig {ref_contig!r}"
        )
    if len(query_sequence) != query_end - query_start:
        raise ValueError(
            f"{query_id}: query sequence length {len(query_sequence)} does "
            f"not match interval {query_start}-{query_end}"
        )
    if strand not in {"+", "-"}:
        raise ValueError(f"{query_id}: invalid strand {strand!r}")

    reference_length = reference.length(ref_contig)
    marker = ">" if strand == "+" else "<"
    rendered: List[str] = []
    current_ops: List[Tuple[int, str]] = []
    last_query_end = query_start
    last_oriented_ref_end: Optional[int] = None

    def close_segment() -> None:
        nonlocal current_ops, last_oriented_ref_end
        if last_oriented_ref_end is None:
            return
        _append_cigar_op(
            current_ops,
            reference_length - last_oriented_ref_end,
            "H",
        )
        rendered.append(
            marker
            + ref_contig
            + ":"
            + "".join(f"{length}{op}" for length, op in current_ops)
        )
        current_ops = []
        last_oriented_ref_end = None

    normalized_blocks = _coalesce_collinear_blocks(
        paired_blocks, strand, reference_length,
    )
    for q0, q1, ref0, ref1 in normalized_blocks:
        if (
            q0 < query_start
            or q1 > query_end
            or q1 <= q0
            or ref0 < 0
            or ref1 > reference_length
            or ref1 <= ref0
            or q1 - q0 != ref1 - ref0
        ):
            raise ValueError(
                f"{query_id}: invalid original-reference paired block "
                f"{q0}-{q1}:{ref0}-{ref1}"
            )
        if q0 < last_query_end:
            raise ValueError(
                f"{query_id}: conflicting accepted query blocks overlap in "
                "graphic liftover; the blocks map the same query bases to "
                "different reference coordinates"
            )
        if strand == "+":
            oriented_ref_start, oriented_ref_end = ref0, ref1
            reference_piece = reference.fetch(ref_contig, ref0, ref1)
        else:
            oriented_ref_start = reference_length - ref1
            oriented_ref_end = reference_length - ref0
            reference_piece = _reverse_complement(
                reference.fetch(ref_contig, ref0, ref1)
            )

        query_piece = query_sequence[
            q0 - query_start:q1 - query_start
        ]
        exact_ops = _exact_base_operations(query_piece, reference_piece)
        query_gap = q0 - last_query_end
        if (
            last_oriented_ref_end is None
            or oriented_ref_start < last_oriented_ref_end
        ):
            close_segment()
            _append_cigar_op(current_ops, oriented_ref_start, "H")
            _append_cigar_op(current_ops, query_gap, "I")
        else:
            _append_cigar_op(
                current_ops,
                oriented_ref_start - last_oriented_ref_end,
                "D",
            )
            _append_cigar_op(current_ops, query_gap, "I")
        for length, operation in exact_ops:
            _append_cigar_op(current_ops, length, operation)
        last_query_end = q1
        last_oriented_ref_end = oriented_ref_end

    _append_cigar_op(current_ops, query_end - last_query_end, "I")
    close_segment()
    return "".join(rendered) if rendered else "."


def build_multilocus_graphic_path(
    *,
    query_id: str,
    query_start: int,
    query_end: int,
    query_sequence: str,
    components: Sequence[object],
    references: Mapping[Tuple[str, str], IndexedFasta],
    exact_bases: bool = True,
) -> str:
    """Encode an ordered query path that can switch reference locations.

    Each component must expose ``ref_haplotype``, ``ref_contig``, ``strand``,
    ``query_start``, ``query_end`` and ``paired_blocks`` attributes. Query
    bases between selected paired blocks or target switches are represented by
    ``I``. A target-coordinate reversal closes the current graph segment and
    starts another marker, preserving inversions, translocations, and other
    structural transitions instead of rejecting them. Unlike the legacy
    graphcigarlight output, the first leading query flank is also emitted as
    ``I``, making the graphic CIGAR self-contained without a separate
    qpositions column. When ``exact_bases`` is false, aligned paired blocks
    are emitted as ``M`` directly from their intervals; reference/query bases
    are not fetched or compared.
    """
    if not components:
        return "."
    if exact_bases and len(query_sequence) != query_end - query_start:
        raise ValueError(
            f"{query_id}: query sequence length {len(query_sequence)} does "
            f"not match interval {query_start}-{query_end}"
        )
    ordered = sorted(
        components,
        key=lambda row: (
            int(getattr(row, "query_start")),
            int(getattr(row, "query_end")),
            str(getattr(row, "ref_haplotype")),
            str(getattr(row, "ref_contig")),
            str(getattr(row, "strand")),
        ),
    )
    rendered: List[str] = []
    query_cursor = query_start
    for component_index, component in enumerate(ordered):
        haplotype = str(getattr(component, "ref_haplotype"))
        contig = str(getattr(component, "ref_contig"))
        strand = str(getattr(component, "strand"))
        if strand not in {"+", "-"}:
            raise ValueError(
                f"{query_id}: invalid component strand {strand!r}"
            )
        try:
            reference = references[(haplotype, contig)]
        except KeyError as error:
            raise ValueError(
                f"{query_id}: no reference reader for {haplotype}:{contig}"
            ) from error
        reference_length = reference.length(contig)
        blocks = sorted(
            tuple(getattr(component, "paired_blocks")),
            key=lambda block: (block[0], block[1], block[2], block[3]),
        )
        if not blocks:
            continue
        component_start = int(getattr(component, "query_start"))
        component_end = int(getattr(component, "query_end"))
        if component_start < query_cursor or component_end <= component_start:
            raise ValueError(
                f"{query_id}: overlapping/out-of-order multi-locus component "
                f"{component_start}-{component_end} after {query_cursor}"
            )
        marker = ">" if strand == "+" else "<"
        operations: List[Tuple[int, str]] = []
        last_query = query_cursor
        last_oriented_ref: Optional[int] = None

        def close_segment() -> None:
            nonlocal operations, last_oriented_ref
            if last_oriented_ref is None:
                return
            _append_cigar_op(
                operations, reference_length - last_oriented_ref, "H",
            )
            rendered.append(
                marker + contig + ":"
                + "".join(
                    f"{length}{operation}"
                    for length, operation in operations
                )
            )
            operations = []
            last_oriented_ref = None

        for q0, q1, ref0, ref1 in blocks:
            if (
                q0 < component_start
                or q1 > component_end
                or q1 <= q0
                or ref0 < 0
                or ref1 > reference_length
                or ref1 <= ref0
                or q1 - q0 != ref1 - ref0
            ):
                raise ValueError(
                    f"{query_id}: invalid multi-locus paired block "
                    f"{q0}-{q1}:{ref0}-{ref1}"
                )
            if strand == "+":
                oriented_ref_start, oriented_ref_end = ref0, ref1
                reference_piece = (
                    reference.fetch(contig, ref0, ref1)
                    if exact_bases else None
                )
            else:
                oriented_ref_start = reference_length - ref1
                oriented_ref_end = reference_length - ref0
                reference_piece = (
                    _reverse_complement(
                        reference.fetch(contig, ref0, ref1)
                    )
                    if exact_bases else None
                )
            if last_oriented_ref is None:
                _append_cigar_op(operations, oriented_ref_start, "H")
            else:
                if oriented_ref_start < last_oriented_ref:
                    close_segment()
                    _append_cigar_op(
                        operations, oriented_ref_start, "H",
                    )
                else:
                    _append_cigar_op(
                        operations,
                        oriented_ref_start - last_oriented_ref,
                        "D",
                    )
            _append_cigar_op(operations, q0 - last_query, "I")
            if exact_bases:
                query_piece = query_sequence[
                    q0 - query_start:q1 - query_start
                ]
                for length, operation in _exact_base_operations(
                    query_piece, reference_piece,
                ):
                    _append_cigar_op(operations, length, operation)
            else:
                _append_cigar_op(operations, q1 - q0, "M")
            last_query = q1
            last_oriented_ref = oriented_ref_end
        if component_index == len(ordered) - 1:
            _append_cigar_op(operations, query_end - last_query, "I")
            last_query = query_end
        close_segment()
        query_cursor = last_query
    return "".join(rendered) if rendered else "."
