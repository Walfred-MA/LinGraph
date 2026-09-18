#!/usr/bin/env python3
"""Discover and index local references for graphs lacking the chosen genome.

Each selected alternative-locus path is an independent reference contig.  Its
public graph/FASTA path ID is retained as the contig name, its coordinates are
local ``0..length`` coordinates, and its sequence remains in that FASTA path's
orientation.  Original assembly coordinates are provenance metadata only.
"""

from __future__ import annotations

import argparse
import dataclasses
import os
from pathlib import Path
import re
import sys
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from align_partition_hotspots import graph_name_from_list_entry
from block_partition_alignments import graph_prefix
from minsetref_core import sanitize_id
from summarize_partition_hotspot_segments import read_alignment_output
from assembly_contigs import accepted_contig


PROTOCOL = "local-reference-templates-v2"
LEGACY_PROTOCOL = "local-reference-templates-v1"
SOURCE_COORDINATE_RE = re.compile(
    r"^([^:]+):(.+):(\d+)-(\d+)([+-])$"
)
LOCAL_COORDINATE_RE = re.compile(r"^(.+):(\d+)-(\d+)([+-])$")


@dataclasses.dataclass(frozen=True)
class LocalTemplate:
    name: str
    graph_name: str
    graph_prefix: str
    graph_path: str
    source_haplotype: str
    source_contig: str
    output_contig: str
    start: int
    end: int
    source_strand: str
    source_start: int = 0
    source_end: int = 0
    # Orientation of the physical checkpoint FASTA record relative to the
    # public alternative-locus contig.  v2 always writes "+".  This remains
    # explicit so a v1 catalog can be consumed without regeneration.
    storage_strand: str = "+"
    sequence: str = ""
    backbone: str = ""
    lift_status: str = ""
    reference: str = "."
    lift_note: str = "."
    fixed: bool = False


def alternative_locus_name(graph_path: str) -> str:
    """Return the independent public contig ID for an alternative locus."""
    if (
        not graph_path
        or any(character.isspace() for character in graph_path)
        or graph_path.startswith((">", "<"))
        or "/" in graph_path
    ):
        raise ValueError(f"invalid alternative-locus FASTA ID {graph_path!r}")
    return graph_path


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return sequence.translate(table)[::-1]


def iter_fasta(path: str) -> Iterator[Tuple[Tuple[str, ...], str]]:
    fields: Optional[Tuple[str, ...]] = None
    pieces: List[str] = []
    with open(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if fields is not None:
                    yield fields, "".join(pieces)
                fields = tuple(raw[1:].strip().split())
                if not fields:
                    raise ValueError(f"{path}: empty FASTA header")
                pieces = []
            elif raw.strip() and not raw.startswith("#"):
                if fields is None:
                    raise ValueError(f"{path}: sequence before FASTA header")
                pieces.append(raw.strip())
    if fields is not None:
        yield fields, "".join(pieces)


def _header_attributes(fields: Sequence[str]) -> Dict[str, str]:
    return {
        key: value
        for field in fields[1:] if "=" in field
        for key, value in (field.split("=", 1),)
    }


def _is_reference_template(fields: Sequence[str]) -> bool:
    attributes = _header_attributes(fields)
    return (
        attributes.get("class", "").lower() == "reference"
        or attributes.get("role", "").lower() == "reference"
        or any(field.lower() == "reference" for field in fields[2:])
    )


def templates_from_graph(
    graph_name: str,
    graph_fasta: str,
) -> List[LocalTemplate]:
    """Read reference-class paths and their original source coordinates."""
    graph_key = sanitize_id(graph_prefix(graph_name))
    templates: List[LocalTemplate] = []
    for fields, sequence in iter_fasta(graph_fasta):
        if not _is_reference_template(fields) or len(fields) < 2:
            continue
        source = SOURCE_COORDINATE_RE.fullmatch(fields[1])
        if source is None:
            continue
        haplotype, source_contig, start_text, end_text, strand = source.groups()
        if not accepted_contig(haplotype, source_contig):
            continue
        start, end = int(start_text), int(end_text)
        if end <= start:
            raise ValueError(
                f"{graph_fasta}: invalid template interval {fields[1]!r}"
            )
        if len(sequence) != end - start:
            raise ValueError(
                f"{graph_fasta}: template path {fields[0]!r} has "
                f"{len(sequence)} bases but source interval spans {end - start}"
            )

        # The first FASTA field is the alternative locus's public identity.
        # The source interval records provenance only; it must not leak back
        # into VCF CHROM/POS or change the locus FASTA's orientation.
        output_contig = alternative_locus_name(fields[0])
        locus_length = len(sequence)
        templates.append(LocalTemplate(
            name=output_contig,
            graph_name=graph_name,
            graph_prefix=graph_key,
            graph_path=fields[0],
            source_haplotype=haplotype,
            source_contig=source_contig,
            output_contig=output_contig,
            start=0,
            end=locus_length,
            source_strand=strand,
            source_start=start,
            source_end=end,
            storage_strand="+",
            sequence=sequence,
            fixed='sequence_role=imported_alternative' in fields,
        ))
    return templates


def read_graph_names(path: str) -> List[str]:
    output: List[str] = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            output.append(graph_name_from_list_entry(
                raw.strip().split()[0], f"{path}:{line_number}",
            ))
    if not output:
        raise ValueError(f"graph list is empty: {path}")
    return output


def discover_templates(
    graph_folder: str,
    graph_list: str,
    reference_alignment: str,
    reference_haplotype: str = "",
) -> List[LocalTemplate]:
    """Return templates only for graphs with no reference alignment row."""
    graph_names = read_graph_names(graph_list)
    reference_hotspots = {
        row.hotspot_index
        for row in read_alignment_output(reference_alignment)
        if (
            1 <= row.hotspot_index <= len(graph_names)
            and row.graph_path
            and row.graph_path != "*"
            and accepted_contig(reference_haplotype, row.contig)
        )
    }
    output: List[LocalTemplate] = []
    missing_graph_fastas: List[str] = []
    for hotspot_index, graph_name in enumerate(graph_names, 1):
        if hotspot_index in reference_hotspots:
            continue
        graph_fasta = os.path.join(
            graph_folder, graph_name, graph_name + ".FA",
        )
        if not os.path.isfile(graph_fasta):
            missing_graph_fastas.append(graph_fasta)
            continue
        output.extend(templates_from_graph(graph_name, graph_fasta))
    if missing_graph_fastas:
        examples = ", ".join(missing_graph_fastas[:5])
        remainder = (
            f" (and {len(missing_graph_fastas) - 5} more)"
            if len(missing_graph_fastas) > 5 else ""
        )
        print(
            "[local_reference_templates] warning: "
            f"{len(missing_graph_fastas)} graph path(s) selected by -L do "
            f"not exist under -G and will be skipped: {examples}{remainder}",
            file=sys.stderr,
        )
    output.sort(key=lambda row: (
        row.output_contig, row.start, row.end, row.graph_name, row.graph_path,
    ))
    return output


def _metadata_fields(template: LocalTemplate) -> str:
    values = {
        "protocol": PROTOCOL,
        "graph": template.graph_name,
        "graph_prefix": template.graph_prefix,
        "graph_path": template.graph_path,
        "source_haplotype": template.source_haplotype,
        "source_contig": template.source_contig,
        "source_start": str(template.source_start),
        "source_end": str(template.source_end),
        "source_strand": template.source_strand,
        "source": f"{template.source_haplotype}:{template.source_contig}:{template.source_start}-{template.source_end}{template.source_strand}",
    }
    if template.fixed:
        values['sequence_role'] = 'imported_alternative'
    if template.backbone:
        values.update(backbone=template.backbone, lift_status=template.lift_status,
                      reference=template.reference, lift_note=template.lift_note.replace(' ', '_'))
    return " ".join(f"{key}={value}" for key, value in values.items())


def write_templates(path: str, templates: Sequence[LocalTemplate]) -> None:
    """Write interval FASTA plus a faidx-compatible index atomically."""
    destination = Path(path).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + f".tmp.{os.getpid()}")
    index_temporary = Path(str(temporary) + ".fai")
    try:
        with temporary.open("wb") as output, index_temporary.open(
            "wt", encoding="utf-8",
        ) as index:
            output.write(f"#{PROTOCOL}\n".encode("utf-8"))
            index.write(f"#{PROTOCOL}\n")
            for template in templates:
                if template.storage_strand != "+":
                    raise ValueError(
                        f"{template.name}: v2 template storage must use the "
                        "alternative-locus FASTA orientation"
                    )
                if len(template.sequence) != template.end - template.start:
                    raise ValueError(
                        f"{template.name}: template sequence length does not "
                        "match its local coordinate span"
                    )
                coordinate = (
                    f"{template.output_contig}:{template.start}-"
                    f"{template.end}+"
                )
                header = (
                    f">{template.name} {coordinate} "
                    f"{_metadata_fields(template)}\n"
                ).encode("utf-8")
                output.write(header)
                offset = output.tell()
                sequence = template.sequence.encode("ascii")
                output.write(sequence + b"\n")
                index.write(
                    f"{template.name}\t{len(sequence)}\t{offset}\t"
                    f"{max(1, len(sequence))}\t{max(1, len(sequence) + 1)}\n"
                )
        os.replace(temporary, destination)
        os.replace(index_temporary, Path(str(destination) + ".fai"))
    finally:
        for leftover in (temporary, index_temporary):
            try:
                leftover.unlink()
            except FileNotFoundError:
                pass


def read_templates(path: str, load_sequences: bool = False) -> List[LocalTemplate]:
    output: List[LocalTemplate] = []
    for fields, sequence in iter_fasta(path):
        attributes = _header_attributes(fields)
        protocol = attributes.get("protocol")
        if protocol not in {PROTOCOL, LEGACY_PROTOCOL}:
            raise ValueError(f"{path}: unsupported local-template protocol")
        coordinate = LOCAL_COORDINATE_RE.fullmatch(fields[1])
        if coordinate is None:
            raise ValueError(f"{path}: malformed template coordinate {fields[1]!r}")
        output_contig, start_text, end_text, strand = coordinate.groups()
        if strand != "+":
            raise ValueError(f"{path}: stored template sequence must be forward")
        stored_start, stored_end = int(start_text), int(end_text)
        graph_path = attributes["graph_path"]
        if protocol == PROTOCOL:
            output_contig = alternative_locus_name(output_contig)
            start, end = stored_start, stored_end
            source_start = int(attributes["source_start"])
            source_end = int(attributes["source_end"])
            storage_strand = "+"
            expected_contig = alternative_locus_name(graph_path)
            if (
                fields[0] != expected_contig
                or output_contig != expected_contig
                or start != 0
            ):
                raise ValueError(
                    f"{path}: v2 template {fields[0]!r} does not use its "
                    "alternative-locus FASTA ID and local coordinates"
                )
        else:
            # v1 used the original assembly contig and offsets as the public
            # coordinate system.  Translate it at read time so old checkpoint
            # catalogs immediately acquire the v2 alternative-locus contract.
            output_contig = alternative_locus_name(graph_path)
            start, end = 0, stored_end - stored_start
            source_start, source_end = stored_start, stored_end
            storage_strand = attributes["source_strand"]
        if source_end <= source_start or end <= start:
            raise ValueError(f"{path}: invalid template coordinate {fields[1]!r}")
        if source_end - source_start != end - start:
            raise ValueError(
                f"{path}: source provenance span disagrees with alternative "
                f"locus {fields[0]!r}"
            )
        if len(sequence) != end - start:
            raise ValueError(
                f"{path}: template {fields[0]!r} sequence length does not "
                "match its local coordinate span"
            )
        exposed_sequence = sequence
        if load_sequences and storage_strand == "-":
            exposed_sequence = reverse_complement(sequence)
        output.append(LocalTemplate(
            name=fields[0],
            graph_name=attributes["graph"],
            graph_prefix=attributes["graph_prefix"],
            graph_path=graph_path,
            source_haplotype=attributes["source_haplotype"],
            source_contig=attributes["source_contig"],
            output_contig=output_contig,
            start=start,
            end=end,
            source_strand=attributes["source_strand"],
            source_start=source_start,
            source_end=source_end,
            storage_strand=storage_strand,
            sequence=exposed_sequence if load_sequences else "",
            backbone=attributes.get("backbone", ""),
            lift_status=attributes.get("lift_status", ""),
            reference=attributes.get("reference", "."),
            lift_note=attributes.get("lift_note", "."),
            fixed=attributes.get('sequence_role') == 'imported_alternative',
        ))
    return output


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-folder", required=True)
    parser.add_argument("--graph-list", required=True)
    parser.add_argument("--align-ref", required=True)
    parser.add_argument("--reference-haplotype", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--disabled",
        action="store_true",
        help="write an empty checkpoint catalog and disable fallback",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    # The reference name is intentionally explicit in the CLI and retained for
    # auditability even though absence is determined from its extracted align.
    if not args.reference_haplotype.strip():
        raise ValueError("--reference-haplotype is empty")
    templates = (
        [] if args.disabled else discover_templates(
            args.graph_folder, args.graph_list, args.align_ref, args.reference_haplotype,
        )
    )
    write_templates(args.output, templates)
    Path(str(args.output)+'.selection.complete').write_text(
        f'local-template-selection-v4-fixed-alternatives:{args.reference_haplotype}:{args.align_ref}\n')
    print(
        f"[local_reference_templates] wrote {len(templates)} template "
        f"interval(s) to {args.output}",
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"local_reference_templates.py: error: {error}", file=__import__("sys").stderr)
        raise SystemExit(1)
