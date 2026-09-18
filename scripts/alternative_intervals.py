#!/usr/bin/env python3
"""Original-template BED intervals and calling restrictions for graph CIGARs.

Alignment context is retained. A compact eighth TSV field carries the projected
calling intervals to the VCF caller, which can reject an entire boundary-crossing
event instead of creating a shorter, artificial indel by clipping its CIGAR.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile

BED_NAME = "alternative_intervals.bed"
TAG = "ALTERNATIVE_INTERVALS:Z:"
PROTOCOL = "original-template-call-intervals-v1"


def partition_prefix(name):
    fields = name.split("_")
    return fields[1] if fields[0] == "merged" and len(fields) > 1 else fields[0]


@dataclass(frozen=True)
class Interval:
    contig: str
    start: int
    end: int
    name: str

    @property
    def prefix(self):
        return partition_prefix(self.name)


def read_bed(path):
    intervals = []
    seen = set()
    with open(path) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith(("#", "track ", "browser ")):
                continue
            fields = raw.split()
            if len(fields) != 4:
                raise ValueError(f"{path}:{number}: expected exactly four BED columns")
            contig, start, end, name = fields
            try:
                start, end = int(start), int(end)
            except ValueError as error:
                raise ValueError(f"{path}:{number}: BED coordinates must be integers") from error
            if start < 0 or end <= start or any(c in name for c in "/:;<>"):
                raise ValueError(f"{path}:{number}: invalid original-template BED interval")
            value = Interval(contig, start, end, name)
            if value not in seen:
                seen.add(value)
                intervals.append(value)
    return intervals


def write_bed(summary, output, reference_haplotype="CHM13_h1", all_templates=False):
    """Stream either the legacy nine-column or current graph summary.

    Only original template rows contribute coordinates. Named alternative/novel
    partitions are included even if their template came from the backbone; other
    templates are included when their source haplotype differs from the backbone.
    """
    from assembly_contigs import accepted_contig
    rows = set()
    with open(summary) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"partition", "type", "source_haplotype", "source_contig",
                    "source_start", "source_end"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"{summary}: missing columns: {sorted(required - set(reader.fieldnames or ())) }")
        for number, row in enumerate(reader, 2):
            if row["type"] != "original":
                continue
            prefix = partition_prefix(row["partition"])
            haplotype = row["source_haplotype"]
            contig = row["source_contig"]
            if not accepted_contig(haplotype, contig):
                continue
            if not (all_templates or prefix.startswith(("alternative", "novel"))
                    or haplotype != reference_haplotype):
                continue
            start, end = int(row["source_start"]), int(row["source_end"])
            if start < 0 or end <= start:
                raise ValueError(f"{summary}:{number}: invalid original interval")
            name = f"{prefix}_{contig}_{start}_{end}"
            rows.add((contig, start, end, name))
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".alternative_intervals.", dir=destination.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            for row in sorted(rows, key=lambda row: (partition_prefix(row[3]), *row)):
                handle.write("\t".join(map(str, row)) + "\n")
        # Do not invalidate downstream calls on an unchanged summary rerun.
        if destination.is_file() and destination.read_bytes() == Path(temporary).read_bytes():
            os.unlink(temporary)
        else:
            os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return len(rows)


def merge_intervals(intervals):
    result = []
    for start, end in sorted(set(intervals)):
        if end <= start:
            continue
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result


def contains(regions, chrom, start, end):
    """Full containment; an insertion breakpoint follows BED [start,end)."""
    return any(c == chrom and lo <= start < hi and start <= end <= hi
               for c, lo, hi in regions)


def read_tag(fields):
    values = [value[len(TAG):] for value in fields[7:] if value.startswith(TAG)]
    if not values:
        return None
    if len(values) != 1:
        raise ValueError("duplicate alternative interval annotations")
    value = json.loads(values[0])
    if value.get("protocol") != PROTOCOL or not isinstance(value.get("regions"), list):
        raise ValueError("invalid alternative interval annotation")
    for region in value["regions"]:
        if (len(region) != 3 or not isinstance(region[0], str)
                or not isinstance(region[1], int) or not isinstance(region[2], int)
                or region[1] < 0 or region[2] <= region[1]):
            raise ValueError("invalid projected alternative calling interval")
    return value


def append_tag(row, value):
    if value is None:
        return row
    fields = row.split("\t")
    fields = fields[:7] + [field for field in fields[7:] if not field.startswith(TAG)]
    fields.append(TAG + json.dumps(value, separators=(",", ":"), sort_keys=True))
    return "\t".join(fields)


def filter_result(result, value):
    if value is None:
        return result
    regions = tuple(tuple(region) for region in value['regions'])
    retained = []
    for event in result.svs:
        if contains(regions, event.ref_contig, event.ref_start, event.ref_end):
            event.alternative_regions = regions
            retained.append(event)
    # Nested events are retained only under an insertion that itself passed.
    parents = [event for event in retained if event.ins_size > 0]
    retained_ids = {id(event) for event in retained}
    for event in result.svs:
        if id(event) in retained_ids or not event.type_prefix.startswith('INS_'):
            continue
        if any(parent.qry_contig == event.qry_contig
               and parent.qry_start <= event.qry_start <= event.qry_end <= parent.qry_end
               for parent in parents):
            retained.append(event)
    result.alternative_filtered_events += len(result.svs) - len(retained)
    result.svs = retained
    snps = [event for event in result.snps
            if contains(regions, event.chrom, event.pos0, event.pos0 + 1)
            or (event.type_prefix.startswith('INS_') and any(
                parent.qry_contig == event.qry_contig
                and parent.qry_start <= event.qry_pos < parent.qry_end
                for parent in parents))]
    result.alternative_filtered_events += len(result.snps) - len(snps)
    result.snps = snps
    # Outside the calling interval there is no reference/no-call evidence.
    for attribute in ('coverage_spans', 'filtered_spans'):
        cropped = []
        for query, chrom, start, end in getattr(result, attribute):
            for target, lo, hi in regions:
                left, right = max(start, lo), min(end, hi)
                if chrom == target and left < right:
                    cropped.append((query, chrom, left, right))
        setattr(result, attribute, cropped)
    return result


class AlternativeIntervals:
    def __init__(self, path, templates=()):
        self.intervals = read_bed(path)
        self.by_prefix = {}
        for interval in self.intervals:
            self.by_prefix.setdefault(interval.prefix, []).append(interval)
        canonical = sorted((v.contig, v.start, v.end, v.name) for v in self.intervals)
        self.digest = hashlib.sha256(json.dumps(canonical).encode()).hexdigest()
        self.template_sources = {
            template.graph_path: (template.source_contig, template.source_start,
                                  template.source_end, template.source_strand)
            for template in templates
        }

    def tag(self, regions):
        grouped = {}
        for chrom, start, end in regions:
            grouped.setdefault(chrom, []).append((start, end))
        return {"protocol": PROTOCOL, "bed_sha256": self.digest,
                "regions": [[chrom, start, end] for chrom in sorted(grouped)
                            for start, end in merge_intervals(grouped[chrom])]}

    def _path_ranges(self, path, length, coordinates):
        import graphcigartoref as core
        name = core.path_key(path)
        intervals = self.by_prefix.get(partition_prefix(name))
        if intervals is None:
            return None
        source = self.template_sources.get(name)
        if source is None:
            coord = coordinates.get(path) or coordinates.get(name)
            if coord is None:
                raise ValueError(f"--alternative: graph path {path!r} lacks source coordinates")
            source = (coord.chrom, coord.start, coord.end, coord.strand)
        contig, start, end, strand = source
        if end - start != length or strand not in {"+", "-"}:
            raise ValueError(f"--alternative: graph path {path!r} source span/orientation is invalid")
        ranges = []
        for interval in intervals:
            if contig != interval.contig and not contig.endswith(":" + interval.contig):
                continue
            lo, hi = max(start, interval.start), min(end, interval.end)
            if hi > lo:
                ranges.append((lo - start, hi - start) if strand == "+"
                              else (end - hi, end - lo))
        return merge_intervals(ranges)

    def project(self, cigar, coordinates):
        """Project original source intervals onto a reference alignment's query axis."""
        import graphcigartoref as core
        # Most cohort rows are ordinary backbone partitions. Avoid parsing
        # their operation payloads merely to discover that no BED applies.
        if not any(partition_prefix(core.path_key(match.group(2))) in self.by_prefix
                   for match in core._SEG_RE.finditer(cigar)):
            return None
        segments = core.parse_graphic_segments(cigar, "alternative reference")
        query = 0
        projected = []
        for segment in segments:
            ranges = self._path_ranges(segment.path, segment.path_len, coordinates) or []
            if segment.direction == "<":
                ranges = [(segment.path_len - hi, segment.path_len - lo) for lo, hi in ranges]
            ref = 0
            for op in segment.ops:
                qspan = core.query_consume(op.op, op.n)
                rspan = core.ref_consume(op.op, op.n)
                if op.op in {"=", "X"}:
                    for lo, hi in ranges:
                        lo, hi = max(ref, lo), min(ref + op.n, hi)
                        if lo < hi:
                            projected.append((query + lo - ref, query + hi - ref))
                elif op.op == "I" and any(lo <= ref < hi for lo, hi in ranges):
                    projected.append((query, query + op.n))
                query += qspan
                ref += rspan
        return merge_intervals(projected)

    def annotate_pairs(self, pairs, cigars, coordinates, record_coordinates=None):
        """Pair attributes survive all multiprocessing modes, without global state."""
        import graphcigartoref as core
        count = 0
        for pair in pairs:
            reference_cigar = cigars.get(pair.ref_name, "")
            ranges = self.project(reference_cigar, coordinates)
            if ranges is None:
                # A restricted partition without an original-template anchor
                # provides no supported calling interval for this comparison.
                if partition_prefix(pair.query_name or "") not in self.by_prefix:
                    continue
                ranges = []
            coord = (record_coordinates or {}).get(pair.ref_name) or pair.ref_coord
            if coord is None:
                raise ValueError(f"--alternative: {pair.label}: missing reference coordinate")
            span = core.graph_cigar_query_span(reference_cigar)
            if span != coord.end - coord.start:
                raise ValueError(f"--alternative: {pair.label}: reference CIGAR/coordinate span mismatch")
            regions = []
            for lo, hi in ranges:
                start, end = ((coord.start + lo, coord.start + hi) if coord.strand == "+"
                              else (coord.end - hi, coord.end - lo))
                regions.append((coord.chrom, start, end))
            pair.alternative_intervals = self.tag(regions)
            count += 1
        return count

    def direct_tag(self, chrom):
        """Independent local-template coordinates, for explicit deletion rows."""
        source = self.template_sources.get(chrom)
        if source is None:
            if partition_prefix(chrom) in self.by_prefix:
                raise ValueError(f"--alternative: missing source provenance for {chrom!r}")
            return None
        lo, hi = source[1:3]
        ranges = self._path_ranges(chrom, hi - lo, {})
        return None if ranges is None else self.tag((chrom, a, b) for a, b in ranges)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Write original alternative/novel intervals from local_graphs.tsv; no sequences or alignments are read.")
    parser.add_argument("--summary", default="graph/summary/local_graphs.tsv")
    parser.add_argument("-o", "--output", default="")
    parser.add_argument("--reference-haplotype", default="CHM13_h1")
    parser.add_argument("--all-templates", action="store_true", help="include every original row (template-defined backbone)")
    args = parser.parse_args(argv)
    output = args.output or str(Path(args.summary).parent / BED_NAME)
    count = write_bed(args.summary, output, args.reference_haplotype, args.all_templates)
    print(f"[alternative-bed] wrote {count} original intervals to {output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
