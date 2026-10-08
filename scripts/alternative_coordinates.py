#!/usr/bin/env python3
"""Report a per-sample VCF's reference-free calls on its run's alternative loci.

graphreftovcf calls a reference-free graph against its local template (the
graph's reference-class path, its original block plus anchor flanks) and names
the contig after the template. alternative_loci.py cuts the templates at the
original blocks into the run's alternative loci. This rewrites one VCF from
template to locus coordinates with alternative_loci.pieces.tsv:

  - rows: CHROM template -> locus, POS/END shifted, IDs renamed (and every
    reference to the old ID in the row); a row not entirely inside a piece
    (a template flank, or across a cut) goes to --outside instead;
  - ##contig and ##alternativeLocus: one line per locus, with the locus's
    original source interval;
  - ##referenceCoverage: clipped to the pieces and shifted;
  - ##pseudoLinearMapping: Reference shifted; a mapping that leaves the
    pieces is clipped, and its Query clipped to the same bases through the
    sample's own calls on the template (a cut never falls inside a call), so
    the clipped region still rebuilds exactly. Path=alt lines outside the
    pieces keep their ID with Reference ".".

Rows and header lines on other contigs pass through unchanged.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import defaultdict
import gzip
import os
from pathlib import Path
import re
import sys

from alternative_loci import read_loci, read_pieces

MARKER = "##alternativeLociCoordinates="
_INTERVAL = re.compile(r'(.+):(\d+)-(\d+)([+-])')
_ID = re.compile(r'(?:^|[<,])ID=([^,>]+)')
_COORD = re.compile(r'(\d+)(?:-(\d+))?([+-]?)')
SV_ALLELES = {'<INS>', '<DEL>', '<SUB>', '<DUP>', '<INV>'}


def _open(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def _meta_value(line, key):
    match = re.search(r'(?:<|,)' + key + r'=("(?:[^"\\]|\\.)*"|[^,>]*)', line)
    if match is None:
        return None
    value = match.group(1)
    return value[1:-1] if value.startswith('"') else value


def _replace_meta(line, key, value, quoted=True):
    text = f'"{value}"' if quoted else str(value)
    return re.sub(r'((?:<|,)' + key + r'=)("(?:[^"\\]|\\.)*"|[^,>]*)',
                  lambda match: match.group(1) + text, line, count=1)


class Pieces:
    """template -> sorted (t0, t1, locus, l0)."""

    def __init__(self, pieces):
        self.pieces = pieces

    def __contains__(self, template):
        return template in self.pieces

    def inside(self, template, start, end):
        """(locus, offset) when [start, end) lies in one piece (a point at
        start counts as inside, at the piece end as outside)."""
        for t0, t1, locus, l0 in self.pieces.get(template, ()):
            if t0 <= start < t1 and end <= t1:
                return locus, l0 - t0
        return None

    def overlaps(self, template, start, end):
        for t0, t1, locus, l0 in self.pieces.get(template, ()):
            lo, hi = max(start, t0), min(end, t1)
            if lo < hi:
                yield lo, hi, locus, l0 - t0


class Event:
    __slots__ = ('r0', 'r1', 'contig', 'q0', 'q1')

    def __init__(self, r0, r1, contig, q0, q1):
        self.r0, self.r1, self.contig, self.q0, self.q1 = r0, r1, contig, q0, q1


def _unescape(text):
    return (text.replace('@0A', '\n').replace('@09', '\t').replace('@2C', ',')
            .replace('@3B', ';').replace('@3A', ':').replace('@40', '@'))


def row_interval(fields):
    """0-based reference interval a row replaces (lossless-checker rules)."""
    pos = int(fields[1])
    if fields[4] in SV_ALLELES:
        info = dict(item.split('=', 1) if '=' in item else (item, '')
                    for item in fields[7].split(';'))
        return pos, int(info.get('END', pos))
    return pos - 1, pos - 1 + max(1, len(fields[3]))


def row_events(fields, sample_column=9):
    """Query-placed observations (GT=1) of one row, for mapping clips."""
    if len(fields) <= sample_column:
        return []
    keys = fields[8].split(':')
    values = fields[sample_column].split(':')
    if len(keys) != len(values):
        return []
    sample = dict(zip(keys, values))
    if sample.get('GT') != '1':
        return []
    r0, r1 = row_interval(fields)
    contigs = sample.get('ASSEMBLYCONTIG', '.').split(',')
    events = []
    for index, coord in enumerate(sample.get('QUERYCOORD', '.').split(',')):
        match = _COORD.fullmatch(coord)
        if match is None:
            continue
        q0 = int(match.group(1))
        q1 = int(match.group(2)) if match.group(2) else q0
        contig = _unescape(contigs[index if len(contigs) > 1 else 0])
        if fields[4] not in SV_ALLELES:
            if match.group(3) == '-':
                q0 -= 1
            q1 = q0 + 1
        events.append(Event(r0, r1, contig, q0, q1))
    return events


def clip_mapping(query, reference, pieces, events):
    """[(query, reference interval on a locus)] for one mapping line."""
    qcontig, qs, qe, qstrand = query
    template, rs, re_, rstrand = reference
    same = qstrand == rstrand
    mine = [event for event in events
            if event.contig == qcontig and qs <= event.q0 and event.q1 <= qe
            and rs <= event.r0 and event.r1 <= re_]
    mine.sort(key=lambda event: (event.r1, event.r0))
    ends = [event.r1 for event in mine]

    def query_at(x):
        # The last call ending at or before x (an insertion at x is after x).
        index = bisect_right(ends, x) - 1
        while index >= 0 and not mine[index].r0 < x:
            index -= 1
        if index < 0:
            return qs + (x - rs) if same else qe - (x - rs)
        event = mine[index]
        return event.q1 + (x - event.r1) if same else event.q0 - (x - event.r1)

    def outside_call(x, left):
        for event in mine:
            if event.r0 < x < event.r1:
                return event.r1 if left else event.r0
        return None

    output = []
    for lo, hi, locus, shift in pieces.overlaps(template, rs, re_):
        if (lo, hi) == (rs, re_):
            output.append(((qcontig, qs, qe, qstrand), (locus, lo + shift, hi + shift, rstrand)))
            continue
        while lo > rs and (moved := outside_call(lo, True)) is not None:
            lo = moved
        while hi < re_ and (moved := outside_call(hi, False)) is not None:
            hi = moved
        if hi <= lo:
            continue
        a, b = query_at(lo), query_at(hi)
        q0, q1 = (a, b) if same else (b, a)
        output.append(((qcontig, q0, q1, qstrand), (locus, lo + shift, hi + shift, rstrand)))
    return output


def _interval_text(interval):
    contig, start, end, strand = interval
    return f"{contig}:{start}-{end}{strand}"


def _parse_interval(text):
    match = _INTERVAL.fullmatch(text or '')
    if match is None:
        return None
    return match.group(1), int(match.group(2)), int(match.group(3)), match.group(4)


def _rename_id(identifier, template, pos, locus, new_pos):
    token = f"_{template}_{pos}_"
    index = identifier.find(token)
    if index <= 0:
        return identifier
    return f"{identifier[:index]}_{locus}_{new_pos}_{identifier[index + len(token):]}"


def _locus_header(locus):
    return ("##alternativeLocus=<"
            f"ID={locus['name']},Length={locus['length']},"
            f"SourceHaplotype=\"{locus['haplotype']}\",SourceContig=\"{locus['contig']}\","
            f"SourceStart={locus['start']},SourceEnd={locus['end']},"
            f"SourceStrand={locus['strand']},SourceCoordinates=ZeroBasedHalfOpen>\n")


def convert(vcf, loci_prefix, output, outside=None):
    loci = read_loci(str(loci_prefix) + '.tsv')
    pieces = Pieces(read_pieces(str(loci_prefix) + '.pieces.tsv'))
    header, template_rows = [], []
    declared, contig_ids = set(), []
    with _open(vcf) as handle:
        for raw in handle:
            if raw.startswith('#'):
                if raw.startswith(MARKER):
                    raise ValueError(f'{vcf}: already in alternative-locus coordinates')
                header.append(raw)
                if raw.startswith('##alternativeLocus=<'):
                    declared.add(_meta_value(raw, 'ID'))
                elif raw.startswith('##contig=<'):
                    contig_ids.append(_meta_value(raw, 'ID'))
                continue
            chrom = raw.split('\t', 1)[0]
            if chrom in declared or chrom in pieces:
                template_rows.append(raw.rstrip('\n').split('\t'))
    templates = {name for name in declared if name}
    templates.update(fields[0] for fields in template_rows)
    unknown = sorted(name for name in templates if name not in pieces)
    if unknown:
        raise ValueError(f'{vcf}: templates absent from {loci_prefix}.pieces.tsv: {unknown[:5]}')
    used_loci = []
    for template in sorted(templates):
        for _t0, _t1, locus, _l0 in pieces.pieces[template]:
            if locus not in used_loci:
                used_loci.append(locus)
    clash = sorted(set(used_loci) & (set(contig_ids) - templates))
    if clash:
        raise ValueError(f'{vcf}: alternative locus names collide with contigs: {clash[:5]}')

    # Rows: inside one piece -> locus; anything else -> outside.
    events = defaultdict(list)
    converted, rejected = [], []
    dropped = defaultdict(list)
    for index, fields in enumerate(template_rows):
        template = fields[0]
        for event in row_events(fields):
            events[template].append(event)
        r0, r1 = row_interval(fields)
        placed = pieces.inside(template, r0, r1)
        if placed is None:
            rejected.append(fields)
            # A call cut by a locus edge leaves those locus bases without the
            # sample's evidence: they are not covered.
            for lo, hi, locus, shift in pieces.overlaps(template, r0, r1):
                dropped[locus].append((lo + shift, hi + shift))
            continue
        locus, shift = placed
        new = list(fields)
        old_pos = int(fields[1])
        new[0], new[1] = locus, str(old_pos + shift)
        new[2] = _rename_id(fields[2], template, old_pos, locus, old_pos + shift)
        if fields[4] in SV_ALLELES:
            new[7] = ';'.join(
                f"END={int(item[4:]) + shift}" if item.startswith('END=') else item
                for item in fields[7].split(';'))
        if new[2] != fields[2]:
            for column in range(7, len(new)):
                new[column] = new[column].replace(fields[2], new[2])
        converted.append((used_loci.index(locus), int(new[1]), index, new))
    converted.sort(key=lambda item: item[:3])
    seen_ids = set()
    for *_key, fields in converted:
        if fields[2] != '.' and fields[2] in seen_ids:
            raise ValueError(f'{vcf}: duplicate converted ID {fields[2]}')
        seen_ids.add(fields[2])

    # Header.
    output_header, coverage, coverage_at = [], defaultdict(list), None
    counts = defaultdict(int)
    locus_lines_at = contig_lines_at = None
    for raw in header:
        if raw.startswith('##alternativeLocus=<') and _meta_value(raw, 'ID') in templates:
            locus_lines_at = len(output_header) if locus_lines_at is None else locus_lines_at
            continue
        if raw.startswith('##contig=<') and _meta_value(raw, 'ID') in templates:
            contig_lines_at = len(output_header) if contig_lines_at is None else contig_lines_at
            continue
        if raw.startswith('##referenceCoverage=<'):
            chrom = _meta_value(raw, 'Chrom')
            if chrom in templates:
                coverage_at = len(output_header) if coverage_at is None else coverage_at
                sample = _meta_value(raw, 'Sample')
                start, end = int(_meta_value(raw, 'Start')), int(_meta_value(raw, 'End'))
                for lo, hi, locus, shift in pieces.overlaps(chrom, start, end):
                    coverage[(sample, locus)].append((lo + shift, hi + shift))
                continue
        if raw.startswith('##referenceCoverageIntervalCount='):
            output_header.append(None)          # filled in below
            continue
        if raw.startswith('##pseudoLinearMapping=<'):
            reference = _parse_interval(_meta_value(raw, 'Reference'))
            if reference is not None and reference[0] in templates:
                output_header.extend(_convert_mapping(raw, reference, pieces, events, counts))
                continue
        if raw.startswith('#CHROM'):
            output_header.append(
                MARKER + '<Description="Reference-free calls on the alternative loci of '
                'alternative_loci.py (original graph blocks, merged when touching), '
                '0-based on each locus; calls in template flanks are excluded">\n')
        output_header.append(raw)
    coverage_lines = []
    for (sample, locus), intervals in sorted(coverage.items(), key=lambda item: (item[0][0], used_loci.index(item[0][1]))):
        merged = []
        for start, end in sorted(intervals):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        for cut_start, cut_end in dropped.get(locus, ()):
            kept = []
            for start, end in merged:
                if cut_end <= start or end <= cut_start:
                    kept.append([start, end])
                    continue
                if start < cut_start:
                    kept.append([start, cut_start])
                if cut_end < end:
                    kept.append([cut_end, end])
            merged = kept
        coverage_lines.extend(
            f'##referenceCoverage=<Sample="{sample}",Chrom="{locus}",Start={start},End={end}>\n'
            for start, end in merged)
    locus_lines = [_locus_header(loci[name]) for name in used_loci]
    contig_lines = [f"##contig=<ID={name},length={loci[name]['length']}>\n" for name in used_loci]
    for position, lines in sorted(
            ((coverage_at, coverage_lines), (contig_lines_at, contig_lines),
             (locus_lines_at, locus_lines)),
            key=lambda item: -1 if item[0] is None else item[0], reverse=True):
        if not lines:
            continue
        if position is None:
            position = next(i for i, line in enumerate(output_header)
                            if line is not None and line.startswith('#CHROM')) - 1
        output_header[position:position] = lines
    total = sum(1 for line in output_header
                if line is not None and line.startswith('##referenceCoverage=<'))
    output_header = [f"##referenceCoverageIntervalCount={total}\n" if line is None else line
                     for line in output_header]

    output = Path(output)
    temporary = output.with_name(output.name + f'.tmp.{os.getpid()}')
    with _open(vcf) as handle, temporary.open('w') as out:
        out.writelines(output_header)
        for raw in handle:
            if raw.startswith('#'):
                continue
            chrom = raw.split('\t', 1)[0]
            if chrom in templates:
                continue
            out.write(raw)
        for *_key, fields in converted:
            out.write('\t'.join(fields) + '\n')
    os.replace(temporary, output)
    outside = Path(outside) if outside else output.with_name(output.name.replace('.vcf', '') + '.outside.vcf')
    if rejected:
        with outside.open('w') as out:
            out.writelines(header)
            for fields in rejected:
                out.write('\t'.join(fields) + '\n')
    elif outside.exists():
        outside.unlink()
    counts.update(template_rows=len(template_rows), rows_on_loci=len(converted),
                  rows_outside=len(rejected), loci=len(used_loci))
    return dict(counts)


def _convert_mapping(raw, reference, pieces, events, counts):
    path_kind = _meta_value(raw, 'Path')
    query = _parse_interval(_meta_value(raw, 'Query'))
    template, start, end, strand = reference
    if start == end:                                     # insertion point
        placed = pieces.inside(template, start, start)
        if placed is None:
            counts['mapping_lines_outside'] += 1
            return []
        locus, shift = placed
        return [_replace_meta(raw, 'Reference', _interval_text((locus, start + shift, end + shift, strand)))]
    if path_kind == 'alt':                               # a sequence source: places no bases
        placed = pieces.inside(template, start, end)
        if placed is None:
            counts['alt_lines_unanchored'] += 1
            return [_replace_meta(raw, 'Reference', '.')]
        locus, shift = placed
        return [_replace_meta(raw, 'Reference', _interval_text((locus, start + shift, end + shift, strand)))]
    if query is None or query[1] == query[2]:            # deletion: query point
        lines = [_replace_meta(raw, 'Reference', _interval_text((locus, lo + shift, hi + shift, strand)))
                 for lo, hi, locus, shift in pieces.overlaps(template, start, end)]
        counts['mapping_lines_outside' if not lines else 'mapping_lines_kept'] += 1
        return lines
    clipped = clip_mapping(query, reference, pieces, events[template])
    if not clipped:
        counts['mapping_lines_outside'] += 1
        return []
    counts['mapping_lines_clipped' if len(clipped) > 1 or clipped[0][1][2] - clipped[0][1][1] != end - start
           else 'mapping_lines_kept'] += 1
    return [_replace_meta(_replace_meta(raw, 'Query', _interval_text(q)), 'Reference', _interval_text(r))
            for q, r in clipped]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--vcf', required=True, help='per-sample VCF on local templates')
    parser.add_argument('-l', '--loci', required=True,
                        help='alternative_loci.py output prefix (PREFIX.tsv, PREFIX.pieces.tsv)')
    parser.add_argument('-o', '--output', required=True, help='VCF on alternative loci')
    parser.add_argument('--outside', help='rows outside the loci (default: OUTPUT.outside.vcf)')
    args = parser.parse_args(argv)
    counts = convert(args.vcf, args.loci, args.output, args.outside)
    print('[alternative-coordinates] ' + ', '.join(f'{key}={value}' for key, value in sorted(counts.items())),
          file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
