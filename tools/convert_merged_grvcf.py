#!/usr/bin/env python3
"""Convert merged grVCFs: split them into individual grVCFs, or exact -> CIGAR.

  split     merged grVCFs (exact or CIGAR version) -> one grVCF per sample
  to-cigar  merged exact grVCFs -> merged CIGAR grVCFs (split, then merge
            again without realignment with tools/merge_grvcfs.py)

Give the merged files of one cohort together (cohort.sv.vcf, cohort.indel.vcf
and cohort.snp.vcf): a sample's insertion is rebuilt from its row, the nested
rows it carries (CHROM = parent row ID) and its insertion SNPs.

How each sample's own call is rebuilt:
  CIGAR version: the sample keeps its own breakpoint (row POS + TEMPLATEOFFSET),
      its own deletion size, and its insertion bases (the representative
      edited by its nested rows and insertion SNPs).
  exact version: the sample sits on the representative's breakpoint; the
      bases it deletes beyond it are <row>_F<n> rows and a substitution there
      is a SNP row. A row, its _F rows and the sample's SNPs between them are
      joined back into one call (separated <row>_S<n> rows stay their own).
The version is read from the ##graphvcfmergeVersion header line; files written
before that line existed need --input-version.

Each sample's ##referenceCoverage and ##pseudoLinearMapping lines come from
the merge's cohort.samples.headers.gz next to the merged files. Without it
the individual files have no coverage: re-merged, those samples are '.'
instead of '0' outside their own calls, and --exact cannot place them.

Examples:
  python tools/convert_merged_grvcf.py split -v merged/cohort.{sv,indel,snp}.vcf \\
      -r chm13.fa -o individual
  python tools/convert_merged_grvcf.py to-cigar -v merged/cohort.{sv,indel,snp}.vcf \\
      -r chm13.fa -O merged_cigar -t 16
"""
import argparse
import bisect
from collections import Counter, defaultdict
import gzip
from pathlib import Path
import re
import shutil
import sys
import tempfile

TOOLS = Path(__file__).resolve().parent
# Pipeline scripts: ROOT/scripts (repository layout) or ROOT (flat copy).
SCRIPTS = (TOOLS.parent / 'scripts' if (TOOLS.parent / 'scripts' / 'cohort_vcf_merge.py').is_file()
           else TOOLS.parent)
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(SCRIPTS / 'tools'))

from allele_edges import edge_owned  # noqa: E402
from check_vcf_lossless import FastaIndex, _open, orient, parse_coord, vcf_unescape  # noqa: E402

VERSION_KEY = '##graphvcfmergeVersion='
SAMPLE_HEADERS = 'cohort.samples.headers.gz'
SV_FORMAT = ('GT:TYPE:SIZE:EXTENDGRAPHCIGAR:ASSEMBLYCONTIG:QUERYCOORD:TEMPLATEOFFSET:'
             'LIFTCATEGORY:ALLELENAME:LABEL_H')
SNP_FORMAT = 'GT:TYPE:SIZE:BASE:ASSEMBLYCONTIG:QUERYCOORD:ALLELENAME:LABEL_H'
FLANK = re.compile(r'_F\d+$')
NAMED_PATH = re.compile(r'[<>]([^:<>]+):\d+=')


def warn(message):
    print(f'[convert_merged_grvcf] warning: {message}', file=sys.stderr, flush=True)


def escape(text):
    return (text.replace('@', '@40').replace(':', '@3A').replace(';', '@3B')
            .replace(',', '@2C').replace('\t', '@09').replace('\n', '@0A'))


class Record:
    """One merged row."""
    __slots__ = ('id', 'chrom', 'pos', 'end', 'svtype', 'svlen', 'seq', 'alt', 'path', 'paclass')

    def __init__(self, fields, info):
        self.id, self.chrom, self.pos = fields[2], fields[0], int(fields[1])
        named = NAMED_PATH.fullmatch(vcf_unescape(info.get('EXTENDGRAPHCIGAR', '')))
        self.path = named.group(1) if named else self.id
        self.alt = fields[4].upper()
        self.svtype = info.get('SVTYPE', 'SNP')
        self.end = int(info['END']) if 'END' in info else self.pos
        self.svlen = int(info['SVLEN']) if 'SVLEN' in info else 0
        seq = info.get('SEQ', '.')
        self.seq = '' if seq in ('', '.') else vcf_unescape(seq)
        self.paclass = info.get('PACLASS', '')


class Obs:
    """One observation of a sample on a merged row."""
    __slots__ = ('kind', 'contig', 'q0', 'q1', 'strand', 'name', 'offset', 'size',
                 'lift', 'label')

    def __init__(self, kind, contig, q0, q1, strand, name, offset, size, lift, label):
        (self.kind, self.contig, self.q0, self.q1, self.strand, self.name, self.offset,
         self.size, self.lift, self.label) = (kind, contig, q0, q1, strand, name, offset,
                                              size, lift, label)

    def key(self):
        return (self.kind, self.contig, self.q0, self.q1, self.strand, self.name)


def _int(text, default=0):
    try:
        return int(text)
    except (TypeError, ValueError):
        return default


def observations(keys, text, record):
    """Every observation of one sample column (comma-joined lists)."""
    values = dict(zip(keys, text.split(':')))
    if values.get('GT') != '1':
        return []

    def column(name, default='.'):
        return values.get(name, default).split(',')
    types, contigs, coords = column('TYPE'), column('ASSEMBLYCONTIG'), column('QUERYCOORD')
    offsets, sizes, names = column('TEMPLATEOFFSET', '0'), column('SIZE'), column('ALLELENAME')
    lifts, labels = column('LIFTCATEGORY'), column('LABEL_H')

    def pick(items, number):
        return items[number if len(items) > 1 else 0]
    output = []
    for number, coord in enumerate(coords):
        parsed = parse_coord(coord)
        if parsed is None:
            continue
        low, high, strand = parsed
        kind = pick(types, number)
        if record.svtype == 'SNP' and low == high:
            # A SNP point on '-' is the traversal boundary: the base is before it.
            low = low - 1 if strand == '-' else low
            high = low + 1
        output.append(Obs(kind, vcf_unescape(pick(contigs, number)), low, high, strand,
                          pick(names, number), _int(pick(offsets, number)),
                          abs(_int(pick(sizes, number))), pick(lifts, number),
                          pick(labels, number)))
    return output


def read_header(vcfs):
    """(version, contig lines, INFO/FORMAT lines, sample names) of the merged files."""
    versions, contigs, definitions, samples = set(), {}, {}, None
    for vcf in vcfs:
        with _open(vcf) as handle:
            for line in handle:
                if line.startswith(VERSION_KEY):
                    versions.add(line.strip()[len(VERSION_KEY):])
                elif line.startswith('##contig=<ID='):
                    contigs.setdefault(line.split('ID=', 1)[1].split(',')[0].rstrip('>\n'), line.rstrip('\n'))
                elif line.startswith(('##INFO=<ID=', '##FORMAT=<ID=')):
                    kind, rest = line.split('=<ID=', 1)
                    definitions.setdefault((kind, rest.split(',')[0]), line.rstrip('\n'))
                elif line.startswith('#CHROM'):
                    names = line.rstrip('\n').split('\t')[9:]
                    if samples is not None and names != samples:
                        raise ValueError(f'{vcf}: samples differ from the other merged files')
                    samples = names
                    break
    return versions, contigs, definitions, samples or []


def merged_version(versions, requested):
    if len(versions) > 1:
        raise ValueError(f'merged files of different versions: {sorted(versions)}')
    found = next(iter(versions), None)
    if requested and found and requested != found:
        raise ValueError(f'the files say {found}, --input-version says {requested}')
    if not (requested or found):
        raise ValueError('these merged files have no ##graphvcfmergeVersion line; '
                         'give --input-version exact or cigar')
    return requested or found


def read_sample_headers(path, wanted):
    lines = defaultdict(list)
    with gzip.open(path, 'rt') as handle:
        for raw in handle:
            sample, _tab, line = raw.rstrip('\n').partition('\t')
            if sample in wanted:
                lines[sample].append(line)
    return lines


class Sample:
    """One sample's rows during a split."""

    def __init__(self, name):
        self.name = name
        self.top = []                              # (record, obs) on reference contigs
        self.children = defaultdict(list)          # parent path -> [(record, obs)]
        self.parents = defaultdict(list)           # path -> [(record, obs)]
        self.snps = []                             # (record, obs) top-level SNPs
        self.windows = []                          # (w0, w1, family) of joined calls
        self.absorbed = defaultdict(list)          # family -> [(record, obs)]


def collect(vcfs, reference, wanted, samples):
    """Rows the wanted samples carry, from all merged files."""
    for vcf in vcfs:
        columns = {}
        with _open(vcf) as handle:
            for line in handle:
                if line.startswith('##'):
                    continue
                if line.startswith('#CHROM'):
                    header = line.rstrip('\n').split('\t')
                    columns = {index: samples[name] for index, name in enumerate(header)
                               if index >= 9 and name in wanted}
                    continue
                fields = line.rstrip('\r\n').split('\t')
                picked = [(index, fields[index]) for index in columns
                          if index < len(fields) and fields[index].startswith('1')]
                if not picked:
                    continue
                info = dict(item.split('=', 1) if '=' in item else (item, '')
                            for item in fields[7].split(';'))
                record = Record(fields, info)
                keys = fields[8].split(':')
                for index, text in picked:
                    sample = columns[index]
                    for obs in observations(keys, text, record):
                        item = (record, obs)
                        sample.parents[record.path].append(item)
                        if record.chrom not in reference:
                            sample.children[record.chrom].append(item)
                        elif record.svtype == 'SNP':
                            sample.snps.append(item)
                        else:
                            sample.top.append(item)


class Rebuild:
    """A sample's alleles from its merged rows (cf. scripts/tools/check_merge_lossless.py,
    with each observation's own breakpoint and size in the CIGAR version)."""

    def __init__(self, sample, exact, stats):
        self.sample, self.exact, self.stats = sample, exact, stats
        self.used = set()
        self.owner = {}
        for path, items in sample.children.items():
            for child, child_obs in items:
                candidates = [(parent, obs) for parent, obs in sample.parents.get(path, ())
                              if self.inside(child_obs, obs)]
                if len(candidates) < 2:
                    continue
                named = [(parent.id, obs) for parent, obs in candidates
                         if child_obs.name not in ('', '.') and obs.name == child_obs.name]
                if len(named) == 1:
                    self.owner[(child.id, child_obs.key())] = (named[0][0], named[0][1].key())

    @staticmethod
    def inside(obs, parent):
        return obs.contig == parent.contig and parent.q0 <= obs.q0 and obs.q1 <= parent.q1

    def deleted(self, record, obs):
        """Reference bases a deletion observation removes."""
        return abs(record.svlen) if self.exact or not obs.size else obs.size

    def allele(self, record, obs, depth=0):
        """Parent sequence edited by the sample's nested rows, recursively."""
        sequence = record.seq
        rows, starts, ends, kept_spans = [], [], [], []
        current = len(sequence)
        sample = self.sample

        def extent(child, child_obs):
            first = child.pos - 1 + (0 if self.exact else child_obs.offset)
            if child.svtype == 'DEL':
                return first, first + self.deleted(child, child_obs)
            if child.svtype == 'SNP':
                return first, first + 1
            return first, (first + child.end - child.pos if child.end > child.pos else first)
        for child, child_obs in sample.children.get(record.path, ()):
            key = (child.id, child_obs.key())
            if not self.inside(child_obs, obs) or key in self.used:
                continue
            if self.owner.get(key, (record.id, obs.key())) != (record.id, obs.key()):
                continue
            first, last = extent(child, child_obs)
            if (child_obs.q0 == child_obs.q1 and obs.q0 < obs.q1 and record.svtype != 'DEL'
                    and child_obs.q0 in (obs.q0, obs.q1)):
                edge = starts if child_obs.q0 == (obs.q1 if obs.strand == '-' else obs.q0) else ends
                edge.append((first, last, (child, child_obs)))
            else:
                rows.append((child, child_obs))
                kept_spans.append((first, last))
                current += (child_obs.q1 - child_obs.q0) - (last - first)
        if record.svtype != 'DEL':
            rows.extend(edge_owned(starts, ends, kept_spans, current, obs.q1 - obs.q0,
                                   len(sequence)))
        edits = []
        for child, child_obs in rows:
            self.used.add((child.id, child_obs.key()))
            first, last = extent(child, child_obs)
            order = child_obs.q0 if child_obs.strand != '-' else -child_obs.q1
            if child_obs.kind == 'INS_SNP' or child.svtype == 'SNP':
                edits.append((first, first + 1, child.alt, order))
            else:
                edits.append((first, last, self.allele(child, child_obs, depth + 1), order))
        limit = self.deleted(record, obs) if record.svtype == 'DEL' else len(sequence)
        pieces, cursor = [], 0
        for start, end, bases, _order in sorted(edits, key=lambda e: (e[0], e[1], e[3])):
            if start < cursor or end > limit:
                self.stats['nested_overlap_or_out_of_bounds'] += 1
                continue
            pieces.append(sequence[cursor:start])
            pieces.append(bases)
            cursor = end
        pieces.append(sequence[cursor:])
        return ''.join(pieces)

    def event(self, record, obs):
        """(r0, r1, bases in reference orientation) of a top-level observation."""
        start = record.pos + (0 if self.exact else obs.offset)
        bases = self.allele(record, obs)
        if record.svtype == 'DEL':
            return start, start + self.deleted(record, obs), bases
        return start, start + (record.end - record.pos if record.end > record.pos else 0), bases


def trim(r0, r1, bases, q0, q1, strand, reference, chrom):
    """Drop reference bases a replacement keeps at either end."""
    if r1 <= r0 or not bases:
        return r0, r1, bases, q0, q1
    ref = reference.fetch(chrom, r0, r1).upper()
    left = 0
    while left < min(len(ref), len(bases)) and ref[left] == bases[left].upper():
        left += 1
    right = 0
    while (right < min(len(ref), len(bases)) - left
           and ref[len(ref) - 1 - right] == bases[len(bases) - 1 - right].upper()):
        right += 1
    bases = bases[left:len(bases) - right]
    if strand == '-':
        q1, q0 = q1 - left, q0 + right
    else:
        q0, q1 = q0 + left, q1 - right
    return r0 + left, r1 - right, bases, q0, q1


def closest_shift(r0, bases, q0, q1, strand, representative, reference, chrom):
    """An insertion slides through reference bases equal to its own ends. The
    exact version keeps no sample breakpoint, so take the placement whose
    sequence best matches the row's representative (fewest edits on re-merge)."""
    def score(sequence):
        prefix = 0
        while prefix < min(len(sequence), len(representative)) and \
                sequence[prefix] == representative[prefix]:
            prefix += 1
        suffix = 0
        while suffix < min(len(sequence), len(representative)) - prefix and \
                sequence[-1 - suffix] == representative[-1 - suffix]:
            suffix += 1
        return prefix + suffix
    representative = representative.upper()
    best = (score(bases.upper()), 0, r0, bases)
    length = reference.length(chrom)
    for step in (1, -1):
        position, sequence, shift = r0, bases, 0
        while abs(shift) < len(bases):
            if step == 1:
                if position >= length or reference.fetch(chrom, position, position + 1).upper() \
                        != sequence[0].upper():
                    break
                sequence, position = sequence[1:] + sequence[0], position + 1
            else:
                if position <= 0 or reference.fetch(chrom, position - 1, position).upper() \
                        != sequence[-1].upper():
                    break
                sequence, position = sequence[-1] + sequence[:-1], position - 1
            shift += step
            candidate = (score(sequence.upper()), -abs(shift), position, sequence)
            if candidate[:2] > best[:2]:
                best = candidate
    shift = best[2] - r0
    move = -shift if strand == '-' else shift
    return best[2], best[3], q0 + move, q1 + move


def join_family(pieces, reference, chrom, stats):
    """One call from a row and its _F pieces (and SNPs between them), if their
    query coordinates line up; else None."""
    pieces = sorted(pieces, key=lambda piece: (piece[0], piece[1]))
    strand, contig = pieces[0][4].strand, pieces[0][4].contig
    w0, w1 = pieces[0][0], max(piece[1] for piece in pieces)
    hap, cursor = [], w0
    first = pieces[0][4]
    origin = first.q1 if strand == '-' else first.q0
    length = 0
    for r0, r1, bases, _record, obs in pieces:
        if r0 < cursor or obs.contig != contig or obs.strand != strand:
            stats['joins_refused'] += 1
            return None
        gap = reference.fetch(chrom, cursor, r0)
        hap.append(gap)
        length += len(gap)
        expected = (origin - length, origin - length - len(bases)) if strand == '-' else \
            (origin + length, origin + length + len(bases))
        actual = (obs.q1, obs.q0) if strand == '-' else (obs.q0, obs.q1)
        if actual != expected:
            stats['joins_refused'] += 1
            return None
        hap.append(bases)
        length += len(bases)
        cursor = r1
    hap = ''.join(hap)
    q0, q1 = (origin - len(hap), origin) if strand == '-' else (origin, origin + len(hap))
    stats['joined_calls'] += 1
    return w0, w1, hap, q0, q1


def sample_calls(sample, exact, reference, stats):
    """The sample's SV calls: (chrom, r0, r1, bases, contig, q0, q1, strand, obs, record)."""
    rebuild = Rebuild(sample, exact, stats)
    families = defaultdict(list)
    for record, obs in sample.top:
        r0, r1, bases = rebuild.event(record, obs)
        family = FLANK.sub('', record.id) if exact else f'{record.id}\t{id(obs)}'
        families[(record.chrom, family)].append((r0, r1, bases, record, obs))
    for (chrom, family), snps in sample.absorbed.items():
        for record, obs in snps:
            families[(chrom, family)].append((record.pos - 1, record.pos, record.alt, record, obs))
    calls = []
    for (chrom, _family), pieces in families.items():
        main = next((piece for piece in pieces if piece[3].svtype != 'SNP'
                     and not FLANK.search(piece[3].id)),
                    next((piece for piece in pieces if piece[3].svtype != 'SNP'), pieces[0]))
        record, obs = main[3], main[4]
        joined = join_family(pieces, reference, chrom, stats) if len(pieces) > 1 else None
        if joined is None:
            if len(pieces) > 1:
                for r0, r1, bases, piece_record, piece_obs in pieces:
                    if piece_record.svtype == 'SNP':
                        sample.snps.append((piece_record, piece_obs))
                    else:
                        calls.append((chrom, r0, r1, bases, piece_obs.contig, piece_obs.q0,
                                      piece_obs.q1, piece_obs.strand, piece_obs, piece_record))
                continue
            r0, r1, bases = main[0], main[1], main[2]
            q0, q1 = obs.q0, obs.q1
        else:
            r0, r1, bases, q0, q1 = joined
        r0, r1, bases, q0, q1 = trim(r0, r1, bases, q0, q1, obs.strand, reference, chrom)
        if exact and r1 == r0 and bases and record.svtype == 'INS' and record.seq:
            moved = closest_shift(r0, bases, q0, q1, obs.strand, record.seq, reference, chrom)
            if moved[0] != r0:
                stats['insertions_reslid'] += 1
            r0, bases, q0, q1 = moved
            r1 = r0
        calls.append((chrom, r0, r1, bases, obs.contig, q0, q1, obs.strand, obs, record))
    for (chrom_parent), items in sample.children.items():
        for child, child_obs in items:
            if (child.id, child_obs.key()) not in rebuild.used:
                stats['nested_without_carried_parent'] += 1
    return calls


def sv_rows(call, reference, used):
    """Individual grVCF row(s): one call, or a deletion and an insertion when it
    replaces bases with others. Rows keep their merged row's ID (unique within
    the file), so merging them again keeps the cohort's row names."""
    chrom, r0, r1, bases, contig, q0, q1, strand, obs, record = call
    parts = []
    if r1 > r0 and bases:
        point = q0 if strand == '+' else q1
        parts.append(('DEL', r0, r1, '', point, point))
        parts.append(('INS', r1, r1, bases, q0, q1))
    elif r1 > r0:
        parts.append(('DEL', r0, r1, '', q0, q0))
    elif bases:
        parts.append(('INS', r0, r0, bases, q0, q1))
    rows = []
    for svtype, start, end, seq, a, b in parts:
        size = end - start if svtype == 'DEL' else len(seq)
        anchor = reference.fetch(chrom, start - 1, start).upper() if start > 0 else 'N'
        cigar = f'>{size}{svtype[0]}'
        coord = f'{a}{strand}' if a == b else f'{a}-{b}{strand}'
        info = [f'SVTYPE={svtype}', f'END={end}', f'SVLEN={-size if svtype == "DEL" else size}',
                f'MAXSIZE={size}', f'SEQ={seq or "."}', f'EXTENDGRAPHCIGAR={cigar}']
        if record.paclass == 'fulllocusdup':
            info.append('PACLASS=fulllocusdup')
        row_id, number = record.id, 1
        while row_id in used:
            number += 1
            row_id = f'{record.id}_{number}'
        used.add(row_id)
        fields = ['1', svtype, str(size), cigar, escape(contig), coord, '0',
                  obs.lift, obs.name, obs.label]
        rows.append((chrom, start, '\t'.join([
            chrom, str(start), row_id, anchor,
            f'<{svtype}>', '.', 'PASS', ';'.join(info), SV_FORMAT, ':'.join(fields)])))
    return rows


def snp_row(record, obs):
    point = obs.q0 + 1 if obs.strand == '-' else obs.q0
    fields = ['1', 'SNP', '1', record.alt, escape(obs.contig), f'{point}{obs.strand}',
              obs.name, obs.label]
    return (record.chrom, record.pos, '\t'.join([
        record.chrom, str(record.pos), f'S_{record.chrom}_{record.pos}_{record.alt}',
        reference_base(record), record.alt, '.', 'PASS', 'NSUP=1', SNP_FORMAT, ':'.join(fields)]))


REFERENCE = None


def reference_base(record):
    return REFERENCE.fetch(record.chrom, record.pos - 1, record.pos).upper()


def assign_windows(sample, exact, reference, stats):
    """Exact version: windows of multi-piece families, so the sample's SNPs
    inside them join the call."""
    if not exact:
        return
    spans = defaultdict(list)
    for record, obs in sample.top:
        if record.svtype == 'DEL':
            r0, r1 = record.pos, record.pos + abs(record.svlen)
        else:
            r0 = record.pos
            r1 = record.pos + (record.end - record.pos if record.end > record.pos else 0)
        spans[(record.chrom, FLANK.sub('', record.id))].append((r0, r1, record.id))
    windows = defaultdict(list)
    for (chrom, family), items in spans.items():
        if len(items) > 1:
            windows[chrom].append((min(a for a, _b, _i in items), max(b for _a, b, _i in items),
                                   family))
    kept = []
    for record, obs in sample.snps:
        here = sorted(windows.get(record.chrom, ()))
        index = bisect.bisect_right(here, (record.pos - 1, float('inf'), '')) - 1
        if index >= 0 and here[index][0] <= record.pos - 1 and record.pos <= here[index][1]:
            sample.absorbed[(record.chrom, here[index][2])].append((record, obs))
            stats['snps_joined'] += 1
        else:
            kept.append((record, obs))
    sample.snps = kept


def split(vcfs, reference_path, output, *, sample_names=None, input_version=None,
          sample_headers=None):
    global REFERENCE
    vcfs = [str(Path(path)) for path in vcfs]
    versions, contigs, definitions, samples = read_header(vcfs)
    version = merged_version(versions, input_version)
    exact = version == 'exact'
    wanted = list(sample_names) if sample_names else samples
    missing = [name for name in wanted if name not in samples]
    if missing:
        raise ValueError(f'samples not in the merged files: {missing[:5]}')
    reference = REFERENCE = FastaIndex([reference_path])
    headers_path = Path(sample_headers) if sample_headers else Path(vcfs[0]).parent / SAMPLE_HEADERS
    if headers_path.is_file():
        header_lines = read_sample_headers(headers_path, set(wanted))
    else:
        header_lines = {}
        warn(f'{headers_path} not found: the individual grVCFs get no ##referenceCoverage or '
             '##pseudoLinearMapping lines. Re-merged, these samples are "." instead of "0" '
             'outside their own calls, and --exact cannot place them.')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    state = {name: Sample(name) for name in wanted}
    collect(vcfs, reference, set(wanted), state)
    order = {name: index for index, name in enumerate(contigs)}
    stats = Counter()
    written = []
    for name in wanted:
        sample = state[name]
        assign_windows(sample, exact, reference, stats)
        used = set()
        rows = []
        for call in sample_calls(sample, exact, reference, stats):
            rows.extend(sv_rows(call, reference, used))
        rows.extend(snp_row(record, obs) for record, obs in sample.snps)
        rows.sort(key=lambda row: (order.get(row[0], len(order)), row[0], row[1]))
        path = output / f'{name}.vcf'
        with open(path, 'w') as handle:
            handle.write('##fileformat=VCFv4.2\n')
            for chrom, line in contigs.items():
                if chrom in reference:
                    handle.write(line + '\n')
            for line in definitions.values():
                handle.write(line + '\n')
            for line in header_lines.get(name, ()):
                handle.write(line + '\n')
            handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + name + '\n')
            for _chrom, _pos, line in rows:
                handle.write(line + '\n')
        written.append(path)
    (output / 'vcfs.list').write_text(''.join(f'{path.resolve()}\n' for path in written))
    print(f'[convert_merged_grvcf] {version} version -> {len(written)} individual grVCF(s) in '
          f'{output}; ' + ', '.join(f'{key}={value}' for key, value in sorted(stats.items())),
          file=sys.stderr, flush=True)
    return written


def to_cigar(vcfs, reference_path, output, *, threads=1, input_version=None,
             sample_headers=None, keep_individual=False, merge_options=()):
    versions, _contigs, _definitions, _samples = read_header(vcfs)
    if merged_version(versions, input_version) == 'cigar':
        warn('the merged files are already the CIGAR version; merging them again anyway')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    individual = output / 'individual' if keep_individual else Path(
        tempfile.mkdtemp(prefix='individual.', dir=output))
    try:
        split(vcfs, reference_path, individual, input_version=input_version,
              sample_headers=sample_headers)
        import merge_grvcfs
        merge_grvcfs.main(['-I', str(individual / 'vcfs.list'), '-O', str(output),
                           '-t', str(threads), *merge_options])
    finally:
        if not keep_individual:
            shutil.rmtree(individual, ignore_errors=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    actions = parser.add_subparsers(dest='action', required=True)
    for name, text in (('split', 'merged grVCFs -> one grVCF per sample'),
                       ('to-cigar', 'merged exact grVCFs -> merged CIGAR grVCFs')):
        action = actions.add_parser(name, help=text, description=text)
        action.add_argument('-v', '--vcf', nargs='+', required=True, metavar='VCF',
                            help='merged grVCFs of one cohort (sv, indel and snp)')
        action.add_argument('-r', '--reference', required=True, metavar='FASTA',
                            help='indexed reference FASTA of the merge')
        action.add_argument('--input-version', choices=('exact', 'cigar'),
                            help='for files without a ##graphvcfmergeVersion line')
        action.add_argument('--sample-headers', metavar='FILE',
                            help=f'per-sample header lines (default: {SAMPLE_HEADERS} '
                                 'next to the first merged file)')
        if name == 'split':
            action.add_argument('-o', '--output', required=True, metavar='DIR')
            action.add_argument('--samples', nargs='+', metavar='NAME',
                                help='only these samples (default: all)')
        else:
            action.add_argument('-O', '--output', required=True, metavar='DIR')
            action.add_argument('-t', '--threads', type=int, default=1)
            action.add_argument('--keep-individual', action='store_true',
                                help='keep the split grVCFs in OUTPUT/individual')
    args, merge_options = parser.parse_known_args(argv)
    if args.action == 'split':
        if merge_options:
            parser.error(f'unrecognized arguments: {" ".join(merge_options)}')
        split(args.vcf, args.reference, args.output, sample_names=args.samples,
              input_version=args.input_version, sample_headers=args.sample_headers)
    else:
        sys.path.insert(0, str(TOOLS))
        to_cigar(args.vcf, args.reference, args.output, threads=args.threads,
                 input_version=args.input_version, sample_headers=args.sample_headers,
                 keep_individual=args.keep_individual, merge_options=merge_options)
    return 0


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError) as error:
        print(f'[convert_merged_grvcf] ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
