#!/usr/bin/env python3
"""Compare each haplotype's variants in the merged VCFs with its GAF walks.

Per sample, from the merged VCFs (carrier columns, GT 1):
  ins_bp   inserted bases: INS SVLEN (other rows: INFO/SEQ length)
  del_bp   deleted bases: DEL |SVLEN| (other rows: END - POS)
  snps     SNP rows (only with the SNP VCF)
and from its GAF records walked on the GFA (rGFA SN/SO tags; each record cut
to its path start/end):
  ins_bp   bases the walks cover on variant sequences (a non-path SN, TP not
           snp), plus the bases they jump over there
  skip_bp  bases a walk jumps over on any sequence (a deletion, or the base a
           SNP replaces): between consecutive visits of one sequence, at the
           start of a variant sequence entered past its first base, and at its
           end when left early. Records separated only by an excluded span
           (SAMPLE.gaf.breaks.tsv excluded_row) are read as one walk.
  snps     SNP sequences (TP:Z:snp) entered
del_bp from the GAF is skip_bp - snps. Each sample's values must be equal;
rows left out of the graph (cohort.sv.10Ksingleton.vcf) are not given as
--vcf, and the GAF skips them as well.

Limits: a deletion or partly walked variant at a GAF record break (contig end,
unplaced interval, excluded row) is seen only from the VCF; full-locus
duplication alleles walk a shared template path, which is not a variant
sequence. Such samples are reported as differing.

    python tools/compare_vcf_gaf_stats.py -g cohort.gfa -a gaf/*.gaf \\
        -v cohort.sv.vcf cohort.indel.vcf [cohort.snp.vcf] -o vcf_gaf_stats.tsv
"""
import argparse
from array import array
from collections import defaultdict
import gzip
from pathlib import Path
import re
import sys

STEP = re.compile(r'([<>])(\d+)')


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def read_gfa(path):
    """Per segment: sequence index, offset, length; per sequence: name, length,
    whether it is a path (P line), whether it is a SNP allele."""
    names, index = [], {}
    seq, offset, length = array('l'), array('q'), array('q')
    snp, paths = set(), set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith('S\t'):
                fields = line.rstrip('\n').split('\t')
                tags = dict(field.split(':', 2)[::2] for field in fields[3:] if field.count(':') >= 2)
                name = tags.get('SN')
                if name is None:
                    raise ValueError(f'{path}: S line without SN tag (not a GrvcfGraph rGFA)')
                if name not in index:
                    index[name] = len(names)
                    names.append(name)
                number = int(fields[1])
                while len(seq) <= number:
                    seq.append(-1), offset.append(0), length.append(0)
                seq[number] = index[name]
                offset[number] = int(tags.get('SO', 0))
                length[number] = int(tags['LN']) if 'LN' in tags else len(fields[2])
                if tags.get('TP') == 'snp':
                    snp.add(index[name])
            elif line.startswith('P\t'):
                paths.add(line.split('\t', 2)[1])
    size = [0] * len(names)
    for number in range(len(seq)):
        if seq[number] >= 0:
            size[seq[number]] = max(size[seq[number]], offset[number] + length[number])
    is_path = [name in paths for name in names]
    return dict(names=names, seq=seq, offset=offset, length=length, size=size,
                is_path=is_path, snp=snp)


def excluded_joins(gaf_path):
    """(contig, end, start) of the excluded spans GrvcfGraph cut the walks at
    (SAMPLE.gaf.breaks.tsv, outcome excluded_row): the records on either side
    are read as one walk."""
    joins = set()
    breaks = Path(str(gaf_path) + '.breaks.tsv')
    if breaks.is_file():
        for line in open(breaks):
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 6 and fields[5] == 'excluded_row':
                joins.add((fields[0], int(fields[1]), int(fields[2])))
    return joins


def walk_stats(gfa, gaf_path):
    """ins_bp, skip_bp, snps, entered rows and records of one sample's GAF.

    Inserted bases are what the walks cover on variant sequences plus what
    they jump over there (a nested deletion), so a variant cut over two
    records is counted once; records separated only by an excluded span are
    read as one walk."""
    seq, offset, length, size = gfa['seq'], gfa['offset'], gfa['length'], gfa['size']
    is_path, snp = gfa['is_path'], gfa['snp']
    out = dict(ins_bp=0, skip_bp=0, snps=0, rows=0, records=0, mixed_steps=0)
    joins = excluded_joins(gaf_path)
    state = dict(last_end={}, entered=[], last_seq=None, contig=None, qend=None)

    def skipped(s, bases):
        out['skip_bp'] += bases
        if not is_path[s] and s not in snp:
            out['ins_bp'] += bases

    def finish():
        # Variant sequences left before their end (a deletion at their end).
        for s in state['entered']:
            if s != state['last_seq'] and state['last_end'][s][1] < size[s]:
                skipped(s, size[s] - state['last_end'][s][1])
        state.update(last_end={}, entered=[], last_seq=None)

    with open_text(gaf_path) as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or not fields[5]:
                continue
            steps = [(o, int(n)) for o, n in STEP.findall(fields[5])]
            if not steps:
                continue
            out['records'] += 1
            contig, qs, qe = fields[0], int(fields[2]), int(fields[3])
            path_length, path_start, path_end = int(fields[6]), int(fields[7]), int(fields[8])
            if not (state['contig'] == contig and (contig, state['qend'], qs) in joins):
                finish()
            state['contig'], state['qend'] = contig, qe
            if sum(o == '<' for o, _ in steps) * 2 > len(steps):
                steps = [('>' if o == '<' else '<', n) for o, n in reversed(steps)]
                path_start, path_end = path_length - path_end, path_length - path_start
            first_step = not state['last_end']
            done = 0
            for position, (orientation, number) in enumerate(steps):
                s = seq[number]
                begin, done = done, done + length[number]
                low, high = max(begin, path_start), min(done, path_end)   # walked part of this step
                if orientation == '<':
                    out['mixed_steps'] += 1
                    continue
                if high <= low:
                    continue
                start, end = offset[number] + (low - begin), offset[number] + (high - begin)
                if not is_path[s] and s not in snp:
                    out['ins_bp'] += end - start
                state['last_seq'] = s
                if s in state['last_end']:
                    # Covered so far on s: [lo, hi). Records joined across an
                    # excluded span on the reverse strand come in decreasing
                    # order, so a gap can lie on either side.
                    lo, hi = state['last_end'][s]
                    if start > hi:
                        skipped(s, start - hi)
                    elif end < lo:
                        skipped(s, lo - end)
                    state['last_end'][s] = (min(lo, start), max(hi, end))
                    continue
                state['last_end'][s] = (start, end)
                if is_path[s]:
                    continue
                state['entered'].append(s)
                out['rows'] += 1
                if s in snp:
                    out['snps'] += 1
                if start > 0 and not (first_step and position == 0):
                    skipped(s, start)                   # entered past its first base
            first_step = False
    finish()
    return out


def info_value(info, key):
    match = re.search(r'(?:^|;)' + key + r'=([^;]*)', info)
    return match.group(1) if match else None


def vcf_stats(paths, wanted):
    """Per sample: ins_bp, del_bp, snps, rows from the carrier columns."""
    stats = defaultdict(lambda: dict(ins_bp=0, del_bp=0, snps=0, rows=0))
    has_snp_file = False
    for path in paths:
        names = []
        with open_text(path) as handle:
            for line in handle:
                if line.startswith('##'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if line.startswith('#CHROM'):
                    names = fields[9:]
                    continue
                if len(fields) < 10:
                    continue
                info = fields[7]
                svtype = info_value(info, 'SVTYPE')
                if svtype is None:
                    kind, ins, dele = 'snp', 0, 0
                    has_snp_file = True
                else:
                    size = info_value(info, 'SVLEN') or '0'
                    size = abs(int(size)) if size.lstrip('-').isdigit() else 0
                    if svtype == 'INS':
                        kind, ins, dele = 'sv', size, 0
                    elif svtype == 'DEL':
                        kind, ins, dele = 'sv', 0, size
                    else:
                        seq = info_value(info, 'SEQ') or ''
                        end = info_value(info, 'END')
                        kind = 'sv'
                        ins = len(seq) if re.fullmatch(r'[ACGTNacgtn]+', seq) else 0
                        dele = int(end) - int(fields[1]) if end and end.isdigit() else 0
                for column, field in enumerate(fields[9:]):
                    if field[:1] != '1' or field.split(':', 1)[0] != '1':
                        continue
                    sample = names[column]
                    if wanted is not None and sample not in wanted:
                        continue
                    entry = stats[sample]
                    entry['rows'] += 1
                    if kind == 'snp':
                        entry['snps'] += 1
                    else:
                        entry['ins_bp'] += ins
                        entry['del_bp'] += dele
    return stats, has_snp_file


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-g', '--gfa', required=True, help='rGFA from GrvcfGraph (S lines with SN/SO/LN)')
    parser.add_argument('-a', '--gaf', nargs='+', required=True, help='SAMPLE.gaf files')
    parser.add_argument('-v', '--vcf', nargs='+', required=True,
                        help='the merged VCFs the graph was built from (cohort.sv.vcf, cohort.indel.vcf, '
                             'optionally cohort.snp.vcf)')
    parser.add_argument('-o', '--output', help='per-sample TSV (default: stdout)')
    args = parser.parse_args(argv)

    samples = {Path(path).name[:-len('.gaf')]: path for path in args.gaf if Path(path).name.endswith('.gaf')}
    print(f'[vcf-gaf] reading {args.gfa}', file=sys.stderr, flush=True)
    gfa = read_gfa(args.gfa)
    print(f'[vcf-gaf] reading {len(args.vcf)} VCF(s) for {len(samples)} sample(s)', file=sys.stderr, flush=True)
    vcf, has_snp = vcf_stats(args.vcf, set(samples))
    columns = ['sample', 'vcf_ins_bp', 'gaf_ins_bp', 'vcf_del_bp', 'gaf_del_bp', 'vcf_snps', 'gaf_snps',
               'gaf_records', 'gaf_mixed_orientation_steps', 'match']
    lines, matched = ['\t'.join(columns)], 0
    totals = defaultdict(int)
    for sample in sorted(samples):
        walk = walk_stats(gfa, samples[sample])
        expect = vcf.get(sample, dict(ins_bp=0, del_bp=0, snps=0, rows=0))
        gaf_del = walk['skip_bp'] - walk['snps']
        ok = (expect['ins_bp'] == walk['ins_bp'] and expect['del_bp'] == gaf_del
              and (not has_snp or expect['snps'] == walk['snps']))
        matched += ok
        row = [sample, expect['ins_bp'], walk['ins_bp'], expect['del_bp'], gaf_del,
               expect['snps'] if has_snp else 'NA', walk['snps'], walk['records'], walk['mixed_steps'],
               'yes' if ok else 'no']
        for key, value in zip(columns[1:7], row[1:7]):
            if value != 'NA':
                totals[key] += value
        lines.append('\t'.join(map(str, row)))
    lines.append('\t'.join(map(str, ['TOTAL', *(totals[key] if key in totals else 'NA' for key in columns[1:7]),
                                     '', '', f'{matched}/{len(samples)}'])))
    text = '\n'.join(lines) + '\n'
    if args.output:
        Path(args.output).write_text(text)
    sys.stdout.write(text)
    print(f'[vcf-gaf] {matched} of {len(samples)} sample(s) match', file=sys.stderr)
    return 0 if matched == len(samples) else 1


if __name__ == '__main__':
    raise SystemExit(main())
