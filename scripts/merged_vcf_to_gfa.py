#!/usr/bin/env python3
"""Convert merged variant VCFs to a reference-anchored rGFA graph."""
import argparse
import os
import resource
import sys


def _fasta_values(groups):
    for group in groups or ():
        yield from group


def build_parser():
    parser = argparse.ArgumentParser(
        description='Convert merged SNP/indel/SV VCFs and query FASTA intervals to rGFA/GFA1.')
    parser.add_argument('--gfa-mode', choices=('rgfa', 'query'), default='rgfa',
                        help='rgfa: reference/parent graph anchors (default); '
                             'query: legacy local query-flank GFA, not rGFA')
    parser.add_argument('--nested-pos-base', type=int, choices=(0, 1), default=1,
                        help='POS origin on insertion parents (default: 1, as emitted '
                             'by graphvcfmerge); root POS/END are breakpoint offsets')
    parser.add_argument('--max-node-length', type=int, default=1024,
                        help='maximum rGFA segment length (default: 1024 for VG indexing)')
    parser.add_argument('-v', '--vcf', required=True, action='append', nargs='+',
                        help='merged VCFs or VCF.gz files; may be repeated')
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--all', dest='variant_mode', action='store_const', const='all',
                      help='include SNPs and SVs of every size (default)')
    mode.add_argument('--svonly', dest='variant_mode', action='store_const', const='svonly',
                      help='include only SVs at or above --svcutoff')
    parser.set_defaults(variant_mode='all')
    parser.add_argument('--svcutoff', type=int, default=20, help='SV minimum for --svonly [20]')
    parser.add_argument('--insertion-only', nargs='?', const=50, type=int, metavar='SIZE',
                        help='export only insertions, with minimum size SIZE (default when flag is present: 50)')
    parser.add_argument('-q', '--query-fasta-list',
                        help='NAME FASTA [FAI] per line; NAME must match a VCF sample')
    parser.add_argument('--graph-folder', default='graph', help='graph root (default: graph)')
    parser.add_argument('-r', '--reference-fasta', help='backbone FASTA; default: GRAPH/inputs/reference_alternatives_novels.fa')
    parser.add_argument('--reference-haplotype', help='selected backbone sample/haplotype; inferred from source metadata when available')
    parser.add_argument('--reference-fai', help='optional index for --reference-fasta')
    parser.add_argument('--local-reference-templates', action='append', default=[], metavar='FASTA',
                        help='fallback templates freshly lifted to this backbone; may be repeated')
    parser.add_argument('-a', '--alternatives-fasta', action='append', nargs='+', default=[],
                        metavar='FASTA', help='additional alternative/novel loci to include as complete paths; may be repeated')
    parser.add_argument('--local-path-fasta', action='append', nargs='+', default=[],
                        metavar='FASTA', help='local graph catalog used only to translate VCF targets '
                        'onto included reference/alternative/novel paths; never exported as extra paths')
    parser.add_argument('-o', '--output', help='output GFA or GFA.gz')
    parser.add_argument('--anchor', type=int, default=50,
                        help='query/graph flank on each side (default: 50)')
    parser.add_argument('--anchor-bed', help='sorted BED output (default: OUTPUT.anchors.bed)')
    parser.add_argument('--size-cutoff', type=int, default=0,
                        help='minimum retained variant size (default: 0, include all sizes)')
    parser.add_argument('--unique-size-cutoff', type=int, default=0,
                        help='additional minimum for #unique paths')
    parser.add_argument('-t', '--processes', type=int,
                        default=int(os.environ.get('SLURM_CPUS_PER_TASK', '1')),
                        help='parallel query FASTA readers (default: Slurm CPUs or 1)')
    parser.add_argument('--bed-only', action='store_true',
                        help='write BED/mapping outputs without GFA')
    parser.add_argument('--validate-only', action='store_true',
                        help='resolve the graph without writing GFA')
    return parser


def _log(message):
    maximum = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    mib = maximum / (1024 * 1024 if sys.platform == 'darwin' else 1024)
    print(f'[vcf-to-gfa] {message}; parent peak RSS={mib:.1f} MiB',
          file=sys.stderr, flush=True)


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    args.vcf = list(_fasta_values(args.vcf))
    if args.svcutoff < 1:
        parser.error('--svcutoff must be positive')
    if args.variant_mode == 'svonly':
        args.size_cutoff = max(args.size_cutoff, args.svcutoff)
    if args.insertion_only is not None:
        if args.insertion_only < 0:
            parser.error('--insertion-only size must be nonnegative')
        args.size_cutoff = max(args.size_cutoff, args.insertion_only)
    if args.anchor < 0 or args.size_cutoff < 0 or args.unique_size_cutoff < 0:
        parser.error('anchor and size cutoffs must be nonnegative')
    if not args.query_fasta_list:
        parser.error('-q/--query-fasta-list is required')
    if args.processes < 1:
        parser.error('--processes must be positive')
    if args.max_node_length < 1:
        parser.error('--max-node-length must be positive')
    if args.gfa_mode == 'rgfa' and not args.reference_fasta:
        args.reference_fasta = os.path.join(args.graph_folder, 'inputs', 'reference_alternatives_novels.fa')
    if args.reference_fai and not args.reference_fasta:
        parser.error('--reference-fai requires --reference-fasta')
    if (args.local_path_fasta or args.local_reference_templates) and args.gfa_mode != 'rgfa':
        parser.error('local path/template catalogs require --gfa-mode rgfa')
    if not args.bed_only and not args.validate_only and not args.output:
        parser.error('--output is required unless --bed-only or --validate-only is used')
    paths = [*args.vcf, args.query_fasta_list, args.reference_fasta,
             args.reference_fai, *_fasta_values(args.alternatives_fasta),
             *_fasta_values(args.local_path_fasta), *args.local_reference_templates]
    for path in paths:
        if path and not os.path.isfile(path):
            parser.error(f'input file does not exist: {path}')
    from gfa_interval_pipeline import run
    return run(args, log=_log)


if __name__ == '__main__':
    raise SystemExit(main())
