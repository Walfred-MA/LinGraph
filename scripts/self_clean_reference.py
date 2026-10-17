#!/usr/bin/env python3
"""Expose the existing first-genome cleaning stage without cohort discovery."""
import argparse
import logging
from pathlib import Path
import re

import minsetref


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-s', '--sample', required=True)
    parser.add_argument('-t', '--threads', type=int, default=8)
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--keep-work', action='store_true')
    args = parser.parse_args(argv)
    source = Path(args.input).resolve()
    if (not source.is_file() or source.suffix.lower() in {'.gz', '.bgz', '.bgzf'}
            or not Path(str(source) + '.fai').is_file()):
        parser.error('input must be an uncompressed FASTA with an adjacent .fai')
    if not re.fullmatch(r'[A-Za-z0-9.-]+_h[1-9][0-9]*', args.sample):
        parser.error('--sample must have the form Species_h1')
    output = Path(args.output).resolve()
    if output == source or output in source.parents:
        parser.error('output must be a separate directory that does not contain the input')
    # Use the orchestrator's actual defaults and algorithm. No novel-locus or
    # cohort stages are needed for this single-reference entry point.
    options = minsetref.parse_args([
        '-l', str(source), '-o', str(output), '-t', str(args.threads),
        *(['--resume'] if args.resume else []),
        *(['--keep-work'] if args.keep_work else []),
    ])
    minsetref.check_dependencies(options)
    output.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    minsetref.build_main(minsetref.AssemblyEntry(args.sample, str(source)),
                        options, str(output), logging.getLogger('self_clean_reference'))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
