#!/usr/bin/env python3
"""Filter mapper BED by score >200, target coverage, and local competition."""
from __future__ import annotations

import argparse
import sys

from map_assemblies_to_block_reference import (
    default_filtered_bed_path,
    filter_competing_mappings_bed,
)


MIN_SCORE_EXCLUSIVE = 200
DEFAULT_MINIMUM_BLOCK_SIZE = 15_000


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("bed", help="query-coordinate-sorted mapper BED")
    parser.add_argument("-o", "--output", help="default: INPUT.filtered.bed")
    parser.add_argument(
        "--minimum-block-size", type=int,
        default=DEFAULT_MINIMUM_BLOCK_SIZE, metavar="INT",
        help=(
            "recursively absorb shorter query blocks into an immediate "
            "neighbor when their gap is no greater than half the smaller "
            "block size (default: 15000; 0 disables)"
        ),
    )
    args = parser.parse_args()
    if args.minimum_block_size < 0:
        parser.error("--minimum-block-size must be nonnegative")
    output = args.output or default_filtered_bed_path(args.bed)
    count = filter_competing_mappings_bed(
        args.bed,
        output,
        minimum_score_exclusive=MIN_SCORE_EXCLUSIVE,
        scale_partial_scores=True,
        minimum_block_size=args.minimum_block_size,
    )
    print(f"Wrote {count} filtered mappings to {output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
