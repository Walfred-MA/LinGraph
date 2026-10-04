"""Indexed, unstranded BED exclusions on assembly/query coordinates."""
from bisect import bisect_left
from collections import defaultdict
import gzip


class QueryErrorBed:
    def __init__(self, intervals):
        self.intervals = {}
        self.ends = {}
        for chrom, spans in intervals.items():
            merged = []
            for start, end in sorted(spans):
                if start == end:
                    continue
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
                else:
                    merged.append((start, end))
            self.intervals[chrom] = tuple(merged)
            self.ends[chrom] = tuple(end for _, end in merged)

    @classmethod
    def read(cls, path):
        intervals = defaultdict(list)
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as handle:
            for line_no, raw in enumerate(handle, 1):
                fields = raw.split()
                if not fields or fields[0].startswith('#') or fields[0] in {'track', 'browser'}:
                    continue
                try:
                    if len(fields) < 3:
                        raise ValueError()
                    start, end = int(fields[1]), int(fields[2])
                    if start < 0 or end < start:
                        raise ValueError()
                except ValueError:
                    raise ValueError(
                        f'{path}:{line_no}: expected BED query-contig, start, end '
                        '(0-based, half-open; 0 <= start <= end)'
                    ) from None
                intervals[fields[0]].append((start, end))
        return cls(intervals)

    def overlaps(self, chrom, start, end):
        if end <= start:
            return False
        spans = self.intervals.get(chrom, ())
        # Skip regions ending at or before the query's start. Adjacent
        # intervals do not share a base under the half-open convention.
        index = bisect_left(self.ends.get(chrom, ()), start + 1)
        return index < len(spans) and spans[index][0] < end
