"""Load a chromosome's compact SNP shards and sort its observations in RAM."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

from graphvcfmerge_snp_compact import DTYPE, LABEL_GROUP, remap_records, source_sections

# Only the loading buffers are bounded. All chromosome columns and the complete
# sort remain in RAM; no sorted binary runs or pair-merge rounds are written.
LOAD_RECORDS = 1_000_000


def _load_source(task, columns):
    binary, spans, destination, query_map, allele_map, group_map, grouped = task
    if not spans:
        return 0
    loaded = 0
    with open(binary, 'rb') as handle:
        for offset, count in spans:
            handle.seek(offset)
            while count:
                take = min(count, LOAD_RECORDS)
                data = np.fromfile(handle, dtype=DTYPE, count=take)
                if len(data) != take:
                    raise ValueError(f'truncated SNP shard: {binary}')
                remap_records(data, query_map, allele_map, group_map, grouped, binary)
                for name in DTYPE.names:
                    columns[name][destination:destination + take] = data[name]
                destination += take
                loaded += take
                count -= take
    return loaded


def sorted_locus(manifest, key, processes):
    """Return one sorted 21-byte array plus its remapped dictionaries/coverage."""
    # Polars fixes its pool size when first initialized. Production chromosome
    # jobs are separate Python processes, so set this before importing it.
    os.environ['POLARS_MAX_THREADS'] = str(max(1, processes))
    try:
        import polars as pl
    except ImportError as error:
        raise RuntimeError('In-memory SNP merging requires polars; install with: python -m pip install polars') from error
    chrom = manifest['loci'][key]
    started = time.monotonic()

    def progress(message):
        print(f'[merge:snp] {chrom}: {message} ({time.monotonic() - started:.1f}s)',
              file=sys.stderr, flush=True)

    paths = list(dict.fromkeys(manifest['sources_by_chrom'].get(key, [])))
    progress(f'preparing {len(paths)} source dictionaries for in-memory merge')
    queries, alleles, groups = {}, {}, {}
    coverage = defaultdict(list)
    tasks, total = [], 0
    for path in paths:
        data = json.loads(Path(path).read_text())
        binary, spans = source_sections(data, key)
        query_map = np.array([queries.setdefault(tuple(q), len(queries)) for q in data['queries']], dtype='<u4')
        allele_map = np.array([alleles.setdefault(a, len(alleles)) for a in data['alleles']], dtype='<u4')
        group_map = np.array([groups.setdefault(g, len(groups)) for g in data.get('label_groups', [])], dtype='<u4')
        grouped = np.array([q[5] == LABEL_GROUP for q in data['queries']], dtype=bool)
        tasks.append((str(binary), spans, total, query_map, allele_map, group_map, grouped))
        total += sum(count for _, count in spans)
        for sample, start, end in data['coverage'].get(chrom, []):
            coverage[sample].append((start, end))
    columns = {name: np.empty(total, dtype=DTYPE.fields[name][0]) for name in DTYPE.names}
    progress(f'loading {total:,} observations; compact data {total * DTYPE.itemsize / 2**30:.2f} GiB '
             f'plus dictionaries and sort workspace; up to {processes} loading threads')
    completed = 0
    last_report = time.monotonic()
    with ThreadPoolExecutor(max_workers=max(1, min(processes, len(tasks)))) as pool:
        futures = [pool.submit(_load_source, task, columns) for task in tasks]
        for future in as_completed(futures):
            completed += future.result()
            if time.monotonic() - last_report >= 10:
                progress(f'loaded {completed:,}/{total:,} observations')
                last_report = time.monotonic()
    del tasks, futures
    frame = pl.DataFrame(columns)
    del columns
    threads = pl.thread_pool_size() if processes > 1 else 1
    progress(f'loaded all observations; sorting in RAM with {threads} Polars threads')
    frame = frame.sort(list(DTYPE.names), multithreaded=processes > 1)
    ordered = np.empty(total, dtype=DTYPE)
    for name in DTYPE.names:
        ordered[name] = frame.get_column(name).to_numpy()
    del frame
    progress('in-memory sort complete')
    return ordered, list(queries), list(alleles), list(groups), coverage
