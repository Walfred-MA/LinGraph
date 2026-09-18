"""Bounded, restartable NumPy merge sort of 21-byte SNP observations.

Each leaf is sorted in memory, then pairs are merged with bounded buffers.
Chromosomes/insertion loci are indexed outside the records throughout. A
completion JSON is committed after each binary, allowing interrupted jobs to
reuse completed leaves and rounds without reading the original VCFs again.
"""
from collections import defaultdict
import hashlib
import json
import multiprocessing
import os
from pathlib import Path

import numpy as np

from graphvcfmerge_snp_compact import DTYPE, LABEL_GROUP, PROTOCOL, remap_records, source_sections
from graphvcfmerge_checkpoints import write_json, file_stamp

RUN_RECORDS = 1_000_000
MERGE_BLOCK = 65_536


def _valid(path, count):
    path = Path(path)
    marker = path.with_suffix('.json')
    if not path.is_file() or not marker.is_file():
        return False
    return (path.stat().st_size == count * DTYPE.itemsize and
            json.loads(marker.read_text()) == dict(count=count, stamp=file_stamp(path)))


def _commit(temporary, output, count):
    os.replace(temporary, output)
    write_json(Path(output).with_suffix('.json'), dict(count=count, stamp=file_stamp(output)))


def _sort_leaf(task):
    output, segments, mapping, count = task
    if _valid(output, count):
        return output
    with np.load(mapping) as maps:
        query_map, allele_map = maps['query'], maps['allele']
        group_map, grouped_queries = maps['label_group'], maps['grouped_query']
    ordered = np.empty(count, dtype=DTYPE)
    start = 0
    # All spans of a leaf belong to one worker's source and dictionary.
    with open(segments[0][0], 'rb') as handle:
        for source, offset, size in segments:
            handle.seek(offset)
            data = np.fromfile(handle, dtype=DTYPE, count=size)
            if len(data) != size:
                raise ValueError(f'truncated SNP shard: {source}')
            remap_records(data, query_map, allele_map, group_map, grouped_queries, source)
            ordered[start:start+size] = data
            start += size
    ordered.sort(kind='mergesort', order=DTYPE.names)
    temporary = str(output) + '.tmp'
    ordered.tofile(temporary)
    _commit(temporary, output, count)
    return output


def _merge_pair(task):
    output, left_path, right_path, count = task
    if _valid(output, count):
        return output
    left = np.memmap(left_path, mode='r', dtype=DTYPE)
    right = np.memmap(right_path, mode='r', dtype=DTYPE)
    i = j = 0
    temporary = str(output) + '.tmp'
    with open(temporary, 'wb') as handle:
        while i < len(left) and j < len(right):
            a, b = left[i:i+MERGE_BLOCK], right[j:j+MERGE_BLOCK]
            # Emit a full block from one side and the matching prefix from
            # the other. Structured search compares every field numerically.
            k = int(np.searchsorted(b, a[-1], side='right'))
            if k < len(b):
                b = b[:k]
            else:
                a = a[:int(np.searchsorted(a, b[-1], side='right'))]
            merged = np.empty(len(a) + len(b), dtype=DTYPE)
            merged[np.arange(len(a)) + np.searchsorted(b, a, side='left')] = a
            merged[np.arange(len(b)) + np.searchsorted(a, b, side='right')] = b
            merged.tofile(handle)
            i += len(a)
            j += len(b)
        for array, start in ((left, i), (right, j)):
            while start < len(array):
                array[start:start+MERGE_BLOCK].tofile(handle)
                start += MERGE_BLOCK
    del left, right
    if Path(temporary).stat().st_size != count * DTYPE.itemsize:
        raise ValueError('SNP pair merge changed the observation count')
    _commit(temporary, output, count)
    return output


def _parallel(pool, function, tasks):
    if pool is None:
        for task in tasks:
            function(task)
    else:
        for _ in pool.imap_unordered(function, tasks):
            pass


def sorted_runs(manifest, key, processes, *, run_records=RUN_RECORDS):
    loci = manifest['job_loci'][key]
    locus_indexes = {locus: i for i, locus in enumerate(loci)}
    source_paths = list(dict.fromkeys(path for locus in loci
                                    for path in manifest['sources_by_chrom'].get(locus, [])))
    sources = {path: json.loads(Path(path).read_text()) for path in source_paths}
    stamps = []
    for path, data in sources.items():
        binaries = set()
        for locus in data['chroms']:
            if locus in locus_indexes:
                binary, _ = source_sections(data, locus)
                if binary.exists():
                    binaries.add(binary)
        stamps.append([path, data, [(str(p), file_stamp(p)) for p in sorted(binaries)]])
    signature = hashlib.sha256(json.dumps([PROTOCOL, run_records, loci, stamps], sort_keys=True).encode()).hexdigest()
    directory = Path(manifest['directory']) / 'merge_runs' / key / signature
    directory.mkdir(parents=True, exist_ok=True)
    queries, alleles, label_groups = {}, {}, {}
    coverage = defaultdict(list)
    leaves, current = [], [[] for _ in loci]
    for index, (path, data) in enumerate(sources.items()):
        query_map = np.array([queries.setdefault(tuple(q), len(queries)) for q in data['queries']], dtype='<u4')
        allele_map = np.array([alleles.setdefault(a, len(alleles)) for a in data['alleles']], dtype='<u4')
        group_map = np.array([label_groups.setdefault(g, len(label_groups))
                              for g in data.get('label_groups', [])], dtype='<u4')
        grouped_queries = np.array([q[5] == LABEL_GROUP for q in data['queries']], dtype=bool)
        mapping = directory / f'map-{index}.npz'
        if not mapping.exists():
            temporary = str(mapping) + '.tmp.npz'
            np.savez(temporary, query=query_map, allele=allele_map,
                     label_group=group_map, grouped_query=grouped_queries)
            os.replace(temporary, mapping)
        for locus in data['chroms']:
            if locus not in locus_indexes:
                continue
            locus_index = locus_indexes[locus]
            for sample, start, end in data['coverage'].get(manifest['loci'][locus], []):
                coverage[(locus_index, sample)].append((start, end))
            binary, spans = source_sections(data, locus)
            segments, size = [], 0

            def leaf():
                output = str(directory / f'leaf-{locus_index}-{len(current[locus_index])}.bin')
                leaves.append((output, list(segments), str(mapping), size))
                current[locus_index].append((output, size))

            for offset, count in spans:
                while count:
                    take = min(count, run_records - size)
                    segments.append((str(binary), offset, take))
                    size += take
                    offset += take * DTYPE.itemsize
                    count -= take
                    if size == run_records:
                        leaf()
                        segments, size = [], 0
            if size:
                leaf()
    write_json(directory / 'dictionaries.json', dict(queries=list(queries), alleles=list(alleles),
                                                   label_groups=list(label_groups)))
    pool = multiprocessing.get_context('spawn').Pool(processes) if processes > 1 and len(leaves) > 1 else None
    try:
        _parallel(pool, _sort_leaf, leaves)
        level = 0
        while any(len(runs) > 1 for runs in current):
            tasks, following = [], []
            for locus_index, runs in enumerate(current):
                merged = []
                for start in range(0, len(runs), 2):
                    if start + 1 == len(runs):
                        merged.append(runs[start])
                        continue
                    left, right = runs[start:start+2]
                    count = left[1] + right[1]
                    output = str(directory / f'round-{level}-{locus_index}-{start//2}.bin')
                    tasks.append((output, left[0], right[0], count))
                    merged.append((output, count))
                following.append(merged)
            _parallel(pool, _merge_pair, tasks)
            current = following
            level += 1
    except BaseException:
        if pool is not None:
            pool.terminate()
            pool.join()
        raise
    else:
        if pool is not None:
            pool.close()
            pool.join()
    return ([runs[0][0] if runs else None for runs in current],
            list(queries), list(alleles), list(label_groups), coverage)
