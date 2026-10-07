#!/usr/bin/env python3
"""Binary index of the cohort GFA holding only what gfa_sample_gaf.py reads.

Built once per GFA; every gfa_sample_gaf.py job memory-maps it read-only
(``--index DIR``) instead of parsing the GFA text and the variant sidecar, so
jobs and their forked workers share one copy through the page cache.

Kept (raw little-endian arrays, dtypes in index.json):

* ``segment_length``   int64 per segment ID (LN tag, else sequence length)
* ``path_offsets``     int64, CSR offsets of P-line walks (paths + 1)
* ``path_steps``       int64, steps encoded as segment ID * 2 + reverse
* ``step_offset``      uint32, each step's base offset within its path, and
  ``path_length``      uint32, each path's bases: walks search these instead of
                       caching running lengths per path (a path must be < 4 Gb)
* ``name_offsets``     int64, offsets of P-line names in ``names``
* ``names``            uint8, P-line names concatenated
* ``name_hash``        uint64, sorted 64-bit name hashes; ``name_order`` int64
                       gives the path of each (lookups compare the name)
* ``edges``            uint64, sorted canonical link keys
* ``var_path``, ``var_parent``, ``var_start``, ``var_end`` int64 and
  ``var_kind``, ``var_reverse`` int8 per merged-VCF data line, from
  GFA.variants.tsv (paths already resolved to numbers)

Dropped: sequences, other tags, link overlaps, and variant IDs. Semantics
match gfa_sample_gaf.load_gfa / load_variant_index: a repeated segment ID or
P name keeps its last occurrence. ``index.json`` is written last; it records
the size and mtime of the GFA and sidecar it was built from.
"""
import argparse
from array import array
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np

KINDS = {'insertion': 0, 'snp': 1, 'substitution': 2, 'deletion': 3}
FLUSH = 1 << 22  # steps / name bytes buffered before writing


def _log(message):
    print(f'[gfa-gaf-index] {message}', file=sys.stderr, flush=True)


def _hash(name):
    return int.from_bytes(hashlib.blake2b(name, digest_size=8).digest(), 'little')


def _edge_key(a, b):
    """Canonical key of a directed link a->b (as gfa_sample_gaf._edge_key)."""
    ra, rb = b ^ 1, a ^ 1
    if (ra, rb) < (a, b):
        a, b = ra, rb
    return (a << 32) | b


def _stamp(path):
    stat = os.stat(path)
    return {'path': str(Path(path).resolve()), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def _save(folder, name, values, meta):
    values = np.ascontiguousarray(values)
    values.tofile(folder / f'{name}.bin')
    meta['arrays'][name] = [values.dtype.str, int(values.size)]


class PathNames:
    """P-line names by path number, read from the mapped name bytes."""

    def __init__(self, names, offsets):
        self.names, self.offsets = names, offsets

    def __len__(self):
        return len(self.offsets) - 1

    def raw(self, path):
        return self.names[self.offsets[path]:self.offsets[path + 1]].tobytes()

    def __getitem__(self, path):
        if path < 0:
            path += len(self)
        if not 0 <= path < len(self):
            raise IndexError(path)
        return self.raw(path).decode()


class PathIndex:
    """name -> path number (last P line of a repeated name), as a dict."""

    def __init__(self, names, hashes, order):
        self.names, self.hashes, self.order = names, hashes, order

    def get(self, name, default=None):
        raw = name.encode() if isinstance(name, str) else name
        key = np.uint64(_hash(raw))
        low = int(np.searchsorted(self.hashes, key, 'left'))
        high = int(np.searchsorted(self.hashes, key, 'right'))
        found = default
        for position in range(low, high):  # ascending path numbers
            path = int(self.order[position])
            if self.names.raw(path) == raw:
                found = path
        return found

    def __contains__(self, name):
        return self.get(name) is not None

    def __getitem__(self, name):
        path = self.get(name)
        if path is None:
            raise KeyError(name)
        return path


def build(gfa, variants, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / 'index.json').unlink(missing_ok=True)
    meta = {'gfa': _stamp(gfa), 'variants': _stamp(variants), 'arrays': {}}

    segment_ids, segment_lengths = array('q'), array('q')
    edges = array('Q')
    path_offsets, name_offsets = array('q', [0]), array('q', [0])
    name_hashes = array('Q')
    steps, name_bytes = array('q'), bytearray()
    step_count = name_count = 0
    with open(gfa, 'rb') as handle, open(output / 'path_steps.bin', 'wb') as step_out, \
            open(output / 'names.bin', 'wb') as name_out:
        for line in handle:
            kind = line[:1]
            if kind == b'S':
                fields = line.rstrip(b'\n').split(b'\t')
                length = next((int(tag[5:]) for tag in fields[3:] if tag.startswith(b'LN:i:')),
                              len(fields[2]))
                segment_ids.append(int(fields[1]))
                segment_lengths.append(length)
            elif kind == b'L':
                fields = line.split(b'\t', 6)
                a = int(fields[1]) * 2 + (fields[2] == b'-')
                b = int(fields[3]) * 2 + (fields[4] == b'-')
                edges.append(_edge_key(a, b))
            elif kind == b'P':
                fields = line.split(b'\t', 3)
                name = fields[1]
                tokens = fields[2].split(b',')
                for token in tokens:
                    steps.append(int(token[:-1]) * 2 + (token[-1:] == b'-'))
                step_count += len(tokens)
                path_offsets.append(step_count)
                if len(steps) >= FLUSH:
                    steps.tofile(step_out)
                    steps = array('q')
                name_bytes += name
                name_count += len(name)
                name_offsets.append(name_count)
                name_hashes.append(_hash(name))
                if len(name_bytes) >= FLUSH:
                    name_out.write(name_bytes)
                    name_bytes = bytearray()
        steps.tofile(step_out)
        name_out.write(name_bytes)
    meta['arrays']['path_steps'] = ['<i8', step_count]
    meta['arrays']['names'] = ['|u1', name_count]

    ids = np.frombuffer(segment_ids, dtype=np.int64)
    segment_length = np.zeros(int(ids.max(initial=0)) + 1, dtype=np.int64)
    segment_length[ids] = np.frombuffer(segment_lengths, dtype=np.int64)
    segments = int(np.unique(ids).size)
    del segment_ids, segment_lengths, ids
    _save(output, 'segment_length', segment_length, meta)
    paths = len(path_offsets) - 1
    path_offsets = np.frombuffer(path_offsets, dtype=np.int64)
    # Walk pieces store (path, start, end) as uint32 step numbers within a path.
    if paths >= 1 << 32 or (paths and int(np.diff(path_offsets).max()) >= 1 << 32):
        raise SystemExit(f'over {1 << 32} P lines or steps in one P line; uint32 pieces cannot hold it')
    _save(output, 'path_offsets', path_offsets, meta)
    _step_offsets(output, meta, segment_length, path_offsets)
    del path_offsets
    _save(output, 'name_offsets', np.frombuffer(name_offsets, dtype=np.int64), meta)
    hashes = np.frombuffer(name_hashes, dtype=np.uint64)
    order = np.argsort(hashes, kind='stable')
    _save(output, 'name_hash', hashes[order], meta)
    _save(output, 'name_order', order.astype(np.int64), meta)
    del name_hashes, hashes, order
    links = np.unique(np.frombuffer(edges, dtype=np.uint64))
    del edges
    _save(output, 'edges', links, meta)
    _log(f'GFA: {segments} segments, {paths} paths, {len(links)} links')
    del links

    index = PathIndex(*_path_tables(output, meta))
    rows = [array('q') for _ in range(5)] + [array('b'), array('b')]
    count = 0
    missing = 0
    with open(variants) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            vcf_index, _identifier, name, variant_kind, parent_name, low, high, strand = \
                line.rstrip('\n').split('\t')
            i = int(vcf_index)
            count = max(count, i + 1)
            if i < 0 or parent_name == '.':
                continue
            parent = index.get(parent_name)
            path = index.get(name, -1)
            if parent is None or (variant_kind != 'deletion' and path < 0):
                missing += 1
                continue
            for column, value in zip(rows, (i, path, parent, int(low), int(high),
                                            KINDS.get(variant_kind, 0), strand == '-')):
                column.append(value)
    kept = np.frombuffer(rows[0], dtype=np.int64)
    # A repeated VCF index keeps its last row (as load_variant_index does).
    last = len(kept) - 1 - np.unique(kept[::-1], return_index=True)[1]
    kept = kept[last]
    for name, column, dtype, empty in (('var_path', rows[1], np.int64, -1),
                                       ('var_parent', rows[2], np.int64, -1),
                                       ('var_start', rows[3], np.int64, 0),
                                       ('var_end', rows[4], np.int64, 0),
                                       ('var_kind', rows[5], np.int8, -1),
                                       ('var_reverse', rows[6], np.int8, 0)):
        values = np.full(count, empty, dtype=dtype)
        values[kept] = np.frombuffer(column, dtype=dtype)[last]
        if name == 'var_parent':
            exported = int((values >= 0).sum())
        _save(output, name, values, meta)
    if missing:
        _log(f'warning: {missing} indexed variants name a path absent from the GFA')
    _log(f'variant index: {exported} exported of {count} VCF records')
    meta['paths'] = paths
    temporary = output / 'index.json.tmp'
    temporary.write_text(json.dumps(meta, indent=1) + '\n')
    temporary.replace(output / 'index.json')
    _log(f'index written to {output}')


def _step_offsets(output, meta, segment_length, path_offsets):
    """``step_offset``: each step's base offset within its path, and
    ``path_length``: each path's bases (both uint32), in chunks of steps."""
    steps = _array(output, meta, 'path_steps')
    total = len(steps)
    # Exclusive running base count over all steps, at every path start.
    prefix = np.zeros(len(path_offsets), dtype=np.int64)
    carry = 0
    chunk = 4 * FLUSH
    with open(output / 'step_offset.bin', 'wb') as out:
        for low in range(0, total, chunk):
            high = min(total, low + chunk)
            nodes = np.asarray(steps[low:high]) >> 1
            if nodes.size and int(nodes.max()) >= len(segment_length):
                raise SystemExit(f'a P line uses segment {int(nodes.max())}, which has no S line')
            lengths = segment_length[nodes]
            before = np.cumsum(lengths) - lengths + carry
            starting = np.searchsorted(path_offsets, [low, high], side='left')
            prefix[starting[0]:starting[1]] = before[path_offsets[starting[0]:starting[1]] - low]
            owner = np.searchsorted(path_offsets, np.arange(low, high), side='right') - 1
            offset = before - prefix[owner]
            if offset.size and int(offset.max()) >= 1 << 32:
                raise SystemExit(f'a P path exceeds {1 << 32} bases; uint32 offsets cannot hold it')
            offset.astype(np.uint32).tofile(out)
            carry = int(before[-1] + lengths[-1])
    prefix[np.searchsorted(path_offsets, total, side='left'):] = carry
    meta['arrays']['step_offset'] = ['<u4', total]
    length = np.diff(prefix)
    if length.size and int(length.max()) >= 1 << 32:
        raise SystemExit(f'a P path exceeds {1 << 32} bases; uint32 lengths cannot hold it')
    _save(output, 'path_length', length.astype(np.uint32), meta)


def _array(folder, meta, name):
    dtype, size = meta['arrays'][name]
    if not size:
        return np.empty(0, dtype=np.dtype(dtype))
    return np.memmap(Path(folder) / f'{name}.bin', dtype=np.dtype(dtype), mode='r', shape=(size,))


def _path_tables(folder, meta):
    names = PathNames(_array(folder, meta, 'names'), _array(folder, meta, 'name_offsets'))
    return names, _array(folder, meta, 'name_hash'), _array(folder, meta, 'name_order')


def load(folder, gfa=None, variants=None):
    """The gfa_sample_gaf.G entries, memory-mapped from an index folder."""
    folder = Path(folder)
    try:
        meta = json.loads((folder / 'index.json').read_text())
    except FileNotFoundError:
        raise SystemExit(f'{folder}: no complete index (index.json); run gfa_gaf_index.py')
    for key, path in (('gfa', gfa), ('variants', variants)):
        if path is not None:
            stamp = _stamp(path)
            if (stamp['size'], stamp['mtime_ns']) != (meta[key]['size'], meta[key]['mtime_ns']):
                raise SystemExit(f'{folder} was built from another {key} than {path}; rebuild it')
    names, hashes, order = _path_tables(folder, meta)
    tables = {name: _array(folder, meta, name) for name in meta['arrays']
              if name not in ('names', 'name_offsets', 'name_hash', 'name_order')}
    tables['offsets'] = tables.pop('path_offsets')
    tables['steps'] = tables.pop('path_steps')
    tables.update(path_names=names, path_index=PathIndex(names, hashes, order))
    _log(f'index {folder}: {meta["paths"]} paths, {len(tables["edges"])} links, '
         f'{len(tables["var_parent"])} VCF records (memory-mapped)')
    return tables


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-g', '--gfa', required=True, help='rGFA from merged_vcf_to_gfa.py')
    parser.add_argument('--variant-index', help='default: GFA.variants.tsv')
    parser.add_argument('-o', '--output-folder', required=True, help='index folder')
    args = parser.parse_args(argv)
    build(args.gfa, args.variant_index or args.gfa + '.variants.tsv', args.output_folder)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
