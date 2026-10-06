"""Indexed insertion spools in 32 append-only bundles, readable across workers.

Each committed entry contains metadata, packed SNP records and a checksum.
An advisory file lock serializes writers of one bundle. A footer publishes the
entry; an interrupted tail is ignored by readers and removed by the next writer.
Retries append a new generation. Existing checkpoint entries remain independent
of unrelated appends. Durability is supplied by the chromosome checkpoint flush.

File operations stay per bundle, not per locus: inside ``batch()`` a thread's
small puts are buffered and appended with one open and lock per bundle;
``read_many`` reads many entries through one open per bundle. A lookup first
publishes the calling thread's buffer, so a thread always finds its own puts.
"""
from __future__ import annotations

from collections import OrderedDict
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import struct
import threading

_HEADER = struct.Struct('<8s32sQQ32s')
_FOOTER = struct.Struct('<Q8s')
_MAGIC = b'ISNPv1\0\0'
_COMMIT = b'ISNPdone'
_INDEXES = OrderedDict()
_INDEX_LOCK = threading.RLock()
_LOCAL = threading.local()          # per thread: batch depth and buffered puts
BATCH_RECORD_LIMIT = 1024 * 1024    # larger record blobs are appended at once
BATCH_FLUSH_BYTES = 16 * 1024 * 1024


def _after_fork():
    global _INDEX_LOCK, _LOCAL
    _INDEX_LOCK = threading.RLock()
    _INDEXES.clear()
    # A child must not publish its parent's buffered puts a second time.
    _LOCAL = threading.local()


os.register_at_fork(after_in_child=_after_fork)


def _path(root, key):
    raw = bytes.fromhex(key)
    if len(raw) != 32:
        raise ValueError('insertion key must be a SHA256 digest')
    return Path(root) / f'insertions-{raw[0] % 32:02x}.bundle'


def _index(handle):
    """Refresh a process's header-only index while holding the file lock."""
    stat = os.fstat(handle.fileno())
    identity = (os.getpid(), os.path.abspath(handle.name), stat.st_dev, stat.st_ino)
    end, entries = _INDEXES.pop(identity, (0, {}))
    if stat.st_size < end:
        end, entries = 0, {}
    while end < stat.st_size:
        handle.seek(end)
        header = handle.read(_HEADER.size)
        if len(header) < _HEADER.size:
            break
        magic, key, meta_size, record_size, digest = _HEADER.unpack(header)
        if magic != _MAGIC:
            raise ValueError(f'invalid insertion bundle header: {handle.name}:{end}')
        stop = end + _HEADER.size + meta_size + record_size
        if stop + _FOOTER.size > stat.st_size:
            break
        handle.seek(stop)
        offset, committed = _FOOTER.unpack(handle.read(_FOOTER.size))
        if offset != end or committed != _COMMIT:
            raise ValueError(f'invalid insertion bundle footer: {handle.name}:{end}')
        entries[key.hex()] = (end, digest.hex())
        end = stop + _FOOTER.size
    _INDEXES[identity] = (end, entries)
    while len(_INDEXES) > 64:
        _INDEXES.popitem(last=False)
    return end, entries


def _source(path, key, entry):
    offset, digest = entry
    return dict(bundle=os.path.abspath(path), key=key, offset=offset, digest=digest)


def _append_end(handle):
    """Normal appends inspect only the last entry, not other workers' indexes."""
    size = os.fstat(handle.fileno()).st_size
    if size == 0:
        return 0
    if size >= _HEADER.size + _FOOTER.size:
        start, committed = _FOOTER.unpack(os.pread(handle.fileno(), _FOOTER.size,
                                                   size - _FOOTER.size))
        if committed == _COMMIT and start <= size - _HEADER.size - _FOOTER.size:
            magic, _, meta_size, record_size, _ = _HEADER.unpack(
                os.pread(handle.fileno(), _HEADER.size, start))
            if magic == _MAGIC and start + _HEADER.size + meta_size + record_size + _FOOTER.size == size:
                return size
    # Only a killed/incomplete append needs a scan to the last committed footer.
    with _INDEX_LOCK:
        end, _ = _index(handle)
    handle.truncate(end)
    return end


def _pending():
    pending = getattr(_LOCAL, 'pending', None)
    if pending is None:
        pending = _LOCAL.pending = []
        _LOCAL.depth = 0
        _LOCAL.size = 0
    return pending


@contextmanager
def batch():
    """Buffer this thread's puts; append them per bundle when the outermost
    batch ends. A failed task's buffered puts are dropped (a retry appends a
    new generation)."""
    _pending()
    _LOCAL.depth += 1
    try:
        yield
    except BaseException:
        _LOCAL.depth -= 1
        if not _LOCAL.depth:
            _LOCAL.pending.clear()
            _LOCAL.size = 0
        raise
    _LOCAL.depth -= 1
    if not _LOCAL.depth:
        flush()


def flush():
    """Append this thread's buffered puts: one open and lock per bundle."""
    pending = _pending()
    if not pending:
        return
    entries, _LOCAL.pending, _LOCAL.size = list(pending), [], 0
    by_bundle = OrderedDict()
    for path, key, metadata, records, checksum in entries:
        by_bundle.setdefault(path, []).append((key, metadata, records, checksum))
    for path, items in by_bundle.items():
        handle = _open_append(path)
        with handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            end = _append_end(handle)
            handle.seek(end)
            for key, metadata, records, checksum in items:
                handle.write(_HEADER.pack(_MAGIC, bytes.fromhex(key), len(metadata),
                                          len(records), checksum))
                handle.write(metadata)
                handle.write(records)
                handle.write(_FOOTER.pack(end, _COMMIT))
                end += _HEADER.size + len(metadata) + len(records) + _FOOTER.size
            handle.flush()


def _open_append(path):
    # a+b creates without truncating and permits recovery of an incomplete tail.
    try:
        return path.open('a+b')
    except FileNotFoundError:
        path.parent.mkdir(parents=True, exist_ok=True)
        return path.open('a+b')


def put(root, key, data, records=b''):
    """Publish one locus, including empty loci with coverage but no SNPs.

    Inside ``batch()`` a locus with records up to BATCH_RECORD_LIMIT bytes is
    buffered (returns None) and appended when the batch ends."""
    path = _path(root, key)
    metadata = json.dumps(data, separators=(',', ':'), sort_keys=True).encode()
    streamed = hasattr(records, 'read')
    record_size = records.seek(0, os.SEEK_END) if streamed else len(records)

    def chunks():
        if streamed:
            records.seek(0)
            while blob := records.read(1024 * 1024):
                yield blob
        else:
            yield records

    digest = hashlib.sha256(metadata)
    for blob in chunks():
        digest.update(blob)
    checksum = digest.digest()
    if getattr(_LOCAL, 'depth', 0) and record_size <= BATCH_RECORD_LIMIT:
        blob = b''.join(chunks())
        _pending().append((path, key, metadata, blob, checksum))
        _LOCAL.size += len(metadata) + len(blob)
        if _LOCAL.size >= BATCH_FLUSH_BYTES:
            flush()
        return None
    handle = _open_append(path)
    with handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        end = _append_end(handle)
        handle.seek(end)
        handle.write(_HEADER.pack(_MAGIC, bytes.fromhex(key), len(metadata), record_size, checksum))
        handle.write(metadata)
        for blob in chunks():
            handle.write(blob)
        handle.write(_FOOTER.pack(end, _COMMIT))
        handle.flush()
        return _source(path, key, (end, checksum.hex()))


def find(root, key):
    return find_many(root, (key,)).get(key)


def find_many(root, keys):
    """Resolve a chromosome's loci with at most one open per bundle."""
    flush()     # this thread's buffered puts are findable
    grouped = {}
    for key in keys:
        grouped.setdefault(_path(root, key), []).append(key)
    result = {}
    for path, wanted in grouped.items():
        try:
            handle = path.open('rb')
        except FileNotFoundError:
            continue
        with handle:
            fcntl.flock(handle, fcntl.LOCK_SH)
            with _INDEX_LOCK:
                _, entries = _index(handle)
                for key in wanted:
                    if key in entries:
                        result[key] = _source(path, key, entries[key])
    return result


def sources(root):
    """Discover the latest committed generation of every stored locus."""
    flush()
    result = []
    for path in sorted(Path(root).glob('insertions-*.bundle')):
        with path.open('rb') as handle:
            fcntl.flock(handle, fcntl.LOCK_SH)
            with _INDEX_LOCK:
                _, entries = _index(handle)
                result.extend(_source(path, key, entry) for key, entry in entries.items())
    return result


def read(source, *, with_records=True):
    """Read a specific generation; verify its bytes when loading SNP records."""
    return read_many([source], with_records=with_records)[0]


def read_many(sources, *, with_records=True):
    """Read many generations, in the given order, opening each bundle once."""
    flush()
    by_bundle = OrderedDict()
    for index, source in enumerate(sources):
        by_bundle.setdefault(source['bundle'], []).append(index)
    result = [None] * len(sources)
    for path, indexes in by_bundle.items():
        with open(path, 'rb') as handle:
            fcntl.flock(handle, fcntl.LOCK_SH)
            size = os.fstat(handle.fileno()).st_size
            for index in sorted(indexes, key=lambda item: sources[item]['offset']):
                result[index] = _read_entry(handle, path, size, sources[index], with_records)
    return result


class DirectReader:
    """Lock-free reads of committed entries at known offsets.

    Only for stores no worker writes any more (a backfill after its SV merge):
    each bundle is opened once and every batch is read in file order."""

    def __init__(self):
        self._handles = {}

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        for handle, _size in self._handles.values():
            handle.close()
        self._handles.clear()

    def read_many(self, sources, *, with_records=True):
        result = [None] * len(sources)
        for index in sorted(range(len(sources)),
                            key=lambda item: (sources[item]['bundle'], sources[item]['offset'])):
            path = sources[index]['bundle']
            if path not in self._handles:
                handle = open(path, 'rb')
                self._handles[path] = (handle, os.fstat(handle.fileno()).st_size)
            handle, size = self._handles[path]
            result[index] = _read_entry(handle, path, size, sources[index], with_records)
        return result


def _read_entry(handle, path, size, source, with_records):
    """One committed entry through an open, share-locked bundle handle."""
    start = source['offset']
    handle.seek(start)
    header = handle.read(_HEADER.size)
    if len(header) != _HEADER.size:
        raise ValueError(f'truncated insertion bundle: {path}')
    magic, key, meta_size, record_size, digest = _HEADER.unpack(header)
    if (magic != _MAGIC or key.hex() != source['key']
            or digest.hex() != source['digest']):
        raise ValueError(f'changed insertion bundle entry: {path}:{start}')
    if start + _HEADER.size + meta_size + record_size + _FOOTER.size > size:
        raise ValueError(f'truncated insertion bundle entry: {path}:{start}')
    metadata = handle.read(meta_size)
    records = handle.read(record_size) if with_records else b''
    if not with_records:
        handle.seek(record_size, os.SEEK_CUR)
    footer = handle.read(_FOOTER.size)
    if (len(metadata) != meta_size or len(footer) != _FOOTER.size
            or _FOOTER.unpack(footer) != (start, _COMMIT)):
        raise ValueError(f'truncated insertion bundle entry: {path}:{start}')
    if with_records:
        checksum = hashlib.sha256(metadata)
        checksum.update(records)
        if len(records) != record_size or checksum.digest() != digest:
            raise ValueError(f'insertion bundle checksum mismatch: {path}:{start}')
    return json.loads(metadata), records
