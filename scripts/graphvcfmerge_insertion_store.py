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

Each bundle has an index, ``insertions-XX.idx``: one 16-byte record per
committed entry (the first 8 bytes of its key, its offset), appended by the
writer under the bundle lock. A lookup reads it sequentially instead of
walking every entry header; the entry header still carries the full key and
checksum, which every read verifies. Bundles of older runs have no index and
are walked as before.

With LINGRAPH_SHARED_FS_LOCKS=1 (merge_grvcfs.py --slurm sets it), chromosome
jobs on different nodes share the bundles: flock may then lock one node only
and a node may see an old file size, so every bundle access also takes a
cross-node lock (a directory next to the bundle: mkdir is atomic on shared file
systems), opens the bundle only after taking it, writes at the size it then
sees (no O_APPEND) and syncs before releasing it.
"""
from __future__ import annotations

from collections import OrderedDict
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import socket
import struct
import sys
import threading
import time

_HEADER = struct.Struct('<8s32sQQ32s')
_FOOTER = struct.Struct('<Q8s')
_MAGIC = b'ISNPv1\0\0'
_COMMIT = b'ISNPdone'
_INDEXES = OrderedDict()
_IDX = struct.Struct('<QQ')        # key prefix, entry offset
_TABLES = OrderedDict()            # per process: loaded .idx tables
_INDEX_LOCK = threading.RLock()
_LOCAL = threading.local()          # per thread: batch depth and buffered puts
BATCH_RECORD_LIMIT = 1024 * 1024    # larger record blobs are appended at once
BATCH_FLUSH_BYTES = 16 * 1024 * 1024
SHARED_FS_LOCKS = 'LINGRAPH_SHARED_FS_LOCKS'
STALE_LOCK_SECONDS = 1800           # a cross-node lock this old was left by a killed job


def _shared_fs():
    return os.environ.get(SHARED_FS_LOCKS, '') not in ('', '0')


def _break_stale(lock):
    """Remove LOCK when it is older than STALE_LOCK_SECONDS. Renamed first,
    so of several waiters only one removes it."""
    try:
        age = time.time() - os.stat(lock).st_mtime
    except FileNotFoundError:
        return
    if age < STALE_LOCK_SECONDS:
        return
    stale = f'{lock}.stale.{socket.gethostname()}.{os.getpid()}.{threading.get_ident()}'
    try:
        os.rename(lock, stale)
    except OSError:
        return
    os.rmdir(stale)
    print(f'[insertion-store] removed a stale lock ({age:.0f} s old): {lock}', file=sys.stderr, flush=True)


@contextmanager
def _node_lock(path):
    """Cross-node exclusive lock of one bundle under LINGRAPH_SHARED_FS_LOCKS;
    otherwise nothing (flock alone serializes the workers of one node)."""
    if not _shared_fs():
        yield
        return
    lock = f'{path}.lock'
    delay = 0.005
    while True:
        try:
            os.mkdir(lock)
            break
        except FileExistsError:
            _break_stale(lock)
        except FileNotFoundError:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            continue
        time.sleep(delay)
        delay = min(delay * 2, 0.25)
    try:
        yield
    finally:
        os.rmdir(lock)


def _sync(handle):
    """Shared bundles: written bytes reach the file server before the lock goes."""
    if _shared_fs():
        os.fsync(handle.fileno())


def _after_fork():
    global _INDEX_LOCK, _LOCAL
    _INDEX_LOCK = threading.RLock()
    _INDEXES.clear()
    _TABLES.clear()
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
        with _node_lock(path), _open_append(path) as handle:
            fcntl.flock(handle, fcntl.LOCK_EX)
            start = end = _append_end(handle)
            handle.seek(end)
            written = []
            for key, metadata, records, checksum in items:
                handle.write(_HEADER.pack(_MAGIC, bytes.fromhex(key), len(metadata),
                                          len(records), checksum))
                handle.write(metadata)
                handle.write(records)
                handle.write(_FOOTER.pack(end, _COMMIT))
                written.append((bytes.fromhex(key), end))
                end += _HEADER.size + len(metadata) + len(records) + _FOOTER.size
            handle.flush()
            _sync(handle)
            _index_append(path, handle.fileno(), start, written)


def _open_append(path):
    # Created without truncating; not O_APPEND: entries go at the end the
    # writer found under its locks (on a shared file system O_APPEND writes at
    # the end the node believes, which another node's append may have moved),
    # and an incomplete tail can be recovered.
    try:
        os.close(os.open(path, os.O_RDWR | os.O_CREAT, 0o666))
    except FileNotFoundError:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        os.close(os.open(path, os.O_RDWR | os.O_CREAT, 0o666))
    return open(path, 'r+b')


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
    with _node_lock(path), _open_append(path) as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        end = _append_end(handle)
        handle.seek(end)
        handle.write(_HEADER.pack(_MAGIC, bytes.fromhex(key), len(metadata), record_size, checksum))
        handle.write(metadata)
        for blob in chunks():
            handle.write(blob)
        handle.write(_FOOTER.pack(end, _COMMIT))
        handle.flush()
        _sync(handle)
        _index_append(path, handle.fileno(), end, [(bytes.fromhex(key), end)])
        return _source(path, key, (end, checksum.hex()))


def _index_path(path):
    return Path(path).with_suffix('.idx')


def _prefix(key):
    return int.from_bytes(bytes.fromhex(key)[:8], 'little')


def _entry_end(fd, offset):
    magic, _key, meta_size, record_size, _digest = _HEADER.unpack(os.pread(fd, _HEADER.size, offset))
    if magic != _MAGIC:
        raise ValueError(f'insertion index points to no entry at offset {offset}')
    return offset + _HEADER.size + meta_size + record_size + _FOOTER.size


def _walk(fd, start, stop):
    """(key, offset) of the committed entries from START, headers and footers only."""
    found, end = [], start
    while end < stop:
        header = os.pread(fd, _HEADER.size, end)
        if len(header) < _HEADER.size:
            break
        magic, key, meta_size, record_size, _digest = _HEADER.unpack(header)
        if magic != _MAGIC:
            raise ValueError(f'invalid insertion bundle header at offset {end}')
        footer_at = end + _HEADER.size + meta_size + record_size
        if footer_at + _FOOTER.size > stop:
            break
        offset, committed = _FOOTER.unpack(os.pread(fd, _FOOTER.size, footer_at))
        if offset != end or committed != _COMMIT:
            raise ValueError(f'invalid insertion bundle footer at offset {end}')
        found.append((key, end))
        end = footer_at + _FOOTER.size
    return found, end


def _index_append(path, fd, start, written):
    """Index this append (bundle lock held). Entries committed before START
    but never indexed (a writer killed in between) are indexed first."""
    with _open_append(_index_path(path)) as out:
        size = out.seek(0, os.SEEK_END)
        whole = size - size % _IDX.size
        if whole != size:
            out.truncate(whole)              # a killed writer's partial record
        covered = 0
        if whole:
            _prefix_value, last = _IDX.unpack(os.pread(out.fileno(), _IDX.size, whole - _IDX.size))
            try:
                covered = _entry_end(fd, last)
            except (ValueError, struct.error):
                covered = start + 1
        if covered > start:                  # does not describe this bundle: rebuild
            out.truncate(0)
            covered = 0
        missing = _walk(fd, covered, start)[0] if covered < start else []
        out.seek(0, os.SEEK_END)             # after any truncation above
        out.write(b''.join(_IDX.pack(int.from_bytes(key[:8], 'little'), offset)
                           for key, offset in missing + written))
        out.flush()
        _sync(out)


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
        if _shared_fs() and not path.is_file():
            continue
        with _node_lock(path):
            try:
                handle = path.open('rb')
            except FileNotFoundError:
                continue
            with handle:
                fcntl.flock(handle, fcntl.LOCK_SH)
                for key, offset in _locate(path, handle, wanted).items():
                    result[key] = _source(path, key, offset)
    return result


def _locate(path, handle, wanted):
    """{key: (offset, digest)} of WANTED's latest generations in one bundle.

    Through its .idx (digest None: read_entry verifies the key and checksum);
    a bundle of an older run without one is walked."""
    fd = handle.fileno()
    size = os.fstat(fd).st_size
    wanted = set(wanted)
    table = _index_table(path) if _index_path(path).is_file() else None
    if table is None:
        with _INDEX_LOCK:
            _, entries = _index(handle)
        return {key: entries[key] for key in wanted if key in entries}
    import numpy as np
    by_prefix = {}
    for key in wanted:
        by_prefix.setdefault(_prefix(key), []).append(key)
    found = {}
    if len(table):
        want = np.fromiter(by_prefix, dtype='<u8', count=len(by_prefix))
        for row in np.flatnonzero(np.isin(table['key'], want)):     # file order: latest wins
            for key in by_prefix[int(table['key'][row])]:
                found[key] = (int(table['offset'][row]), None)
        covered = _entry_end(fd, int(table['offset'][-1]))
    else:
        covered = 0
    for key, offset in _walk(fd, covered, size)[0]:     # committed, not yet indexed
        if key.hex() in wanted:
            found[key.hex()] = (offset, None)
    return found


def _index_table(path):
    """A process's copy of PATH's .idx, extended by the records appended since."""
    import numpy as np
    dtype = np.dtype([('key', '<u8'), ('offset', '<u8')])
    index = _index_path(path)
    count = os.stat(index).st_size // _IDX.size
    with _INDEX_LOCK:
        identity = (os.getpid(), str(index))
        table = _TABLES.pop(identity, np.empty(0, dtype=dtype))
        if count < len(table):
            table = np.empty(0, dtype=dtype)     # rebuilt by a writer
        if count > len(table):
            tail = np.fromfile(index, dtype=dtype, count=count - len(table),
                               offset=len(table) * _IDX.size)
            table = np.concatenate([table, tail])
        _TABLES[identity] = table
        while len(_TABLES) > 64:
            _TABLES.popitem(last=False)
        return table


def sources(root):
    """Discover the latest committed generation of every stored locus."""
    flush()
    result = []
    for path in sorted(Path(root).glob('insertions-*.bundle')):
        with _node_lock(path), path.open('rb') as handle:
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
        with _node_lock(path), open(path, 'rb') as handle:
            fcntl.flock(handle, fcntl.LOCK_SH)
            size = os.fstat(handle.fileno()).st_size
            for index in sorted(indexes, key=lambda item: sources[item]['offset']):
                result[index] = _read_entry(handle, path, size, sources[index], with_records)
    return result


class DirectReader:
    """Lock-free reads of committed entries at known offsets.

    Each bundle is opened once and every batch is read in file order. Safe
    while other workers append: committed entries never change (writers only
    append, and trim an unfinished tail after the last commit), every read
    checks the entry's footer and checksum, and offsets found before a bundle
    is opened lie within its size."""

    def __init__(self):
        self._handles = {}

    def __enter__(self):
        flush()    # as read_many: this thread's buffered puts first
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
            or source['digest'] is not None and digest.hex() != source['digest']):
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
