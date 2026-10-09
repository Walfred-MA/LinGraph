"""Small-variant C++ worker bridge; manifests/headers remain Python-owned.

The binary request contains only dictionaries, coverage and input spans, never
expanded observations. Existing 21-byte shards and saved workflows stay usable.
GRAPHVCFMERGE_BACKEND=python selects the reference implementation for comparison.
"""
from __future__ import annotations

from contextlib import contextmanager
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import struct
import subprocess
import tempfile


def enabled():
    backend = os.environ.get('GRAPHVCFMERGE_BACKEND', 'native')
    if backend not in ('native', 'python'):
        raise ValueError('GRAPHVCFMERGE_BACKEND must be native or python')
    return backend == 'native'


def memory_bytes(sample_count, divisor=1):
    value = os.environ.get('GRAPHVCFMERGE_RAM_GB')
    # 128 GiB covers the target 1,000-2,000-sample cohort. Beyond that,
    # grow in 64-GiB steps to 512 GiB at 10,000 samples.
    gib = float(value) if value else (64 if sample_count < 100 else
                                    max(128, math.ceil(sample_count / 1280) * 64))
    if not math.isfinite(gib) or gib <= 0:
        raise ValueError('GRAPHVCFMERGE_RAM_GB must be positive and finite')
    budget = int(gib * (1 << 30))
    # Saved Slurm workflows may still request 64G for a larger cohort. Respect
    # that allocation even when using the newer sample-count default.
    allocated = os.environ.get('SLURM_MEM_PER_NODE')
    if allocated:
        limit = int(allocated) << 20
        if limit:
            budget = min(budget, limit)
    elif os.environ.get('SLURM_MEM_PER_CPU'):
        limit = int(os.environ['SLURM_MEM_PER_CPU']) * int(os.environ.get('SLURM_CPUS_PER_TASK', '1')) << 20
        if limit:
            budget = min(budget, limit)
    return max(1 << 20, budget // max(1, divisor))


_BINARY = None


def executable():
    global _BINARY
    if _BINARY is not None:
        return _BINARY
    scripts = Path(__file__).resolve().parent
    source = scripts / 'src' / 'SmallVariantMerge'
    files = [source / 'SmallVariantMerge.cpp', source / 'support.hpp']
    for candidate in (scripts / 'SmallVariantMerge', source / 'SmallVariantMerge'):
        if (candidate.is_file() and os.access(candidate, os.X_OK)
                and all(candidate.stat().st_mtime_ns >= p.stat().st_mtime_ns for p in files)):
            _BINARY = str(candidate)
            return _BINARY
    digest = hashlib.sha256(b''.join(p.read_bytes() for p in files)).hexdigest()[:20]
    cache = Path(tempfile.gettempdir()) / f'graphvcfmerge-native-{os.getuid()}' / digest
    cache.mkdir(parents=True, exist_ok=True)
    binary = cache / 'SmallVariantMerge'
    # Parallel scan workers compile once. Publish only a successful executable.
    with (cache / 'build.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not binary.is_file():
            temporary = cache / f'build.{os.getpid()}'
            command = [*shlex.split(os.environ.get('CXX', 'c++')), '-O3', '-std=c++17',
                       '-pthread', str(files[0]), '-o', str(temporary), '-lz']
            if platform.machine().lower() in ('x86_64', 'amd64'):
                # Baseline x86-64: AVX from a host-tuned default dies with SIGILL on older nodes.
                command[1:1] = ['-march=x86-64', '-mtune=generic']
            prefix = os.environ.get('CONDA_PREFIX')
            if prefix and (Path(prefix) / 'include' / 'zlib.h').is_file():
                command += ['-I' + str(Path(prefix) / 'include'), '-L' + str(Path(prefix) / 'lib'),
                            '-Wl,-rpath,' + str(Path(prefix) / 'lib')]
            try:
                result = subprocess.run(command, capture_output=True, text=True)
                if result.returncode:
                    raise RuntimeError('cannot build SmallVariantMerge (requires C++17 and zlib):\n'
                                       + result.stderr)
                os.replace(temporary, binary)
            finally:
                temporary.unlink(missing_ok=True)
    _BINARY = str(binary)
    return _BINARY


class Request:
    def __init__(self, handle):
        self.handle = handle
        self.string('small-variant-merge-v1')

    def u32(self, value):
        self.handle.write(struct.pack('<I', value))

    def u64(self, value):
        self.handle.write(struct.pack('<Q', value))

    def string(self, value):
        value = str(value).encode('utf-8')
        self.u32(len(value))
        self.handle.write(value)

    def strings(self, values):
        self.u32(len(values))
        for value in values:
            self.string(value)


@contextmanager
def request_file(directory):
    with tempfile.NamedTemporaryFile(prefix='.native-', dir=directory) as handle:
        yield Request(handle), handle


def invoke(command, request, *outputs):
    request.flush()
    result = subprocess.run([executable(), command, request.name, *map(str, outputs)],
                            capture_output=True, text=True)
    if result.returncode:
        raise ValueError(result.stderr.strip() or f'SmallVariantMerge exited {result.returncode}')


def scan_one(index, paths, directory, cutoff=20, *, compact=True):
    target = Path(directory) / str(index)
    target.mkdir()
    output = target / 'scan-result.json'
    with request_file(target) as (request, handle):
        request.string(target)
        request.u32(cutoff)
        request.strings(paths)
        invoke('scan-snp' if compact else 'scan-small', handle, output)
    result = json.loads(output.read_text())
    output.unlink()
    if compact:
        return result
    return dict(index=index, samples=result['samples'], chroms=result['chroms'], meta=result['metadata'])


def _realignment(root, key):
    if root is None:
        return dict(drops=[], adds=[])
    import graphvcfmerge_insertion_store as store
    path = Path(root) / 'realign' / (key + '.json')
    source = store.find(path.parent, key)
    if source is not None:
        return store.read(source)[0]
    return json.loads(path.read_text()) if path.is_file() else dict(drops=[], adds=[])


def compact_locus(manifest, key, processes, output, root=None):
    from graphvcfmerge_snp_compact import source_sections
    import graphvcfmerge as vcf
    paths = list(dict.fromkeys(manifest['sources_by_chrom'].get(key, [])))
    chrom = manifest['loci'][key]
    realignment = _realignment(root, key)
    with request_file(Path(output).parent) as (request, handle):
        request.string(chrom)
        request.string(vcf._safe_variant_token(chrom))
        request.u64(memory_bytes(len(manifest['samples'])))
        request.u32(max(1, processes))
        request.strings(manifest['samples'])
        for kind in ('drops', 'adds'):
            request.u32(len(realignment[kind]))
            for row in realignment[kind]:
                for index, value in enumerate(row):
                    if index == 1 or index == (4 if kind == 'drops' else 5):
                        request.u32(value)
                    else:
                        request.string(value)
        count_offset = handle.tell()
        request.u64(0)
        # Coverage-only source dictionaries can legitimately have no binary.
        request.u32(len(paths))
        total = 0
        for path in paths:
            data = json.loads(Path(path).read_text())
            binary, spans = source_sections(data, key)
            request.u32(len(data['queries']))
            for query in data['queries']:
                values = list(query)
                values[5] = int(values[5])
                request.strings(values)
            request.strings(data['alleles'])
            request.strings(data.get('label_groups', []))
            request.string(binary)
            request.u32(len(spans))
            for offset, count in spans:
                request.u64(offset)
                request.u64(count)
                total += count
            coverage = data['coverage'].get(chrom, [])
            request.u32(len(coverage))
            for sample, start, end in coverage:
                request.string(sample)
                request.u64(start)
                request.u64(end)
        handle.seek(count_offset)
        request.u64(total)
        handle.seek(0, os.SEEK_END)
        invoke('compact', handle, output)


def compact_chrom(manifest, key, processes, root=None):
    from graphvcfmerge_snp import atomic_output, part_path
    parts = Path(manifest['directory']) / 'parts'
    parts.mkdir(exist_ok=True)
    # Body files stay plain; final gzip handling remains in concat.
    with tempfile.TemporaryDirectory(prefix='.native-chrom-', dir=parts) as work:
        bodies = []
        for index, locus in enumerate(manifest['job_loci'][key]):
            body = Path(work) / str(index)
            compact_locus(manifest, locus, processes, body, root)
            bodies.append(body)
        if len(bodies) == 1:
            os.replace(bodies[0], part_path(manifest, key, 'snp'))
        else:
            with atomic_output(part_path(manifest, key, 'snp')) as handle:
                for body in bodies:
                    with body.open() as source:
                        shutil.copyfileobj(source, handle, 1 << 20)
    with atomic_output(part_path(manifest, key, 'indel')):
        pass
    return key


def small_chrom(manifest, key, memory_divisor=1):
    from graphvcfmerge_snp import part_path
    import graphvcfmerge as vcf
    chrom = manifest['chroms'][key]
    parts = Path(manifest['directory']) / 'parts'
    parts.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.native-small-', dir=parts) as work:
        with request_file(work) as (request, handle):
            request.string(chrom)
            request.string(vcf._safe_variant_token(chrom))
            request.u32(manifest['cutoff'])
            request.u64(memory_bytes(len(manifest['samples']), memory_divisor))
            request.strings(manifest['samples'])
            sources = [s for s in manifest['files'] if key in s['chroms']]
            request.u32(len(sources))
            for source in sources:
                request.string(Path(manifest['directory']) / str(source['index']) / (key + '.vcf'))
                request.strings(source['samples'])
            outputs = [Path(work) / kind for kind in ('snp', 'indel')]
            invoke('small', handle, *outputs)
        for kind, path in zip(('snp', 'indel'), outputs):
            os.replace(path, part_path(manifest, key, kind))
    return key
