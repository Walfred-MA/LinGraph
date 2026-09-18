"""Query FASTA coordinates and bounded parallel interval reading."""
from functools import lru_cache
from pathlib import Path
import multiprocessing
import re
import shlex


READ_CHUNK_BASES = 1024 * 1024


def reverse_complement(sequence):
    return sequence.translate(str.maketrans(
        'ACGTRYMKBDHVNacgtrymkbdhvn',
        'TGCAYRKMVHDBNtgcayrkmvhdbn'))[::-1]


def read_query_list(path):
    root = Path(path).resolve().parent
    sources = {}
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            fields = shlex.split(line, comments=True)
            if not fields:
                continue
            if len(fields) not in (2, 3):
                raise ValueError(f'{path}:{number}: expected NAME FASTA [FAI]')
            sample, value = fields[:2]
            fasta = Path(value).expanduser()
            if not fasta.is_absolute():
                fasta = root / fasta
            fai = Path(fields[2]).expanduser() if len(fields) == 3 else Path(str(fasta)+'.fai')
            if not fai.is_absolute():
                fai = root / fai
            if sample in sources:
                raise ValueError(f'duplicate query assembly: {sample}')
            if not fasta.is_file():
                raise FileNotFoundError(fasta)
            sources[sample] = str(fasta), str(fai) if fai.is_file() else None
    if not sources:
        raise ValueError('empty query FASTA list')
    return sources


def _unescape(text):
    return (text.replace('@0A', '\n').replace('@09', '\t').replace('@2C', ',')
            .replace('@3B', ';').replace('@3A', ':').replace('@40', '@'))


@lru_cache(maxsize=64)
def _format_indexes(format_text):
    return {name: index for index, name in enumerate(format_text.split(':'))}


def parse_sample_observations(sample, text, format_text):
    indexes = _format_indexes(format_text)
    if not text or text == '.' or (indexes.get('GT') == 0 and text.startswith(('0:', '.:'))):
        return
    values = text.split(':')
    if len(values) != len(indexes):
        return

    def column(name):
        index = indexes.get(name)
        return values[index].split(',') if index is not None else ['.']

    types, sizes = column('TYPE'), column('SIZE')
    coords, contigs = column('QUERYCOORD'), column('ASSEMBLYCONTIG')
    cigars, labels = column('EXTENDGRAPHCIGAR'), column('ALLELENAME')

    def field(items, index):
        return _unescape(items[index if len(items) > 1 else 0])

    count = max(map(len, (types, sizes, coords, contigs, cigars, labels)))
    if any(len(items) not in (1, count)
           for items in (types, sizes, coords, contigs, cigars, labels)):
        raise ValueError(f'{sample}: inconsistent per-observation FORMAT fields')
    for index in range(count):
        if field(types, index) != 'INS':
            continue
        coordinate = re.fullmatch(r'(\d+)(?:-(\d+))?([+-]?)', field(coords, index))
        if coordinate is None:
            continue
        start, end = sorted((int(coordinate[1]), int(coordinate[2] or coordinate[1])))
        yield (sample, field(contigs, index), start, end, coordinate[3] or '+',
               field(cigars, index), abs(int(field(sizes, index))), field(labels, index))


def parse_sample_query_intervals(sample, text, format_text, *, include_alignment=False):
    """Yield only fields needed to choose a representative query interval."""
    indexes = _format_indexes(format_text)
    if not text or text == '.' or (indexes.get('GT') == 0 and text.startswith(('0:', '.:'))):
        return
    values = text.split(':')
    if len(values) != len(indexes):
        return

    def column(name):
        index = indexes.get(name)
        return values[index].split(',') if index is not None else ['.']

    types, sizes = column('TYPE'), column('SIZE')
    coords, contigs = column('QUERYCOORD'), column('ASSEMBLYCONTIG')
    cigars, offsets = (column('EXTENDGRAPHCIGAR'), column('TEMPLATEOFFSET')) if include_alignment else ([], [])

    def field(items, index):
        return _unescape(items[index if len(items) > 1 else 0])

    columns = (types, sizes, coords, contigs)
    if include_alignment:
        columns += (cigars, offsets)
    count = max(map(len, columns))
    if any(len(items) not in (1, count) for items in columns):
        raise ValueError(f'{sample}: inconsistent per-observation query fields')
    for index in range(count):
        if field(types, index) != 'INS':
            continue
        coordinate = re.fullmatch(r'(\d+)(?:-(\d+))?([+-]?)', field(coords, index))
        if coordinate is None:
            continue
        start, end = sorted((int(coordinate[1]), int(coordinate[2] or coordinate[1])))
        interval = (sample, field(contigs, index), start, end, coordinate[3] or '+',
                    abs(int(field(sizes, index))))
        yield interval + (field(cigars, index), field(offsets, index)) if include_alignment else interval


def interval_chunks(reader, contig, start, end, strand):
    if strand == '+':
        for offset in range(start, end, READ_CHUNK_BASES):
            yield reader.fetch(contig, offset, min(end, offset+READ_CHUNK_BASES))
    else:
        for stop in range(end, start, -READ_CHUNK_BASES):
            yield reverse_complement(reader.fetch(
                contig, max(start, stop-READ_CHUNK_BASES), stop))


def scaffold_pool(workers):
    # Each worker exits after one batch so its FASTA index and allocator pages
    # are returned to the operating system immediately.
    return multiprocessing.get_context('spawn').Pool(workers, maxtasksperchild=1)
