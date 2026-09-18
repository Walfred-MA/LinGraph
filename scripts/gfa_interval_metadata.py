"""Small, streaming metadata reader for the query-backed GFA converter."""
import gzip
import hashlib
import json
import re

from gfa_query_anchors import _unescape, parse_sample_query_intervals
from assembly_contigs import accepted_contig

_CHUNK = re.compile(r'[<>]([^<>]*)')
_OP = re.compile(r'(\d+)([=MXIDHS])')
_PAYLOAD = re.compile(r'[A-Za-z]+')


def graph_runs(text, size):
    """Return query runs without retaining sequence or CIGAR payload text.

    A run is ``(qstart, qend, operation, target, rstart, rend,
    orientation, unique_ordinal)``. For reverse mappings, rstart/rend are
    distances from the target's right edge and are normalized later.
    """
    if not text or text == '.' or text[0] not in '<>':
        return [(0, size, 'I', None, 0, 0, '+', 0)]

    output = []
    qpos = unique = previous = 0
    for chunk in _CHUNK.finditer(text):
        if chunk.start() != previous:
            raise ValueError('text outside an encoded graph chunk')
        previous = chunk.end()
        start, end = chunk.start(1), chunk.end(1)
        delimiters = []
        literal = text.find(':', start, end)
        escaped = text.find('@3A', start, end)
        if literal >= 0:
            delimiters.append((literal, 1))
        if escaped >= 0:
            delimiters.append((escaped, 3))
        if delimiters:
            delimiter, width = min(delimiters)
            target = _unescape(text[start:delimiter])
            position = delimiter + width
        else:
            target, position = None, start

        orientation = '+' if text[chunk.start()] == '>' else '-'
        rpos = 0
        first = True
        while position < end:
            match = _OP.match(text, position, end)
            if match is None:
                raise ValueError('malformed graph mapping operation')
            length, operation = int(match[1]), match[2]
            operation = '=' if operation == 'M' else operation
            if first and operation == 'H':
                rpos = length
            first = False
            if operation in 'IS':
                unique += 1
            if operation in '=X' and target is None:
                raise ValueError(f'unnamed {operation} operation has no graph target')
            if operation in '=XIS' and length:
                output.append((qpos, qpos + length, operation, target,
                               rpos, rpos + length if operation in '=X' else rpos,
                               orientation, unique if operation in 'IS' else 0))
                qpos += length
            if operation in '=XD':
                rpos += length
            position = match.end()
            payload = _PAYLOAD.match(text, position, end)
            if payload:
                position = payload.end()

    if previous != len(text) or qpos != size:
        raise ValueError(f'encoded query length {qpos} differs from SVLEN {size}')
    return output


def vcf_paths(value):
    if isinstance(value, (list, tuple)):
        for item in value:
            yield from vcf_paths(item)
    else:
        yield value


def read_header_contigs(path):
    """Combine headers in input order, sharing contigs across VCF files."""
    ordinals, lengths = {}, {}
    for item in vcf_paths(path):
        _ordinals, current = _read_header_contigs(item)
        for name, length in current.items():
            if name in lengths and lengths[name] is not None and length is not None and lengths[name] != length:
                raise ValueError(f'conflicting contig lengths for {name!r} in {item}')
            if name not in ordinals:
                ordinals[name] = len(ordinals) + 1
            if name not in lengths or lengths[name] is None:
                lengths[name] = length
    return ordinals, lengths


def _read_header_contigs(path):
    """Return original ##contig ordinals and declared lengths."""
    opener = gzip.open if str(path).endswith('.gz') else open
    ordinal = 0
    ordinals, lengths = {}, {}
    with opener(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#CHROM') or not line.startswith('#'):
                break
            if not line.startswith('##contig=<'):
                continue
            ordinal += 1
            identifier = re.search(r'(?:<|,)ID=("(?:[^"\\]|\\.)*"|[^,>]+)', line)
            if identifier is None:
                raise ValueError(f'missing ID in ##contig entry {ordinal}')
            name = identifier[1]
            if name.startswith('"'):
                name = json.loads(name)
            length = re.search(r'(?:<|,)length=(\d+)(?=[,>])', line, re.IGNORECASE)
            ordinals[name] = ordinal
            lengths[name] = int(length[1]) if length else None
    return ordinals, lengths


def literal_checks(text, size):
    """Hash declared literal bases without keeping cohort sequence in RAM."""
    if not text or text == '.':
        raise ValueError('rGFA requires INFO/SEQ to verify the representative insertion')
    if text[0] not in '<>':
        if len(text) != size or not _PAYLOAD.fullmatch(text):
            raise ValueError('plain INFO/SEQ length/alphabet differs from SVLEN')
        return ((0, size, hashlib.sha256(text.upper().encode('ascii')).digest()),)
    checks = []
    qpos = 0
    for chunk in _CHUNK.finditer(text):
        body = chunk[1]
        # Target names can contain digits: strip them before finding CIGAR ops.
        body = re.split(r':|@3A', body, maxsplit=1)[-1]
        position = 0
        while position < len(body):
            match = _OP.match(body, position)
            if match is None:
                raise ValueError('malformed graph mapping operation')
            length, operation = int(match[1]), match[2]
            position = match.end()
            if operation == 'M':
                raise ValueError('ambiguous M in INFO/SEQ: use explicit =/X operations '
                                 'or plain sequence for rGFA')
            payload = _PAYLOAD.match(body, position)
            if payload:
                if operation in 'XIS':
                    if len(payload[0]) != length:
                        raise ValueError(f'INFO/SEQ {operation} payload length differs from CIGAR')
                    checks.append((qpos, qpos+length,
                                   hashlib.sha256(payload[0].upper().encode('ascii')).digest()))
                position = payload.end()
            elif operation in 'XIS' and length:
                raise ValueError(f'INFO/SEQ {operation} is missing literal bases needed '
                                 'to verify the representative insertion')
            if operation in '=MXIS':
                qpos += length
    return tuple(checks)


def _other_variant(fields, info, samples, samples_available, cutoff,
                   sequence_checks, candidate_sink, event_kinds):
    """SNP/SUB alleles use verified assembly bases; DEL alleles join parent flanks."""
    from graphvcfmerge import split_sample_alleles, parse_hsv_allele

    svtype = info.get('SVTYPE', '')
    if not svtype and len(fields[3]) == len(fields[4]) == 1 and fields[4].isalpha():
        svtype = 'SNP'
    kind = {'SNP': 'snp', 'DEL': 'deletion', 'SUB': 'substitution'}.get(svtype)
    if kind is None:
        raise ValueError(f'unsupported variant type {svtype or fields[4]!r}')
    identifier, pos = fields[2], int(fields[1])
    end = pos if kind == 'snp' else int(info['END'])
    if kind == 'deletion':
        # END in a merged DEL row is the envelope of all cluster members.
        # SVLEN describes the representative at POS, whose graph branch must
        # skip only that many parent bases. Keep its reference consumption in
        # the D run so the resolver can apply the nested POS origin first.
        deleted_length = 0
        if info.get('SVLEN', '.') != '.':
            deleted_length = abs(int(info['SVLEN']))
            if not deleted_length:
                raise ValueError('deletion SVLEN must be nonzero')
        elif end <= pos:
            raise ValueError('deletion END must be greater than POS')
        size, runs, query, reason = 0, [(0, 0, 'D', None, 0, deleted_length, '+', 0)], None, 'no query bases for deletion'
        measure = deleted_length or end - pos
    else:
        sequence = fields[4] if kind == 'snp' else _unescape(info.get('SEQ', ''))
        if not sequence or not sequence.isalpha():
            raise ValueError('SNP/SUB requires a literal ALT/INFO/SEQ allele')
        size = len(sequence)
        measure = max(size, end - pos)
        runs = [(0, size, 'I', None, 0, 0, '+', 0)]
        if sequence_checks is not None:
            sequence_checks[identifier] = literal_checks(sequence, size)
        candidates = []
        if len(fields) == 10:
            for sample, text in zip(samples, fields[9].split('\t')):
                if sample not in samples_available:
                    continue
                for allele in split_sample_alleles(text, fields[8]):
                    values = parse_hsv_allele(allele)
                    if values is None or values[1] != svtype or not accepted_contig(sample, values[4]):
                        continue
                    if kind == 'snp' and values[3].upper() != sequence.upper():
                        continue
                    coordinate = re.fullmatch(r'(\d+)(?:-(\d+))?([+-]?)', values[5])
                    if coordinate is None:
                        continue
                    low = int(coordinate[1])
                    high = int(coordinate[2]) if coordinate[2] else low + (1 if kind == 'snp' else 0)
                    low, high = sorted((low, high))
                    if high - low == size:
                        candidate = (sample, values[4], low, high, coordinate[3] or '+')
                        if candidate not in candidates:
                            candidates.append(candidate)
        query = candidates[0] if candidates else None
        reason = 'ok' if query else f'coordinate_mismatch: no {svtype} allele has a usable query interval'
        if candidate_sink is not None:
            candidate_sink(identifier, candidates)
    if event_kinds is not None:
        event_kinds[identifier] = kind
    return (identifier, fields[0], pos, end, size, runs, query, reason,
            fields[6] in ('PASS', '.') and measure >= cutoff)


def rows(path, samples_available, cutoff, sequence_checks=None, candidate_sink=None,
         *, insertion_only=False, event_kinds=None):
    for item in vcf_paths(path):
        yield from _rows(item, samples_available, cutoff, sequence_checks, candidate_sink,
                         insertion_only=insertion_only, event_kinds=event_kinds)


def _rows(path, samples_available, cutoff, sequence_checks=None, candidate_sink=None,
          *, insertion_only=False, event_kinds=None):
    """Stream candidates, preferring full representative matches in rGFA mode."""
    opener = gzip.open if str(path).endswith('.gz') else open
    samples = ()
    with opener(path, 'rt') as handle:
        for lineno, line in enumerate(handle, 1):
            if line.startswith('#CHROM\t'):
                samples = tuple(line.rstrip('\r\n').split('\t')[9:])
                continue
            if line.startswith('#') or not line.strip():
                continue
            try:
                fields = line.rstrip('\r\n').split('\t', 9)
                if len(fields) < 8:
                    raise ValueError('malformed VCF record')
                info = dict(item.split('=', 1) for item in fields[7].split(';') if '=' in item)
                if info.get('SVTYPE') != 'INS':
                    if not insertion_only:
                        yield _other_variant(fields, info, samples, samples_available, cutoff,
                                             sequence_checks, candidate_sink, event_kinds)
                    continue
                identifier = fields[2]
                size = abs(int(info['SVLEN']))
                if size == 0:
                    continue
                runs = graph_runs(info.get('SEQ', ''), size)
                if sequence_checks is not None:
                    sequence_checks[identifier] = literal_checks(info.get('SEQ', ''), size)
                query = None
                query_priority = 2
                candidates = []
                candidate_priority = {}
                eligible = ins_seen = size_seen = span_seen = False
                observed_sizes, observed_spans = set(), set()
                if len(fields) == 10:
                    offset = 0
                    for sample in samples:
                        stop = fields[9].find('\t', offset)
                        if stop < 0:
                            stop = len(fields[9])
                        if sample in samples_available:
                            eligible = True
                            observations = (parse_sample_query_intervals(
                                sample, fields[9][offset:stop], fields[8], include_alignment=True)
                                if sequence_checks is not None else parse_sample_query_intervals(
                                    sample, fields[9][offset:stop], fields[8]))
                            for observation in observations:
                                if not accepted_contig(sample, observation[1]):
                                    continue
                                ins_seen = True
                                if len(observed_sizes) < 3:
                                    observed_sizes.add(observation[5])
                                if observation[5] != size:
                                    continue
                                size_seen = True
                                span = observation[3] - observation[2]
                                if len(observed_spans) < 3:
                                    observed_spans.add(span)
                                if span != size:
                                    continue
                                span_seen = True
                                candidate = observation[:5]
                                priority = 1
                                if sequence_checks is not None:
                                    cigar, template_offset = observation[6:8]
                                    # FORMAT CIGAR is relative to this row's
                                    # implicit representative. Prefer a full
                                    # forward exact match at template offset 0;
                                    # QUERYCOORD separately orients FASTA bases.
                                    if (cigar in (f'>{size}=', f'>0H{size}=') and
                                            template_offset in ('.', '0')):
                                        priority = 0
                                if candidate_sink is not None or sequence_checks is not None:
                                    if candidate not in candidate_priority:
                                        candidates.append(candidate)
                                    candidate_priority[candidate] = min(
                                        priority, candidate_priority.get(candidate, priority))
                                # The nominated sequence is still verified
                                # against INFO/SEQ during extraction.
                                if priority < query_priority:
                                    query = candidate
                                    query_priority = priority
                        if stop == len(fields[9]):
                            break
                        offset = stop + 1
                if candidates:
                    candidates.sort(key=candidate_priority.__getitem__)
                    if candidate_priority[candidates[0]] == 0:
                        # Preserve representative provenance. If an exact
                        # representative source fails FASTA validation, do not
                        # substitute a member with a non-exact alignment.
                        candidates = [value for value in candidates
                                      if candidate_priority[value] == 0]
                    query = candidates[0]
                if candidate_sink is not None:
                    candidate_sink(identifier, candidates)
                if query is not None:
                    query_reason = 'ok'
                elif not eligible:
                    query_reason = 'sample_missing: no VCF sample is present in query FASTA list'
                elif not ins_seen:
                    query_reason = 'ins_missing: query-list samples have no INS observation'
                elif not size_seen:
                    values = ','.join(map(str, sorted(observed_sizes))) or 'none'
                    query_reason = (f'size_mismatch: no INS observation has SIZE={size}; '
                                    f'observed SIZE={values}')
                elif not span_seen:
                    values = ','.join(map(str, sorted(observed_spans))) or 'none'
                    query_reason = (f'coordinate_mismatch: no SIZE-matched observation has '
                                    f'QUERYCOORD span={size}; '
                                    f'observed span={values}')
                else:
                    query_reason = 'coordinate_mismatch: no template-length query interval is usable'
                pos = int(fields[1])
                yield (identifier, fields[0], pos, int(info.get('END', pos)),
                       size, runs, query, query_reason,
                       fields[6] in ('PASS', '.') and size >= cutoff)
            except (ValueError, KeyError, IndexError) as error:
                raise ValueError(f'{path}:{lineno}: {error}') from error
