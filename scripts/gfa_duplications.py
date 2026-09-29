"""Shared reference-PA paths for full-locus duplication insertions.

A ``PACLASS=fulllocusdup`` insertion copies a reference PA. Instead of giving
every copy its own literal sequence, all copies whose aligned source intervals
overlap on one reference path share a single new path: a copy of the union of
those source intervals (``TP:Z:duplication``). Each copy then walks that shared
path where ``INFO/ALTERNATIVECIGAR`` aligns it to the PA. Mismatches and
inserted bases stay literal (bubbles), deleted PA bases and PA bases the copy
does not cover become jumps (deletions), and insertion bases outside every
aligned piece remain literal pieces.

A piece is used only if its bases, rebuilt from the PA sequence plus its CIGAR
payloads, occur in the insertion's literal ``INFO/SEQ`` after the previous
piece. Otherwise that piece stays literal, so no unverified sequence enters
the graph.
"""
from collections import defaultdict
import re

from gfa_query_anchors import interval_chunks, reverse_complement
from minsetref_core import IndexedFasta

_CHUNK = re.compile(r'([<>])([^<>:]+):([^<>]*)')
_OP = re.compile(r'(\d+)([=XIDHS])([A-Za-z]*)')


def _parse_piece(text, roots, source_aliases):
    """Return one aligned piece or None when it cannot be used."""
    match = _CHUNK.fullmatch(text)
    if match is None:
        return None
    direction, target, body = match.groups()
    orientation = '+' if direction == '>' else '-'
    ops = []
    position = 0
    while position < len(body):
        op = _OP.match(body, position)
        if op is None:
            return None
        ops.append((op[2], int(op[1]), op[3]))
        position = op.end()
    lead = ops[0][1] if ops and ops[0][0] == 'H' else 0
    ops = [op for op in ops if op[0] != 'H']
    consumed = sum(n for op, n, _payload in ops if op in '=XD')
    if not any(op == '=' for op, _n, _payload in ops):
        return None
    if target in source_aliases:
        alias = source_aliases[target]
        length = alias.end - alias.start
    elif target in roots and roots[target].kind != 'duplication':
        alias, length = None, roots[target].length
    else:
        return None
    if lead + consumed > length:
        return None
    # Oriented coordinates run along the piece's direction on the target.
    low, high = ((lead, lead + consumed) if orientation == '+' else
                 (length - lead - consumed, length - lead))
    if alias is not None:
        target, low, high, orientation = alias.interval(low, high, orientation)
    return {'target': target, 'orientation': orientation, 'low': low,
            'high': high, 'ops': ops}


def _piece_sequence(piece, oriented):
    """Rebuild the piece's query bases; None if a payload is missing."""
    output = []
    cursor = 0
    for op, n, payload in piece['ops']:
        if op == '=':
            output.append(oriented[cursor:cursor + n])
            cursor += n
        elif op == 'X':
            if len(payload) == 2 * n:
                payload = payload[n:]
            if len(payload) != n:
                return None
            output.append(payload)
            cursor += n
        elif op in 'IS':
            if len(payload) != n:
                return None
            output.append(payload)
        elif op == 'D':
            cursor += n
    return ''.join(output)


def _clusters(intervals):
    """Union overlapping [low, high) intervals; return sorted clusters."""
    merged = []
    for low, high in sorted(intervals):
        if merged and low < merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], high)
        else:
            merged.append([low, high])
    return merged


def add_duplication_paths(events, roots, dup_alignments, source_aliases, run_type,
                          root_type, log):
    """Add shared PA roots and rewrite duplication insertions onto them."""
    if not dup_alignments:
        return {}
    insertions = {identifier for identifier in dup_alignments
                  if identifier in events and events[identifier].kind == 'insertion'}
    plan = plan_duplication_paths(dup_alignments, roots, source_aliases, root_type,
                                  insertions, events)
    apply_duplication_paths(events, plan, run_type)
    log_duplication_plan(plan, len(dup_alignments), log)
    return plan['shared']


def log_duplication_plan(plan, count, log):
    log(f'Full-locus duplications: {count} insertion(s), '
        f'{sum(len(values) for values in plan["shared"].values())} shared PA path(s); '
        + ', '.join(f'{key}={value}' for key, value in sorted(plan['stats'].items())))
    for reason, identifiers in sorted(plan['examples'].items()):
        log(f'  not_found_{reason} examples: {", ".join(identifiers)}')


def plan_duplication_paths(dup_alignments, roots, source_aliases, root_type,
                           insertions=None, event_names=()):
    """Verify pieces, cluster them per target and add the shared PA roots.

    Global step: needs every duplication insertion's CIGAR and SEQ, but not
    the variants themselves. ``insertions`` limits it to IDs that are
    insertion variants; ``event_names`` guards the new root names.
    """
    stats = defaultdict(int)
    parsed = {}
    readers = {}
    try:
        for identifier, (cigar, literal) in dup_alignments.items():
            if insertions is not None and identifier not in insertions:
                continue
            if literal is None:
                stats['no_literal_seq'] += 1
                continue
            pieces = []
            for text in cigar.split('&'):
                piece = _parse_piece(text, roots, source_aliases)
                if piece is None:
                    stats['unusable_piece'] += 1
                    continue
                root = roots[piece['target']]
                key = root.path, root.fai
                if key not in readers:
                    readers[key] = IndexedFasta(root.path, root.fai)
                forward = ''.join(interval_chunks(readers[key], piece['target'],
                                                  piece['low'], piece['high'], '+'))
                oriented = forward if piece['orientation'] == '+' else reverse_complement(forward)
                sequence = _piece_sequence(piece, oriented)
                if sequence is None:
                    stats['unusable_piece'] += 1
                    continue
                piece['sequence'] = sequence.upper()
                pieces.append(piece)
            if pieces:
                parsed[identifier] = (pieces, literal.upper())
    finally:
        for reader in readers.values():
            reader.close()

    # Place pieces first: only verified pieces define the shared PA extent.
    placed = {}
    examples = defaultdict(list)
    for identifier, (pieces, literal) in parsed.items():
        cursor = 0
        kept = []
        for piece in pieces:
            position = literal.find(piece['sequence'], cursor)
            if position < 0:
                stats['sequence_not_found'] += 1
                reason = _not_found_reason(piece['sequence'], literal, cursor)
                stats[f'not_found_{reason}'] += 1
                if len(examples[reason]) < 5:
                    examples[reason].append(identifier)
                continue
            kept.append((piece, position))
            cursor = position + len(piece['sequence'])
        if kept:
            placed[identifier] = (kept, literal)
        else:
            stats['kept_literal'] += 1

    by_target = defaultdict(list)
    for kept, _literal in placed.values():
        for piece, _position in kept:
            by_target[piece['target']].append((piece['low'], piece['high']))
    shared = {}
    next_order = max((root.order for root in roots.values()), default=-1) + 1
    for target, intervals in sorted(by_target.items()):
        base = roots[target]
        for low, high in _clusters(intervals):
            name = f'dup_{target}_{low}_{high}'
            if name in roots or name in event_names:
                raise ValueError(f'duplication path name {name!r} collides with another path')
            roots[name] = root_type(name, base.path, base.fai, high - low, True,
                                    'duplication', next_order, record=target, base=low)
            next_order += 1
            shared.setdefault(target, []).append((low, high, name))

    return {'shared': shared, 'placed': placed, 'stats': stats, 'examples': examples}


def apply_duplication_paths(events, plan, run_type):
    """Rewrite the runs of every planned duplication insertion in ``events``."""
    shared, stats = plan['shared'], plan['stats']

    def locate(piece):
        for low, high, name in shared[piece['target']]:
            if low <= piece['low'] and piece['high'] <= high:
                return name, low, high
        raise AssertionError('piece outside its duplication cluster')

    for identifier, (kept, literal) in plan['placed'].items():
        event = events.get(identifier)
        if event is None:
            continue
        runs = []
        unique = 0
        cursor = 0

        def literal_run(start, end):
            nonlocal unique
            if end > start:
                unique += 1
                runs.append(run_type(start, end, 'I', None, 0, 0, '+', unique))

        for piece, position in kept:
            literal_run(cursor, position)
            name, low, high = locate(piece)
            length = high - low
            query = position
            # Walk in oriented target coordinates from the piece start.
            oriented = 0
            for op, n, _payload in piece['ops']:
                if op in '=X':
                    if piece['orientation'] == '+':
                        start = piece['low'] - low + oriented
                        rstart, rend = start, start + n
                    else:
                        end = piece['high'] - low - oriented
                        # Reverse runs store distances from the right edge.
                        rstart, rend = length - end, length - end + n
                    runs.append(run_type(query, query + n, op, name, rstart, rend,
                                         piece['orientation'], 0))
                    query += n
                    oriented += n
                elif op in 'IS':
                    unique += 1
                    runs.append(run_type(query, query + n, 'I', None, 0, 0, '+', unique))
                    query += n
                elif op == 'D':
                    oriented += n
            cursor = position + len(piece['sequence'])
            stats['pieces_on_shared_path'] += 1
        literal_run(cursor, len(literal))
        if sum(run.qend - run.qstart for run in runs) != event.length:
            raise ValueError(f'{identifier}: duplication runs cover '
                             f'{sum(run.qend - run.qstart for run in runs)} of {event.length} bases')
        event.runs = runs
        stats['insertions_on_shared_path'] += 1


def _not_found_reason(sequence, literal, cursor):
    """Why a piece's rebuilt bases are absent from SEQ (diagnostic only)."""
    if reverse_complement(sequence) in literal:
        return 'reverse_complement'
    if sequence in literal:
        return 'before_previous_piece'
    if len(sequence) > len(literal):
        if literal in sequence:
            return 'seq_inside_piece'
        if reverse_complement(literal) in sequence:
            return 'seq_inside_piece_reverse_complement'
        return 'piece_longer_than_seq'
    return 'no_exact_match'
