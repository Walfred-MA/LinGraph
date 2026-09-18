"""Query-only sequence payloads at graph-CIGAR file boundaries."""

import re


_CHUNK_RE = re.compile(r"([><])([^><]*)")
_OP_RE = re.compile(r"(\d+)([=MXIDHSN])([A-Za-z]*)")


def strip_graph_cigar_payloads(text: str) -> str:
    """Keep operation lengths and path names, dropping embedded sequences."""
    def body(value):
        return _OP_RE.sub(lambda match: match[1] + match[2], value)

    def chunk(match):
        direction, value = match.groups()
        name, separator, cigar = value.partition(":")
        if not separator:
            name, cigar = "", value
        return direction + name + separator + body(cigar)

    return _CHUNK_RE.sub(chunk, text) if text.startswith((">", "<")) else body(text)


def query_only_graph_cigar(text: str) -> str:
    """Normalize X/D payloads without changing operations or coordinates.

    Internal pairwise alignments can carry 2n bases after X: n reference
    bases followed by n query bases. Files exchanged between steps use only
    the query half. D can carry optional reference bases, including only a
    partial sequence after adjacent deletions are coalesced. Discard those
    bases at the file boundary, retaining the authoritative deletion length.
    Already query-only or payload-free operations are kept.
    Target names (which can themselves contain digits and X) are untouched.
    Operation/coordinate validation remains the caller's responsibility;
    an invalid X payload is rejected rather than truncated or padded.
    """
    if "X" not in text and "D" not in text:
        return text

    def operation(match):
        length_text, op, payload = match.groups()
        if op == "D" and payload:
            return f"{length_text}D"
        if op != "X" or not payload:
            return match.group(0)
        length = int(length_text)
        if len(payload) == length:
            return match.group(0)
        if length > 0 and len(payload) == 2 * length:
            return f"{length_text}X{payload[length:]}"
        raise ValueError(
            f"X payload has {len(payload)} bases, expected {length} "
            f"(query-only) or {2 * length} (reference+query)"
        )

    def chunk(match):
        direction, value = match.groups()
        name, separator, body = value.partition(":")
        if not separator:
            name, body = "", value
        return direction + name + separator + _OP_RE.sub(operation, body)

    return _CHUNK_RE.sub(chunk, text)
