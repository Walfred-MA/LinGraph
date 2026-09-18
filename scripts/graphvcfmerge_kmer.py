#!/usr/bin/env python3
"""Alignment-first insertion refinement with sparse fallback and group bridging.

Each size-compatible query first aligns to the longest template. On failure,
KmerMatch selects between the accepted member with the closest size and the
accepted member with the closest coordinate, then one second alignment decides
the merge. Accepted members then probe every later group through its closest
KmerMatch member, merging whole groups on alignment success. Pairs of at most
50 bases use direct comparison for selection.
"""
from __future__ import annotations

import csv
import functools
import math
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import uuid

DEFAULT_KMERMATCH = str(Path(__file__).resolve().with_name("KmerMatch"))
SHORT_SEQUENCE_LIMIT = 50
KMER_SIMILARITY = 0.5
BATCH_BASES = 64 * 1024 * 1024
KMERMATCH_TASK_SECONDS = 2 * 60 * 60


@functools.lru_cache(maxsize=16)
def _preflight_kmermatch(executable, signature):
    """Exercise the exact CLI/report contract used by the merger once."""
    del signature  # part of the cache key; the executable can be replaced.
    with tempfile.TemporaryDirectory(prefix="kmermatch-preflight-") as directory:
        root = Path(directory)
        sequence = "ACGTTGCATGTCGCATGATGCATGAGAGCT" * 2
        query = root / "query.fa"
        reference = root / "reference.fa"
        report = root / "scores.tsv"
        query.write_text(f">m0\n{sequence}\n", encoding="ascii")
        reference.write_text(f">template\n{sequence}\n", encoding="ascii")
        try:
            completed = subprocess.run(
                [
                    executable, "-i", str(query), "-r", str(reference),
                    "-o", str(report), "-m", "0", "--repeat-score",
                    "--paired",
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise RuntimeError(
                f"KmerMatch preflight could not run {executable!r}: {error}. "
                "Rebuild it from src/KmerMatch/KmerMatch.cpp on this host."
            ) from error
        if completed.returncode != 0:
            raise RuntimeError(
                f"KmerMatch preflight failed for {executable!r} with status "
                f"{completed.returncode}: {completed.stdout[-2000:]}. Rebuild "
                "it from src/KmerMatch/KmerMatch.cpp on this host."
            )
        try:
            with report.open(encoding="utf-8") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                expected_header = [
                    "seqAname", "seqBname", "common_kmers",
                    "kmers_only_in_A", "kmers_only_in_B", "variance_A",
                    "variance_B",
                ]
                if reader.fieldnames != expected_header:
                    raise ValueError(
                        f"unexpected report header {reader.fieldnames!r}"
                    )
                rows = list(reader)
            if len(rows) != 1:
                raise ValueError(f"expected one report row, observed {len(rows)}")
            row = rows[0]
            if row["seqAname"] != "m0" or row["seqBname"] != "template":
                raise ValueError("preflight report names do not match the inputs")
            common, variance_a, variance_b = (
                float(row[name])
                for name in ("common_kmers", "variance_A", "variance_B")
            )
            if not all(
                math.isfinite(value) and value >= 0
                for value in (common, variance_a, variance_b)
            ) or variance_a <= 0 or variance_b <= 0:
                raise ValueError("preflight report contains invalid scores")
            cosine = common / math.sqrt(variance_a * variance_b)
            if not math.isclose(cosine, 1.0, rel_tol=1e-5, abs_tol=1e-6):
                raise ValueError(
                    f"identical-sequence cosine is {cosine}, expected 1"
                )
        except (OSError, KeyError, TypeError, ValueError) as error:
            raise RuntimeError(
                f"KmerMatch preflight produced an incompatible report: {error}. "
                "Rebuild it from src/KmerMatch/KmerMatch.cpp on this host."
            ) from error
    return executable


def resolve_kmermatch(value):
    executable = shutil.which(os.fspath(value))
    if executable is None:
        raise FileNotFoundError(f"KmerMatch is missing or not executable: {value}")
    resolved = str(Path(executable).resolve())
    stat = os.stat(resolved)
    signature = (
        int(stat.st_dev), int(stat.st_ino), int(stat.st_size),
        int(stat.st_mtime_ns), int(stat.st_ctime_ns),
    )
    return _preflight_kmermatch(resolved, signature)


def _read_sequence(handle, locator):
    offset, length = locator
    handle.seek(offset)
    blob = handle.read(length)
    if len(blob) != length:
        raise ValueError("Truncated temporary sequence store")
    return blob.decode("ascii")


def _score_task(args):
    executable, store, template, locators, directory = args
    root = Path(directory)
    with (root / "queries.fa").open("w") as out, open(store, "rb") as sequences:
        for index, locator in locators:
            out.write(f">m{index}\n{_read_sequence(sequences, locator)}\n")
    report = root / "scores.tsv"
    with (root / "kmermatch.log").open("w") as log:
        try:
            completed = subprocess.run(
                [executable, "-i", str(root / "queries.fa"), "-r", template,
                 "-o", str(report), "-m", "0", "--repeat-score"],
                stdout=log, stderr=log, timeout=KMERMATCH_TASK_SECONDS,
            )
        except subprocess.TimeoutExpired as error:
            raise TimeoutError(
                f"KmerMatch task exceeded {KMERMATCH_TASK_SECONDS} seconds"
            ) from error
    if completed.returncode:
        raise RuntimeError(f"KmerMatch failed ({completed.returncode}): "
                           + (root / "kmermatch.log").read_text()[-4000:])
    expected = {f"m{index}": index for index, _ in locators}
    seen, passed = set(), []
    with report.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"seqAname", "seqBname", "common_kmers", "variance_A", "variance_B"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Unrecognized KmerMatch report header")
        for row in reader:
            name = row["seqAname"]
            if name not in expected or name in seen or row["seqBname"] != "template":
                raise ValueError("Unexpected or duplicate pair in KmerMatch report")
            seen.add(name)
            common, va, vb = (float(row[k]) for k in
                              ("common_kmers", "variance_A", "variance_B"))
            if not all(math.isfinite(v) and v >= 0 for v in (common, va, vb)):
                raise ValueError("Invalid KmerMatch score")
            # Zero-norm sequences provide no evidence of k-mer similarity.
            score = common / math.sqrt(va * vb) if va and vb else 0.0
            if score >= KMER_SIMILARITY:
                passed.append(expected[name])
    if seen != set(expected):
        raise ValueError("KmerMatch report is incomplete")
    return passed


def _paired_score_task(args):
    """Score explicitly paired query/reference records with KmerMatch."""
    executable, store, pairs, directory = args
    root = Path(directory)
    queries_path = root / "queries.fa"
    references_path = root / "references.fa"
    with (queries_path.open("w") as queries,
          references_path.open("w") as references,
          open(store, "rb") as sequences):
        for slot, (query_index, query_locator, reference_index,
                   reference_locator) in enumerate(pairs):
            queries.write(
                f">q{slot}\n{_read_sequence(sequences, query_locator)}\n"
            )
            references.write(
                f">r{slot}\n{_read_sequence(sequences, reference_locator)}\n"
            )
    report = root / "scores.tsv"
    with (root / "kmermatch.log").open("w") as log:
        try:
            completed = subprocess.run(
                [
                    executable, "-i", str(queries_path),
                    "-r", str(references_path), "-o", str(report),
                    "-m", "0", "--repeat-score", "--paired",
                ],
                stdout=log, stderr=log, timeout=KMERMATCH_TASK_SECONDS,
            )
        except subprocess.TimeoutExpired as error:
            raise TimeoutError(
                f"KmerMatch paired task exceeded "
                f"{KMERMATCH_TASK_SECONDS} seconds"
            ) from error
    if completed.returncode:
        raise RuntimeError(
            f"KmerMatch paired task failed ({completed.returncode}): "
            + (root / "kmermatch.log").read_text()[-4000:]
        )
    expected = {
        f"q{slot}": (f"r{slot}", query_index, reference_index)
        for slot, (query_index, _query_locator, reference_index,
                   _reference_locator) in enumerate(pairs)
    }
    seen, scores = set(), []
    with report.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "seqAname", "seqBname", "common_kmers",
            "variance_A", "variance_B",
        }
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Unrecognized paired KmerMatch report header")
        for row in reader:
            name = row["seqAname"]
            if name not in expected or name in seen:
                raise ValueError(
                    "Unexpected or duplicate pair in paired KmerMatch report"
                )
            reference_name, query_index, reference_index = expected[name]
            if row["seqBname"] != reference_name:
                raise ValueError("Paired KmerMatch report order differs from input")
            seen.add(name)
            common, variance_a, variance_b = (
                float(row[field])
                for field in ("common_kmers", "variance_A", "variance_B")
            )
            if not all(
                math.isfinite(value) and value >= 0
                for value in (common, variance_a, variance_b)
            ):
                raise ValueError("Invalid paired KmerMatch score")
            score = (
                common / math.sqrt(variance_a * variance_b)
                if variance_a and variance_b else 0.0
            )
            scores.append((query_index, reference_index, score))
    if seen != set(expected):
        raise ValueError("Paired KmerMatch report is incomplete")
    return scores


def _direct_task(args):
    """Compare short pairs once, retaining the traceback for final output."""
    import graphvcfmerge as merger
    store, template, locators, cutoff = args
    passed = []
    with open(store, "rb") as sequences:
        for index, locator in locators:
            query = _read_sequence(sequences, locator)
            if not query or not template:
                continue
            if max(len(query), len(template)) > SHORT_SEQUENCE_LIMIT:
                raise ValueError("Direct comparison exceeds short sequence limit")
            body = merger._fast_python_global_cigar(query, template)
            merger._fast_validate_alignment_body(
                body, len(query), len(template),
                context=f"short KmerMatch member {index}",
            )
            distance = sum(int(n) for n, op in re.findall(r"(\d+)([=XID])", body)
                           if op != "=")
            similarity = 1.0 - distance / float(len(query) + len(template))
            if similarity >= cutoff:
                passed.append((index, body))
    return passed


def _direct_pair_task(args):
    """Score explicitly paired short sequences by direct edit distance."""
    import graphvcfmerge as merger
    store, pairs, _cutoff = args
    scores = []
    with open(store, "rb") as sequences:
        for query_index, query_locator, reference_index, reference_locator in pairs:
            query = _read_sequence(sequences, query_locator)
            reference = _read_sequence(sequences, reference_locator)
            if not query or not reference:
                continue
            if max(len(query), len(reference)) > SHORT_SEQUENCE_LIMIT:
                raise ValueError("Direct paired comparison exceeds short sequence limit")
            body = merger._fast_python_global_cigar(query, reference)
            merger._fast_validate_alignment_body(
                body, len(query), len(reference),
                context=f"short fallback member {query_index}",
            )
            distance = sum(
                int(length)
                for length, operation in re.findall(r"(\d+)([=XID])", body)
                if operation != "="
            )
            similarity = 1.0 - distance / float(len(query) + len(reference))
            scores.append((query_index, reference_index, similarity))
    return scores


def _alignment_task(args):
    """Align a batch of queries to one template and retain every result."""
    import graphvcfmerge as merger
    import graphreftovcf as core
    store, template_locator, locators = args
    aligner_cache = {}
    with open(store, "rb") as sequences:
        template = _read_sequence(sequences, template_locator)
        results = []
        for index, locator in locators:
            query = _read_sequence(sequences, locator)
            similarity, body = merger._fast_pair_alignment(
                query, template, aligner_cache, core,
            )
            if body is not None:
                merger._fast_validate_alignment_body(
                    body, len(query), len(template),
                    context=f"template comparison for member {index}",
                )
            results.append((index, similarity, body))
    return results


def _alignment_pair_task(args):
    """Align each query to its selected fallback member once."""
    import graphvcfmerge as merger
    import graphreftovcf as core
    store, pairs = args
    aligner_cache = {}
    with open(store, "rb") as sequences:
        results = []
        for query_index, query_locator, reference_index, reference_locator in pairs:
            query = _read_sequence(sequences, query_locator)
            reference = _read_sequence(sequences, reference_locator)
            similarity, body = merger._fast_pair_alignment(
                query, reference, aligner_cache, core,
            )
            if body is not None:
                merger._fast_validate_alignment_body(
                    body, len(query), len(reference),
                    context=(
                        f"fallback comparison for member {query_index} "
                        f"against member {reference_index}"
                    ),
                )
            results.append((query_index, reference_index, similarity, body))
    return results


def mappy_body_once(query, template, aligner):
    """One map() invocation, forward primary hit only; no exact rescue.

    H pads unmatched reference ends, I preserves unmatched query ends, as in
    the merger's existing representation. A failed mapping is kept separate.
    """
    if not query or not template:
        return None
    if query.upper() == template.upper():
        return f"{len(query)}="  # Identity needs no mapper invocation.
    best = None
    for hit in aligner.map(query):
        if hit.is_primary and hit.strand == 1 and hit.cigar_str and (
            best is None or hit.mlen > best.mlen
        ):
            best = hit
    if best is None:
        return None
    if not (0 <= best.r_st < best.r_en <= len(template)
            and 0 <= best.q_st < best.q_en <= len(query)):
        raise ValueError("mappy returned invalid alignment coordinates")
    return ((f"{best.r_st}H" if best.r_st else "")
            + (f"{best.q_st}I" if best.q_st else "") + best.cigar_str
            + (f"{len(query) - best.q_en}I" if best.q_en < len(query) else "")
            + (f"{len(template) - best.r_en}H" if best.r_en < len(template) else ""))


def _align_task(args):
    import graphvcfmerge as merger
    import graphreftovcf as core
    store, template_locator, locators = args
    with open(store, "rb") as sequences:
        template = _read_sequence(sequences, template_locator)
        aligner = merger._make_single_thread_insertion_aligner(core, template)
        if not aligner:
            raise RuntimeError("mappy could not build the representative index")
        results = []
        for index, locator in locators:
            query = _read_sequence(sequences, locator)
            body = mappy_body_once(query, template, aligner)
            if body is not None:
                merger._fast_validate_alignment_body(
                    body, len(query), len(template),
                    context=f"KmerMatch member {index}",
                )
            results.append((index, body))
        return results


def _batches(indexes, locators, workers):
    """Fine enough for all workers, bounded by bytes and number of members."""
    budget = min(BATCH_BASES, max(1, sum(locators[i][1] for i in indexes)
                                  // (max(1, workers) * 4)))
    batch, bases = [], 0
    for index in indexes:
        length = locators[index][1]
        if batch and (bases + length > budget or len(batch) >= 1024):
            yield batch
            batch, bases = [], 0
        batch.append((index, locators[index]))
        bases += length
    if batch:
        yield batch


def _pair_batches(pairs, locators, workers):
    """Bound explicit pair tasks by sequence bytes and pair count."""
    total = sum(
        locators[query][1] + locators[reference][1]
        for query, reference in pairs
    )
    budget = min(BATCH_BASES, max(1, total // (max(1, workers) * 4)))
    batch, bases = [], 0
    for query, reference in pairs:
        pair_bases = locators[query][1] + locators[reference][1]
        if batch and (bases + pair_bases > budget or len(batch) >= 1024):
            yield batch
            batch, bases = [], 0
        batch.append((
            query, locators[query], reference, locators[reference],
        ))
        bases += pair_bases
    if batch:
        yield batch


def _fallback_candidates(query, chosen, members, size_threshold, template):
    """Closest-size and closest-coordinate accepted members for one failure."""
    eligible = [
        candidate for candidate in chosen
        if candidate != template
        and min(members[query][2], members[candidate][2])
        / float(max(1, max(members[query][2], members[candidate][2])))
        >= size_threshold
    ]
    if not eligible:
        return []
    size_pick = max(eligible, key=lambda candidate: (
        min(members[query][2], members[candidate][2])
        / float(max(1, max(members[query][2], members[candidate][2]))),
        -abs(members[query][0] - members[candidate][0]),
        -candidate,
    ))
    distance_pick = min(eligible, key=lambda candidate: (
        abs(members[query][0] - members[candidate][0]),
        -min(members[query][2], members[candidate][2])
        / float(max(1, max(members[query][2], members[candidate][2]))),
        candidate,
    ))
    return list(dict.fromkeys((size_pick, distance_pick)))


def _serial_map(function, tasks):
    return [function(task) for task in tasks]


def _closest_group_members(query, groups, members, locators, store, root,
                           executable, size_threshold, sequence_threshold,
                           map_batches, workers):
    """Choose one member per group, applying the size gate before scoring."""
    import graphvcfmerge as merger
    if not locators[query][1]:
        return {}
    candidates = {
        candidate: slot for slot, indexes in groups for candidate in indexes
        if locators[candidate][1] and merger._size_similarity(
            members[query][2], members[candidate][2],
        ) >= size_threshold
    }
    direct_pairs, kmer_pairs = [], []
    for candidate in candidates:
        destination = (
            direct_pairs
            if max(locators[query][1], locators[candidate][1])
            <= SHORT_SEQUENCE_LIMIT else kmer_pairs
        )
        destination.append((query, candidate))
    best = {}

    def record_best(scores, cutoff):
        for result_query, candidate, score in scores:
            if result_query != query or candidate not in candidates:
                raise ValueError("Unexpected group-bridge score pair")
            if score < cutoff:
                continue
            slot = candidates[candidate]
            # Stable ties do not depend on worker completion order.
            rank = (score, -candidate)
            if slot not in best or rank > best[slot][0]:
                best[slot] = (rank, candidate)

    direct_tasks = (
        (store, batch, sequence_threshold)
        for batch in _pair_batches(direct_pairs, locators, workers)
    )
    for scores in map_batches(_direct_pair_task, direct_tasks):
        record_best(scores, sequence_threshold)

    pending = iter(_pair_batches(kmer_pairs, locators, workers))
    while True:
        tasks, directories = [], []
        try:
            for _ in range(max(1, workers) * 2):
                batch = next(pending, None)
                if batch is None:
                    break
                directory = tempfile.mkdtemp(prefix="bridge-score-", dir=root)
                directories.append(directory)
                tasks.append((executable, store, batch, directory))
            if not tasks:
                break
            for scores in map_batches(_paired_score_task, tasks):
                record_best(scores, KMER_SIMILARITY)
        finally:
            for directory in directories:
                shutil.rmtree(directory)
    return {slot: value[1] for slot, value in sorted(best.items())}


class OnlineGroupMatcher:
    """Lazy disk-backed KmerMatch selection for nested online groups."""

    def __init__(self, members, sequences, context):
        self.members = members
        self.context = context
        self.executable = resolve_kmermatch(context.get("kmermatch", DEFAULT_KMERMATCH))
        self.temporary = tempfile.TemporaryDirectory(
            prefix="merge-online-kmer-", dir=context.get("records_dir"),
        )
        self.root = Path(self.temporary.name)
        self.store = str(self.root / "sequences.bin")
        self.locators = []
        try:
            with open(self.store, "wb") as handle:
                for sequence in sequences:
                    blob = sequence.encode("ascii")
                    self.locators.append((handle.tell(), len(blob)))
                    handle.write(blob)
            if len(self.locators) != len(members):
                raise ValueError("Sequence/member count mismatch")
        except BaseException:
            self.close()
            raise

    def closest(self, query, groups, mapped=None, workers=1):
        return _closest_group_members(
            query, groups, self.members, self.locators, self.store, self.root,
            self.executable, float(self.context.get("size_similarity", 0.7)),
            float(self.context.get("sequence_similarity", 0.7)), mapped or _serial_map, workers,
        )

    def close(self):
        self.temporary.cleanup()


def _bridge_later_groups(partitions, members, locators, store, root, executable,
                         size_threshold, sequence_threshold, map_batches,
                         workers):
    """Use each accepted query to test one closest member per later group.

    Initial acceptance order is retained. A group absorbed by an earlier group
    follows that group's ownership, including for its still-pending queries.
    Each query scans later surviving groups once; failed earlier groups are not
    revisited. Only moved members need fresh final-template CIGARs.
    """
    owners = list(range(len(partitions)))
    groups = [list(indexes) for _pick, indexes, _bodies in partitions]
    active = set(range(len(partitions)))
    pass_id = f"{os.getpid()}:{uuid.uuid4().hex[:12]}"
    initial_groups = len(active)
    bridges = [
        (slot, query)
        for slot, (pick, indexes, _bodies) in enumerate(partitions)
        for query in indexes if query != pick
    ]
    print(
        f"[merge:bridge:start] pass={pass_id} groups={initial_groups} "
        f"members={sum(len(group) for group in groups)} queries={len(bridges)}",
        file=sys.stderr, flush=True,
    )

    def owner(slot):
        while owners[slot] != slot:
            owners[slot] = owners[owners[slot]]
            slot = owners[slot]
        return slot

    for original_slot, query in bridges:
        if len(active) <= 1:
            break
        first = owner(original_slot)
        if first not in active:
            raise RuntimeError(f"bridge pass {pass_id}: inactive query owner {first}")
        best = _closest_group_members(
            query,
            [(slot, groups[slot]) for slot in sorted(active) if slot > first],
            members, locators, store, root, executable, size_threshold,
            sequence_threshold, map_batches, workers,
        )
        candidates = {candidate: slot for slot, candidate in best.items()}
        selected = [(query, candidate) for candidate in best.values()]
        accepted = set()
        alignment_tasks = (
            (store, batch)
            for batch in _pair_batches(selected, locators, workers)
        )
        for results in map_batches(_alignment_pair_task, alignment_tasks):
            for result_query, candidate, similarity, body in results:
                if result_query != query or candidate not in candidates:
                    raise ValueError("Unexpected group-bridge alignment pair")
                if similarity >= sequence_threshold and body is not None:
                    accepted.add(candidates[candidate])
        # Validate the whole result before changing group ownership. Never
        # silently count or absorb a group that was already merged.
        for slot in accepted:
            if slot <= first or slot not in active or owner(slot) != slot or not groups[slot]:
                raise RuntimeError(
                    f"bridge pass {pass_id}: attempted to remerge or absorb "
                    f"invalid group {slot} into {first}"
                )
        before = len(active)
        added_members = sum(len(groups[slot]) for slot in accepted)
        for slot in sorted(accepted):
            owners[slot] = first
            groups[first].extend(groups[slot])
            groups[slot] = []
            active.remove(slot)
        if accepted:
            print(
                f"[merge:bridge] member {members[query][4]}: "
                f"merged {len(accepted)} later group(s) into template "
                f"{members[partitions[first][0]][4]}; pass={pass_id} "
                f"groups={before}->{len(active)} "
                f"absorbed_total={initial_groups - len(active)} "
                f"added_members={added_members} template_members={len(groups[first])}",
                file=sys.stderr, flush=True,
            )

    print(
        f"[merge:bridge:end] pass={pass_id} groups={initial_groups}->{len(active)} "
        f"absorbed_total={initial_groups - len(active)}; remapping moved members",
        file=sys.stderr, flush=True,
    )
    merged = []
    for slot, (pick, indexes, original_bodies) in enumerate(partitions):
        if owner(slot) != slot:
            continue
        bodies = {index: original_bodies[index] for index in indexes}
        moved = [index for index in groups[slot] if index not in bodies]
        tasks = (
            (store, locators[pick], batch)
            for batch in _batches(moved, locators, workers)
        )
        for results in map_batches(_alignment_task, tasks):
            for index, _similarity, body in results:
                if body is None:
                    # Preserve the full query when it is linked transitively
                    # but has no direct alignment to the surviving template.
                    body = (
                        (f"{locators[index][1]}I" if locators[index][1] else "")
                        + (f"{locators[pick][1]}H" if locators[pick][1] else "")
                    )
                bodies[index] = body
        merged.append((pick, groups[slot], bodies))
    return merged


def refine_loaded(members, sequences, row_ids, context, pooled_map=None, workers=1):
    """Alignment-first refinement, sparse fallback, then later-group bridging."""
    import graphvcfmerge as merger
    mapped = pooled_map or _serial_map
    executable = resolve_kmermatch(context.get("kmermatch", DEFAULT_KMERMATCH))
    size_threshold = float(context.get("size_similarity", 0.7))
    sequence_threshold = float(context.get("sequence_similarity", 0.7))
    with tempfile.TemporaryDirectory(prefix="merge-kmer-") as directory:
        root = Path(directory)
        store = str(root / "sequences.bin")
        locators = []
        with open(store, "wb") as out:
            for sequence in sequences:
                blob = sequence.encode("ascii")
                locators.append((out.tell(), len(blob)))
                out.write(blob)
        if len(locators) != len(members):
            raise ValueError("Sequence/member count mismatch")
        # Actual extracted length determines longest; stable ties preserve reproducibility.
        selection = sorted(range(len(members)), key=lambda i: (
            locators[i][1], -members[i][0], row_ids.get(members[i][4], ""), -members[i][4],
        ), reverse=True)
        remaining = set(selection)
        partitions = []
        score_number = 0
        alignment_counts = {index: 0 for index in selection}

        def map_batches(function, tasks):
            pending = iter(tasks)
            while True:
                window = []
                for _ in range(max(1, workers) * 2):
                    task = next(pending, None)
                    if task is None:
                        break
                    window.append(task)
                if not window:
                    return
                for result in mapped(function, window):
                    yield result

        for pick in selection:
            if pick not in remaining:
                continue
            remaining.remove(pick)
            representative = members[pick]
            chosen = [pick]
            template_bodies = {pick: None}
            initial_attempted = set()

            # First path: align every size-compatible member directly to the
            # longest template. Keep the body even if its score does not pass.
            initial = [
                index for index in selection
                if index in remaining and locators[index][1]
                and alignment_counts[index] < 2
                and merger._size_similarity(
                    representative[2], members[index][2],
                ) >= size_threshold
            ]
            initial_tasks = (
                (store, locators[pick], batch)
                for batch in _batches(initial, locators, workers)
            )
            for index in initial:
                alignment_counts[index] += 1
            initial_passed = []
            for results in map_batches(_alignment_task, initial_tasks):
                for index, similarity, body in results:
                    initial_attempted.add(index)
                    template_bodies[index] = body
                    if similarity >= sequence_threshold and body is not None:
                        initial_passed.append(index)
            chosen.extend(initial_passed)
            remaining.difference_update(initial_passed)

            # Second path: only initial failures are eligible. KmerMatch picks
            # between the closest-size and closest-coordinate accepted member.
            candidate_pairs = []
            pair_order = {}
            for query in selection:
                if (query not in remaining or not locators[query][1]
                        or alignment_counts[query] >= 2):
                    continue
                candidates = _fallback_candidates(
                    query, chosen, members, size_threshold, pick,
                )
                for order, candidate in enumerate(candidates):
                    pair_order[(query, candidate)] = order
                    candidate_pairs.append((query, candidate))

            direct_pairs = [
                pair for pair in candidate_pairs
                if max(locators[pair[0]][1], locators[pair[1]][1])
                <= SHORT_SEQUENCE_LIMIT
            ]
            direct_pair_set = set(direct_pairs)
            kmer_pairs = [
                pair for pair in candidate_pairs if pair not in direct_pair_set
            ]
            scores = []
            direct_tasks = (
                (store, batch, sequence_threshold)
                for batch in _pair_batches(direct_pairs, locators, workers)
            )
            for result in map_batches(_direct_pair_task, direct_tasks):
                scores.extend(result)

            kmer_batches = iter(_pair_batches(kmer_pairs, locators, workers))
            while True:
                tasks, task_roots = [], []
                for _ in range(max(1, workers) * 2):
                    batch = next(kmer_batches, None)
                    if batch is None:
                        break
                    task_root = root / f"fallback_score_{score_number}"
                    score_number += 1
                    task_root.mkdir()
                    task_roots.append(task_root)
                    tasks.append((
                        executable, store, batch, str(task_root),
                    ))
                if not tasks:
                    break
                try:
                    for result in mapped(_paired_score_task, tasks):
                        scores.extend(result)
                finally:
                    for task_root in task_roots:
                        shutil.rmtree(task_root)

            best = {}
            for query, candidate, score in scores:
                cutoff = (
                    sequence_threshold
                    if (query, candidate) in direct_pair_set
                    else KMER_SIMILARITY
                )
                if score < cutoff:
                    continue
                rank = (score, -pair_order[(query, candidate)], -candidate)
                if query not in best or rank > best[query][0]:
                    best[query] = (rank, candidate)

            selected_pairs = [
                (query, value[1]) for query, value in best.items()
                if alignment_counts[query] < 2
            ]
            for query, _candidate in selected_pairs:
                alignment_counts[query] += 1
            fallback_passed = []
            fallback_tasks = (
                (store, batch)
                for batch in _pair_batches(selected_pairs, locators, workers)
            )
            for results in map_batches(_alignment_pair_task, fallback_tasks):
                for query, _candidate, similarity, body in results:
                    if similarity >= sequence_threshold and body is not None:
                        fallback_passed.append(query)

            # A member skipped against the template because of the 0.7 size
            # gate needs one template alignment for its emitted VCF CIGAR.
            makeup = [
                index for index in fallback_passed
                if index not in initial_attempted
                and alignment_counts[index] < 2
            ]
            for index in makeup:
                alignment_counts[index] += 1
            makeup_tasks = (
                (store, locators[pick], batch)
                for batch in _batches(makeup, locators, workers)
            )
            for results in map_batches(_alignment_task, makeup_tasks):
                for index, _similarity, body in results:
                    template_bodies[index] = body

            for index in fallback_passed:
                if template_bodies.get(index) is None:
                    with open(store, "rb") as sequence_store:
                        query = _read_sequence(sequence_store, locators[index])
                        template = _read_sequence(sequence_store, locators[pick])
                    template_bodies[index] = merger._pure_insertion_alignment_body(
                        query, template,
                    )
            chosen.extend(fallback_passed)
            remaining.difference_update(fallback_passed)
            partitions.append((pick, chosen, template_bodies))
            if len(members) > 1:
                print(
                    f"[merge:kmer] template {representative[4]}: "
                    f"{len(initial_passed)} direct alignment pass, "
                    f"{len(fallback_passed)} fallback pass; "
                    f"{len(remaining)} unassigned",
                    file=sys.stderr, flush=True,
                )

        if len(partitions) > 1:
            partitions = _bridge_later_groups(
                partitions, members, locators, store, root, executable,
                size_threshold, sequence_threshold, map_batches, workers,
            )

        refined = []
        for pick, indexes, bodies in partitions:
            group = [members[index] + (bodies[index],) for index in indexes]
            group.sort(key=lambda member: (
                member[0], member[1], member[2], member[4],
            ))
            refined.append((members[pick], group))
        refined.sort(key=lambda group: (min(m[0] for m in group[1]), min(m[4] for m in group[1])))
        return refined


def refine_kmer_first(members, sequences, row_ids, context,
                      pooled_map=None, workers=1):
    """Timeout fallback: Kmer/direct partitions, template alignment, bridging."""
    import graphvcfmerge as merger
    mapped = pooled_map or _serial_map
    executable = resolve_kmermatch(context.get("kmermatch", DEFAULT_KMERMATCH))
    size_threshold = float(context.get("size_similarity", 0.7))
    sequence_threshold = float(context.get("sequence_similarity", 0.7))
    with tempfile.TemporaryDirectory(prefix="merge-kmer-timeout-") as directory:
        root = Path(directory)
        store = str(root / "sequences.bin")
        locators = []
        with open(store, "wb") as out:
            for sequence in sequences:
                blob = sequence.encode("ascii")
                locators.append((out.tell(), len(blob)))
                out.write(blob)
        if len(locators) != len(members):
            raise ValueError("Sequence/member count mismatch")
        selection = sorted(range(len(members)), key=lambda index: (
            locators[index][1], -members[index][0],
            row_ids.get(members[index][4], ""), -members[index][4],
        ), reverse=True)
        remaining = set(selection)
        partitions = []
        score_number = 0

        def map_batches(function, tasks):
            pending = iter(tasks)
            while True:
                window = []
                for _ in range(max(1, workers) * 2):
                    task = next(pending, None)
                    if task is None:
                        break
                    window.append(task)
                if not window:
                    return
                for result in mapped(function, window):
                    yield result

        for pick in selection:
            if pick not in remaining:
                continue
            remaining.remove(pick)
            eligible = [
                index for index in selection
                if index in remaining and locators[index][1]
                and merger._size_similarity(
                    members[pick][2], members[index][2],
                ) >= size_threshold
            ]
            pairs = [(index, pick) for index in eligible]
            direct_pairs = [
                pair for pair in pairs
                if max(locators[pair[0]][1], locators[pair[1]][1])
                <= SHORT_SEQUENCE_LIMIT
            ]
            direct_pair_set = set(direct_pairs)
            kmer_pairs = [pair for pair in pairs if pair not in direct_pair_set]
            accepted = set()
            direct_tasks = (
                (store, batch, sequence_threshold)
                for batch in _pair_batches(direct_pairs, locators, workers)
            )
            for scores in map_batches(_direct_pair_task, direct_tasks):
                accepted.update(
                    query for query, _reference, score in scores
                    if score >= sequence_threshold
                )
            kmer_batches = iter(_pair_batches(kmer_pairs, locators, workers))
            while True:
                tasks, task_roots = [], []
                for _ in range(max(1, workers) * 2):
                    batch = next(kmer_batches, None)
                    if batch is None:
                        break
                    task_root = root / f"score_{score_number}"
                    score_number += 1
                    task_root.mkdir()
                    task_roots.append(task_root)
                    tasks.append((executable, store, batch, str(task_root)))
                if not tasks:
                    break
                try:
                    for scores in mapped(_paired_score_task, tasks):
                        accepted.update(
                            query for query, _reference, score in scores
                            if score >= KMER_SIMILARITY
                        )
                finally:
                    for task_root in task_roots:
                        shutil.rmtree(task_root)
            chosen = [pick] + [index for index in eligible if index in accepted]
            remaining.difference_update(accepted)
            partitions.append((pick, chosen))
            print(
                f"[merge:kmer-timeout] template {members[pick][4]}: "
                f"{len(chosen)} members; {len(remaining)} unassigned",
                file=sys.stderr, flush=True,
            )

        aligned_partitions = []
        for pick, indexes in partitions:
            bodies = {pick: None}
            candidates = [index for index in indexes if index != pick]
            tasks = (
                (store, locators[pick], batch)
                for batch in _batches(candidates, locators, workers)
            )
            for results in map_batches(_alignment_task, tasks):
                for index, _similarity, body in results:
                    bodies[index] = body
            accepted_indexes = [pick]
            for index in candidates:
                body = bodies.get(index)
                if body is None:
                    aligned_partitions.append((index, [index], {index: None}))
                else:
                    accepted_indexes.append(index)
            aligned_partitions.append((pick, accepted_indexes, bodies))

        selection_rank = {index: rank for rank, index in enumerate(selection)}
        aligned_partitions.sort(key=lambda group: selection_rank[group[0]])
        if len(aligned_partitions) > 1:
            aligned_partitions = _bridge_later_groups(
                aligned_partitions, members, locators, store, root, executable,
                size_threshold, sequence_threshold, map_batches, workers,
            )
        refined = []
        for pick, indexes, bodies in aligned_partitions:
            accepted_group = [members[index] + (bodies[index],) for index in indexes]
            accepted_group.sort(key=lambda member: (
                member[0], member[1], member[2], member[4],
            ))
            refined.append((members[pick], accepted_group))
        refined.sort(key=lambda group: (
            min(member[0] for member in group[1]),
            min(member[4] for member in group[1]),
        ))
        return refined


def refine_cluster(members, context, pooled_map=None, workers=1):
    import graphvcfmerge as merger
    merger._init_fast_worker(context)
    members = sorted(members, key=lambda member: member[3])
    row_ids = {}

    def sequences():
        observations = merger._fast_iter_observations([member[3] for member in members])
        for member, (ref, observation) in zip(members, observations):
            if ref != member[3]:
                raise ValueError("Sequence extraction order differs from member order")
            row_ids[member[4]] = observation.row_id
            yield observation.sequence

    return refine_loaded(members, sequences(), row_ids, context, pooled_map, workers)
