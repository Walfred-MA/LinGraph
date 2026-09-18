#!/usr/bin/env python3
"""Prepare, extract, and assemble inputs for cohort graph calling."""

from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import dataclasses
import json
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import time
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple


HERE = Path(__file__).resolve()
REPOSITORY = HERE.parents[3]
sys.path.insert(0, str(REPOSITORY))

from match_partition_blocks import cohort_sample_for_query
from summarize_partition_hotspot_segments import read_alignment_output
from assembly_contigs import HG38_MAIN_CONTIGS, accepted_contig


REFMATCH_HEADER = (
    "query\tmatch\tquery_locus\tmatch_locus\t"
    "ensembl\tdistance\ttag\tmatches\n"
)
LOCATION_PATTERN = re.compile(r"^(.+):(-?\d+)-(-?\d+)([+-])?$")


def parse_loc(value: str) -> Tuple[str, int, int, str]:
    """Parse a GenomeLift locus without importing the full GenomeLift CLI."""
    match = LOCATION_PATTERN.match(value.strip())
    if not match:
        raise ValueError(f"Malformed location: {value}")
    contig, start_text, end_text, strand = match.groups()
    start, end = int(start_text), int(end_text)
    if end < start:
        start, end = end, start
    return contig, start, end, strand or "+"


def parse_allele_name(value: str) -> Tuple[str, str, str, str, str]:
    """Parse an allele name using the same convention as GenomeLift."""
    parts = value.split("_")
    if len(parts) < 4:
        raise ValueError(f"Unexpected allelename format: {value}")
    group_name = "_".join(parts[:-3])
    sample_name, haplotype, index = parts[-3:]
    return (
        group_name,
        sample_name,
        haplotype,
        index,
        f"{sample_name}_{haplotype}",
    )


def load_config(path: str) -> dict:
    with open(path, "rt", encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("protocol") != "cohort-call-v1":
        raise ValueError(f"unsupported cohort calling protocol in {path}")
    return config


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, path)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def extraction_batches(config: dict) -> List[Tuple[str, List[str]]]:
    size = int(config["extract_batch_size"])
    partitions = list(config["partitions"])
    return [
        (f"batch{index // size + 1:06d}", partitions[index:index + size])
        for index in range(0, len(partitions), size)
    ]


def match_batches(config: dict) -> List[Tuple[str, List[str]]]:
    """Split partitions into a bounded number of balanced contiguous batches."""
    partitions = list(config["partitions"])
    count = min(int(config["match_batches"]), len(partitions))
    base, remainder = divmod(len(partitions), count)
    output = []
    start = 0
    for index in range(count):
        size = base + (1 if index < remainder else 0)
        output.append((
            f"batch{index + 1:06d}",
            partitions[start:start + size],
        ))
        start += size
    return output


def prepare_layout(config: dict) -> None:
    graph = Path(config["graph_folder"])
    output = Path(config["output_folder"])
    source_graphs = graph / "Graphs"
    call_graphs = output / "Graphs"
    samples_root = output / "samples"
    inputs = output / "inputs"
    checkpoints = output / "checkpoints"

    missing: List[str] = []
    for partition in config["partitions"]:
        directory = source_graphs / partition
        for suffix in (
            ".FA", "_align.txt", "_breakpoint_consistency.bed", "cache.json",
        ):
            path = directory / f"{partition}{suffix}"
            if not path.is_file():
                missing.append(str(path))
                if len(missing) >= 20:
                    break
        if len(missing) >= 20:
            break
    if missing:
        raise FileNotFoundError(
            "active partitions are missing completed graph inputs: "
            + ", ".join(missing)
        )

    call_graphs.mkdir(parents=True, exist_ok=True)
    samples_root.mkdir(parents=True, exist_ok=True)
    (output / "logs").mkdir(parents=True, exist_ok=True)
    for partition in config["partitions"]:
        (call_graphs / partition).mkdir(parents=True, exist_ok=True)
    for sample in config["selected_samples"]:
        (samples_root / sample / "pipeline_logs").mkdir(
            parents=True, exist_ok=True,
        )

    normalized_rows = []
    by_name = {}
    for assembly in config["assemblies"]:
        by_name[assembly["name"]] = assembly
        normalized_rows.append(
            f"{assembly['name']}\t{assembly['fasta']}\t{assembly['fai']}\n"
        )
    atomic_write(inputs / "query_paths.normalized.txt", "".join(normalized_rows))
    selected_rows = [
        f"{sample}\t{by_name[sample]['fasta']}\t{by_name[sample]['fai']}\n"
        for sample in config["selected_samples"]
    ]
    atomic_write(inputs / "selected_query_paths.txt", "".join(selected_rows))
    atomic_write(
        inputs / "selected_samples.list",
        "".join(f"{sample}\n" for sample in config["selected_samples"]),
    )
    atomic_write(
        checkpoints / "Graphs.list",
        "".join(f"{partition}\n" for partition in config["partitions"]),
    )

    match_root = inputs / "match_batches"
    match_root.mkdir(parents=True, exist_ok=True)
    for batch, partitions in match_batches(config):
        atomic_write(
            match_root / f"{batch}.list",
            "".join(f"{partition}\n" for partition in partitions),
        )

    extraction_root = inputs / "extraction_batches"
    extraction_root.mkdir(parents=True, exist_ok=True)
    for batch, partitions in extraction_batches(config):
        atomic_write(
            extraction_root / f"{batch}.list",
            "".join(f"{partition}\n" for partition in partitions),
        )

    marker = {
        "protocol": config["protocol"],
        "reference": config["reference"],
        "samples": config["selected_samples"],
        "sample_count": len(config["selected_samples"]),
        "partition_count": len(config["partitions"]),
        "extraction_batch_size": config["extract_batch_size"],
        "extraction_batches": len(extraction_batches(config)),
        "match_batches": len(match_batches(config)),
        "hg38_main_contigs": sorted(HG38_MAIN_CONTIGS),
    }
    atomic_write(
        checkpoints / "layout.json",
        json.dumps(marker, indent=2, sort_keys=True) + "\n",
    )


def read_names(path: str) -> List[str]:
    with open(path, "rt", encoding="utf-8") as handle:
        return [
            raw.strip().split()[0] for raw in handle
            if raw.strip() and not raw.lstrip().startswith("#")
        ]


def allele_genome(name: str) -> Optional[str]:
    try:
        return parse_allele_name(name)[4]
    except (TypeError, ValueError):
        return None


def locus_contig(value: str) -> Optional[str]:
    try:
        return parse_loc(value)[0]
    except (TypeError, ValueError):
        return None


def hg38_alt_allele(allele: str, locus: str) -> bool:
    return (
        allele_genome(allele) == "HG38_h1"
        and not accepted_contig('HG38_h1', locus_contig(locus))
    )


def accepted_alignment_contig(sample: str, contig: str) -> bool:
    """Apply the same contig filter to routing and reference detection."""
    return accepted_contig(sample, contig)


def selected_reference_alignment(row, sample: str, reference: str) -> bool:
    """Whether this row will represent the selected reference downstream."""
    return bool(
        sample == reference
        and accepted_alignment_contig(sample, row.contig)
        and row.graph_path
        and row.graph_path != "*"
    )


class BufferedWriters:
    """Buffer many routes while opening at most one output file at a time."""

    def __init__(
        self,
        max_buffer_bytes: int,
        per_path_bytes: int = 1024 * 1024,
    ):
        self.max_buffer_bytes = max(1, int(max_buffer_bytes))
        self.per_path_bytes = max(1, int(per_path_bytes))
        self.total_bytes = 0
        self.buffers: "collections.OrderedDict[Path, Tuple[List[str], int]]" = (
            collections.OrderedDict()
        )

    def _flush(self, path: Path) -> None:
        texts, size = self.buffers.pop(path)
        # Do not retain per-sample descriptors. ProcessPool workers inherit
        # control pipes, so even an RLIMIT-aware handle cache can exhaust the
        # descriptor limit when --io-jobs is large.
        with path.open(
            "at", encoding="utf-8", buffering=1024 * 1024,
        ) as handle:
            handle.writelines(texts)
        self.total_bytes -= size

    def write(self, path: Path, text: str) -> None:
        texts, size = self.buffers.pop(path, ([], 0))
        texts.append(text)
        size += len(text)
        self.total_bytes += len(text)
        self.buffers[path] = (texts, size)
        if size >= self.per_path_bytes:
            self._flush(path)
        while self.total_bytes > self.max_buffer_bytes and self.buffers:
            self._flush(next(iter(self.buffers)))

    def close(self) -> None:
        while self.buffers:
            self._flush(next(iter(self.buffers)))


def route_buffer_bytes(process_count: int) -> int:
    """Bound aggregate route buffers to roughly two GiB per extraction job."""
    process_count = max(1, int(process_count))
    return max(
        8 * 1024 * 1024,
        min(128 * 1024 * 1024, (2 * 1024**3) // process_count),
    )


def _clean_matches_field(text: str, excluded: set) -> str:
    retained = []
    for token in text.replace(",", ";").split(";"):
        token = token.strip()
        if not token:
            continue
        allele = token.rsplit(":", 1)[0] if ":" in token else token
        if allele not in excluded:
            retained.append(token)
    return ";".join(retained)


def _iter_routed_rows(path: Path) -> Iterable[Tuple[int, str]]:
    if not path.is_file():
        return
    with path.open("rt", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, 1):
            ordinal_text, separator, text = raw.partition("\t")
            if not separator or not ordinal_text.isdigit():
                raise ValueError(
                    f"{path}:{line_number}: malformed routed row"
                )
            yield int(ordinal_text), text


def _merge_routed_rows(
    common_path: Path,
    sample_path: Path,
    output: TextIO,
) -> int:
    """Merge two ordinal-sorted routing streams without loading either."""
    common = iter(_iter_routed_rows(common_path))
    private = iter(_iter_routed_rows(sample_path))
    common_row = next(common, None)
    private_row = next(private, None)
    written = 0
    while common_row is not None or private_row is not None:
        if private_row is None or (
            common_row is not None and common_row[0] < private_row[0]
        ):
            _ordinal, text = common_row
            common_row = next(common, None)
        else:
            _ordinal, text = private_row
            private_row = next(private, None)
        output.write(text if text.endswith("\n") else text + "\n")
        written += 1
    return written


def finalize_batch_sample(
    temporary: Path,
    routes: Sequence[Path],
    sample: str,
) -> Tuple[int, int]:
    """Materialize one sample after every partition in the batch was read."""
    route_directories = (
        [routes] if isinstance(routes, Path) else list(routes)
    )
    refmatch = temporary / f"{sample}.refmatch.tsv"
    alignment = temporary / f"{sample}.align.txt"
    with refmatch.open("wt", encoding="utf-8") as output:
        output.write(REFMATCH_HEADER)
        refmatch_rows = 0
        for route_directory in route_directories:
            refmatch_rows += _merge_routed_rows(
                route_directory / "reference.refmatch.route",
                route_directory / f"{sample}.refmatch.route",
                output,
            )
    alignment_rows = 0
    with alignment.open("wt", encoding="utf-8") as output:
        for route_directory in route_directories:
            alignment_route = route_directory / f"{sample}.align.route"
            if alignment_route.is_file():
                with alignment_route.open("rt", encoding="utf-8") as source:
                    for raw in source:
                        if raw.strip():
                            output.write(
                                raw if raw.endswith("\n") else raw + "\n"
                            )
                            alignment_rows += 1
    return refmatch_rows, alignment_rows


def _partition_chunks(
    indexed_partitions: Sequence[Tuple[int, str]], workers: int,
) -> List[List[Tuple[int, str]]]:
    """Return nonempty contiguous chunks while preserving graph-list order."""
    count = min(max(1, int(workers)), len(indexed_partitions))
    base, remainder = divmod(len(indexed_partitions), count)
    chunks = []
    start = 0
    for index in range(count):
        size = base + (1 if index < remainder else 0)
        chunks.append(list(indexed_partitions[start:start + size]))
        start += size
    return chunks


def _extract_partition_chunk(
    config: dict,
    indexed_partitions: Sequence[Tuple[int, str]],
    routes_text: str,
) -> Tuple[int, int, int]:
    """Read one disjoint partition chunk and write private route streams."""
    source_graphs = Path(config["graph_folder"]) / "Graphs"
    call_graphs = Path(config["output_folder"]) / "Graphs"
    samples = list(config["selected_samples"])
    sample_set = set(samples)
    all_samples = [row["name"] for row in config["assemblies"]]
    reference = config["reference"]
    routes = Path(routes_text)
    routes.mkdir(parents=True, exist_ok=False)
    writers = BufferedWriters(
        route_buffer_bytes(int(config.get("io_jobs", 1)))
    )
    route_ordinal = 0
    alignment_rows = 0
    refmatch_rows = 0

    try:
        for global_partition_index, partition in indexed_partitions:
            alignment_path = source_graphs / partition / f"{partition}_align.txt"
            refmatch_path = call_graphs / partition / f"{partition}_refmatch.tsv"
            if not alignment_path.is_file():
                raise FileNotFoundError(alignment_path)

            alignment_buffers: Dict[str, List[str]] = collections.defaultdict(list)
            partition_alignments = read_alignment_output(str(alignment_path))
            has_reference_alignment = False
            for row in partition_alignments:
                sample = cohort_sample_for_query(row.query_name, all_samples)
                if selected_reference_alignment(row, sample, reference):
                    has_reference_alignment = True
                if sample not in sample_set:
                    continue
                if not accepted_alignment_contig(sample, row.contig):
                    continue
                # Partition-local graph building records hotspot index 1 in
                # every _align.txt. Use the complete active graph-list index.
                row = dataclasses.replace(
                    row, hotspot_index=global_partition_index,
                )
                alignment_buffers[sample].append(row.line())
            for sample, texts in alignment_buffers.items():
                writers.write(
                    routes / f"{sample}.align.route", "".join(texts),
                )
                alignment_rows += len(texts)

            if refmatch_path.is_file():
                with refmatch_path.open(
                    "rt", encoding="utf-8", newline="",
                ) as handle:
                    reader = csv.reader(handle, delimiter="\t")
                    rows = [
                        fields for fields in reader
                        if fields and fields[0] != "query"
                    ]
            elif has_reference_alignment:
                # A graph covered by the requested reference requires its
                # completed RefMatch relationship table. Missing it means the
                # matching stage is genuinely incomplete.
                raise FileNotFoundError(refmatch_path)
            else:
                # Any partition can lack the reference selected for this run:
                # an alternative partition may contain it, while an ordinary
                # partition may not. There is then no reference relationship
                # to route. Sample alignments continue downstream, where
                # local_reference_templates.py supplies the graph's
                # reference-class path for direct variant calling.
                rows = []
                sys.stderr.write(
                    "[cohort_call_inputs] selected-reference-absent partition "
                    f"{partition}: no RefMatch table; using its local "
                    "template\n"
                )
            malformed = [fields for fields in rows if len(fields) != 8]
            if malformed:
                raise ValueError(
                    f"{refmatch_path}: expected eight RefMatch columns"
                )
            hg38_alt_alleles = {
                fields[0] for fields in rows
                if hg38_alt_allele(fields[0], fields[2])
            }

            # The ordinal is local to this contiguous chunk. Finalization
            # merges common/private streams inside each chunk, then appends
            # chunks in complete graph-list order.
            refmatch_buffers: Dict[Path, List[str]] = collections.defaultdict(list)
            for fields in rows:
                route_ordinal += 1
                query_sample = allele_genome(fields[0])
                if query_sample not in sample_set:
                    continue
                if fields[0] in hg38_alt_alleles:
                    continue
                if reference == "HG38_h1":
                    if fields[1] in hg38_alt_alleles:
                        continue
                    if hg38_alt_allele(fields[1], fields[3]):
                        continue
                    fields[7] = _clean_matches_field(
                        fields[7], hg38_alt_alleles,
                    )

                route = (
                    routes / "reference.refmatch.route"
                    if query_sample == reference
                    else routes / f"{query_sample}.refmatch.route"
                )
                text = "\t".join(fields) + "\n"
                refmatch_buffers[route].append(f"{route_ordinal}\t{text}")
                refmatch_rows += 1
            for route, texts in refmatch_buffers.items():
                writers.write(route, "".join(texts))
    finally:
        writers.close()
    return len(indexed_partitions), refmatch_rows, alignment_rows


def _run_process_tasks(function, arguments: Sequence[tuple], jobs: int):
    """Run independent tasks in processes and return results in input order."""
    if not arguments:
        return []
    workers = min(max(1, int(jobs)), len(arguments))
    if workers == 1:
        return [function(*values) for values in arguments]
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(function, *values) for values in arguments]
        return [future.result() for future in futures]


def _finalize_batch_sample_task(
    temporary_text: str, route_texts: Sequence[str], sample: str,
) -> Tuple[str, int, int]:
    refmatch_rows, alignment_rows = finalize_batch_sample(
        Path(temporary_text), [Path(value) for value in route_texts], sample,
    )
    return sample, refmatch_rows, alignment_rows


def extract_batch(
    config: dict,
    batch_list: str,
    output_dir_text: str,
    jobs: Optional[int] = None,
) -> None:
    output_dir = Path(output_dir_text).resolve()
    expected_root = (Path(config["output_folder"]) / "shards").resolve()
    if output_dir.parent != expected_root:
        raise ValueError(
            f"refusing shard output outside {expected_root}: {output_dir}"
        )
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(
        prefix=f".{output_dir.name}.", dir=output_dir.parent,
    ))
    samples = list(config["selected_samples"])
    partitions = read_names(batch_list)
    if not partitions:
        raise ValueError(f"{batch_list}: empty extraction batch")
    requested_jobs = int(config.get("io_jobs", 1) if jobs is None else jobs)
    if requested_jobs < 1:
        raise ValueError("extraction jobs must be positive")
    global_partition_index = {
        partition: index
        for index, partition in enumerate(config["partitions"], 1)
    }
    indexed_partitions = [
        (global_partition_index[partition], partition)
        for partition in partitions
    ]
    worker_config = {
        "graph_folder": config["graph_folder"],
        "output_folder": config["output_folder"],
        "selected_samples": samples,
        "assemblies": [
            {"name": row["name"]} for row in config["assemblies"]
        ],
        "reference": config["reference"],
        "io_jobs": requested_jobs,
    }
    chunks = _partition_chunks(indexed_partitions, requested_jobs)
    routes_root = temporary / ".routes"
    routes_root.mkdir()
    route_directories = [
        routes_root / f"chunk{index:06d}"
        for index in range(len(chunks))
    ]
    counts = {
        sample: {"refmatch": 0, "alignment": 0}
        for sample in samples
    }

    try:
        workers = min(requested_jobs, len(chunks))
        sys.stderr.write(
            "[cohort_call_inputs] extracting "
            f"{output_dir.name}: {len(partitions)} partitions with "
            f"{workers} process{'es' if workers != 1 else ''}\n"
        )
        sys.stderr.flush()
        started = time.monotonic()
        chunk_results = _run_process_tasks(
            _extract_partition_chunk,
            [
                (worker_config, chunk, str(route_directory))
                for chunk, route_directory in zip(chunks, route_directories)
            ],
            requested_jobs,
        )
        elapsed = max(1e-6, time.monotonic() - started)
        extracted_partitions = sum(result[0] for result in chunk_results)
        sys.stderr.write(
            "[cohort_call_inputs] extracted "
            f"{output_dir.name}: {extracted_partitions}/"
            f"{len(partitions)} partitions ({extracted_partitions / elapsed:.2f}/s)\n"
        )
        sys.stderr.flush()

        route_texts = [str(path) for path in route_directories]
        finalized = _run_process_tasks(
            _finalize_batch_sample_task,
            [
                (str(temporary), route_texts, sample)
                for sample in samples
            ],
            requested_jobs,
        )
        for sample_number, (sample, refmatch_rows, alignment_rows) in enumerate(
            finalized, 1,
        ):
            counts[sample]["refmatch"] = refmatch_rows
            counts[sample]["alignment"] = alignment_rows
            if (
                sample_number == 1
                or sample_number == len(samples)
                or sample_number % 25 == 0
            ):
                sys.stderr.write(
                    "[cohort_call_inputs] finalized "
                    f"{output_dir.name}: sample {sample_number}/"
                    f"{len(samples)} ({sample})\n"
                )
                sys.stderr.flush()
        shutil.rmtree(routes_root)

        manifest = {
            "batch": output_dir.name,
            "partitions": len(partitions),
            "reference": config["reference"],
            "counts": counts,
            "hg38_alternatives_removed": True,
            "io_processes": workers,
        }
        (temporary / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        if output_dir.exists():
            shutil.rmtree(output_dir)
        os.replace(temporary, output_dir)
    finally:
        if temporary.exists():
            shutil.rmtree(temporary)


def _iter_without_header(path: Path, header_token: str) -> Iterable[str]:
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            if raw.rstrip("\r\n").split("\t", 1)[0] == header_token:
                continue
            yield raw if raw.endswith("\n") else raw + "\n"


def assemble_sample(
    config: dict,
    sample: str,
    batches: Optional[Sequence[str]] = None,
) -> None:
    if sample not in config["selected_samples"]:
        raise ValueError(f"sample {sample!r} is not selected")
    output = Path(config["output_folder"])
    sample_dir = output / "samples" / sample
    sample_dir.mkdir(parents=True, exist_ok=True)
    refmatch = sample_dir / f"{sample}.refmatch.tsv"
    alignment = sample_dir / f"{sample}.align.txt"
    refmatch_tmp = refmatch.with_name(refmatch.name + f".tmp.{os.getpid()}")
    alignment_tmp = alignment.with_name(alignment.name + f".tmp.{os.getpid()}")
    batch_names = (
        list(batches) if batches is not None
        else [batch for batch, _partitions in extraction_batches(config)]
    )
    refmatch_rows = 0
    alignment_rows = 0
    observed_genomes = set()
    try:
        with refmatch_tmp.open("wt", encoding="utf-8") as ref_output, \
                alignment_tmp.open("wt", encoding="utf-8") as align_output:
            ref_output.write(REFMATCH_HEADER)
            for batch in batch_names:
                shard_dir = output / "shards" / batch
                ref_shard = shard_dir / f"{sample}.refmatch.tsv"
                align_shard = shard_dir / f"{sample}.align.txt"
                for raw in _iter_without_header(ref_shard, "query"):
                    fields = raw.rstrip("\r\n").split("\t")
                    if len(fields) != 8:
                        raise ValueError(f"{ref_shard}: malformed RefMatch row")
                    genome = allele_genome(fields[0])
                    if genome is not None:
                        observed_genomes.add(genome)
                    ref_output.write(raw)
                    refmatch_rows += 1
                for raw in _iter_without_header(align_shard, "hotspot_index"):
                    align_output.write(raw)
                    alignment_rows += 1

        required_genomes = {config["reference"], sample}
        if not required_genomes.issubset(observed_genomes):
            raise ValueError(
                f"{sample}: RefMatch extraction lacks "
                + ", ".join(sorted(required_genomes - observed_genomes))
            )
        if refmatch_rows == 0 or alignment_rows == 0:
            raise ValueError(
                f"{sample}: extracted {refmatch_rows} RefMatch row(s) and "
                f"{alignment_rows} alignment row(s)"
            )
        os.replace(refmatch_tmp, refmatch)
        os.replace(alignment_tmp, alignment)
    finally:
        for temporary in (refmatch_tmp, alignment_tmp):
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass


def _assemble_sample_task(
    config: dict, sample: str, batches: Sequence[str],
) -> str:
    assemble_sample(config, sample, batches=batches)
    return sample


def validate_match_batch(
    config: dict, partition_list: str, manifest_path: str,
) -> None:
    expected = read_names(partition_list)
    expected_set = set(expected)
    if len(expected_set) != len(expected):
        raise ValueError(f"{partition_list}: duplicate partition")
    rows = {}
    with open(manifest_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            partition = row.get("partition", "")
            if not partition:
                raise ValueError(f"{manifest_path}: row without partition")
            if partition in rows:
                raise ValueError(
                    f"{manifest_path}: duplicate partition {partition!r}"
                )
            rows[partition] = row
    if set(rows) != expected_set:
        missing = sorted(expected_set - set(rows))
        extra = sorted(set(rows) - expected_set)
        raise ValueError(
            f"{manifest_path}: batch membership mismatch; "
            f"missing={missing[:10]!r}, extra={extra[:10]!r}"
        )
    call_graphs = Path(config["output_folder"]) / "Graphs"
    for partition in expected:
        row = rows[partition]
        if row.get("status") not in {"built", "reused"}:
            raise ValueError(
                f"{partition}: match status is {row.get('status')!r}: "
                f"{row.get('message', '')}"
            )
        output = call_graphs / partition / f"{partition}_refmatch.tsv"
        if not output.is_file():
            # The extraction stage can route this partition directly to a
            # local template only when no retained row maps the selected
            # reference. Use exactly the same predicate here; a filtered
            # HG38 alternate-contig row is not a usable reference alignment.
            alignment = (
                Path(config["graph_folder"]) / "Graphs" / partition
                / f"{partition}_align.txt"
            )
            all_samples = [entry["name"] for entry in config["assemblies"]]
            if not any(
                selected_reference_alignment(
                    alignment_row,
                    cohort_sample_for_query(alignment_row.query_name, all_samples),
                    config["reference"],
                )
                for alignment_row in read_alignment_output(str(alignment))
            ):
                continue
        if not output.is_file() or output.stat().st_size == 0:
            raise FileNotFoundError(
                f"{partition}: missing completed RefMatch output {output}"
            )


def extract_and_assemble_all(
    config: dict, marker_path: str, jobs: Optional[int] = None,
) -> None:
    output = Path(config["output_folder"])
    requested_jobs = int(config.get("io_jobs", 1) if jobs is None else jobs)
    if requested_jobs < 1:
        raise ValueError("extraction jobs must be positive")
    batches = extraction_batches(config)
    for batch_number, (batch, partitions) in enumerate(batches, 1):
        sys.stderr.write(
            "[cohort_call_inputs] starting extraction batch "
            f"{batch_number}/{len(batches)} ({batch}; "
            f"{len(partitions)} partitions)\n"
        )
        sys.stderr.flush()
        batch_list = output / "inputs" / "extraction_batches" / f"{batch}.list"
        extract_batch(
            config,
            str(batch_list),
            str(output / "shards" / batch),
            jobs=requested_jobs,
        )
    samples = list(config["selected_samples"])
    batch_names = [batch for batch, _partitions in batches]
    assembly_config = {
        "output_folder": config["output_folder"],
        "selected_samples": samples,
        "reference": config["reference"],
    }
    sys.stderr.write(
        "[cohort_call_inputs] concatenating "
        f"{len(samples)} samples with "
        f"{min(requested_jobs, len(samples))} "
        f"process{'es' if min(requested_jobs, len(samples)) != 1 else ''}\n"
    )
    sys.stderr.flush()
    assembled = _run_process_tasks(
        _assemble_sample_task,
        [(assembly_config, sample, batch_names) for sample in samples],
        requested_jobs,
    )
    for sample_number, sample in enumerate(assembled, 1):
        sys.stderr.write(
            "[cohort_call_inputs] concatenated sample "
            f"{sample_number}/{len(samples)}: {sample}\n"
        )
        sys.stderr.flush()
        assemble_sample(config, sample)
    marker = {
        "protocol": config["protocol"],
        "reference": config["reference"],
        "samples": config["selected_samples"],
        "partitions": len(config["partitions"]),
        "batches": len(extraction_batches(config)),
        "batch_size": config["extract_batch_size"],
        "io_processes": requested_jobs,
    }
    atomic_write(
        Path(marker_path),
        json.dumps(marker, indent=2, sort_keys=True) + "\n",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--run-config", required=True)
    extract = subparsers.add_parser("extract-batch")
    extract.add_argument("--run-config", required=True)
    extract.add_argument("--partition-list", required=True)
    extract.add_argument("--output-dir", required=True)
    assemble = subparsers.add_parser("assemble-sample")
    assemble.add_argument("--run-config", required=True)
    assemble.add_argument("--sample", required=True)
    validate = subparsers.add_parser("validate-match-batch")
    validate.add_argument("--run-config", required=True)
    validate.add_argument("--partition-list", required=True)
    validate.add_argument("--manifest", required=True)
    extract_all = subparsers.add_parser("extract-all")
    extract_all.add_argument("--run-config", required=True)
    extract_all.add_argument("--marker", required=True)
    extract_all.add_argument(
        "--jobs", type=int,
        help="worker processes; default: io_jobs from the run configuration",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    config = load_config(args.run_config)
    if args.command == "prepare":
        prepare_layout(config)
    elif args.command == "extract-batch":
        extract_batch(config, args.partition_list, args.output_dir)
    elif args.command == "assemble-sample":
        assemble_sample(config, args.sample)
    elif args.command == "validate-match-batch":
        validate_match_batch(config, args.partition_list, args.manifest)
    elif args.command == "extract-all":
        extract_and_assemble_all(config, args.marker, jobs=args.jobs)
    else:
        raise ValueError(f"unknown command {args.command!r}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"cohort_call_inputs.py: error: {error}", file=sys.stderr)
        raise SystemExit(1)
