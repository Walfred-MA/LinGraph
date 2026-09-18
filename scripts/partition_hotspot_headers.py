#!/usr/bin/env python3
"""Reconstruct local-block FASTA headers from KmerSearcher hotspot tables.

The implementation mirrors KmerBlocking's ``extract_local_blocks.py``:
31-mer coordinate conversion, <30 kb merging, 15 kb anchors, and the optional
5 kb reference-overlap extension.  Only compact numeric arrays are retained.
"""
from __future__ import annotations

import dataclasses
import os
from array import array
from bisect import bisect_left, bisect_right
from collections import defaultdict
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

from minsetref_core import mp_context


@dataclasses.dataclass(frozen=True)
class AssemblySource:
    order: int
    name: str
    fasta: str
    fai: str
    hotspot: str


@dataclasses.dataclass(frozen=True)
class BlockSource:
    index: int
    name: str
    contig: str
    start: int
    end: int


@dataclasses.dataclass
class AssemblyHeaderStore:
    """Compact reconstructed headers owned by exactly one assembly."""

    assembly: AssemblySource
    contigs: List[str]
    partition_id: array
    contig_id: array
    source_start: array
    source_end: array
    sequence_index: array
    by_contig: Dict[int, array]
    hotspot_count: int

    def __len__(self) -> int:
        return len(self.partition_id)

    def partition_range(self, partition_id: int) -> Tuple[int, int]:
        return (
            bisect_left(self.partition_id, partition_id),
            bisect_right(self.partition_id, partition_id),
        )

    def packed_size_bytes(self) -> int:
        columns = (
            self.partition_id,
            self.contig_id,
            self.source_start,
            self.source_end,
            self.sequence_index,
        )
        return (
            sum(column.itemsize * len(column) for column in columns)
            + sum(
                column.itemsize * len(column)
                for column in self.by_contig.values()
            )
        )


class CompactHeaderStore:
    """Independent assembly holders with globally unique source contigs."""

    def __init__(
        self,
        partition_names: Sequence[str],
        sources: Sequence[AssemblySource],
    ) -> None:
        self.partitions = list(partition_names)
        self.partition_ids = {
            name: index for index, name in enumerate(self.partitions)
        }
        self.haplotypes = [source.name for source in sources]
        self.haplotype_ids = {
            name: index for index, name in enumerate(self.haplotypes)
        }
        self.holders: List[Optional[AssemblyHeaderStore]] = [
            None for _source in sources
        ]
        self.contig_owner: Dict[str, int] = {}
        self._header_count = 0

    def __len__(self) -> int:
        return self._header_count

    def packed_size_bytes(self) -> int:
        """Return bytes held by numeric columns/indexes (not Python strings)."""
        return sum(
            holder.packed_size_bytes()
            for holder in self.holders
            if holder is not None
        )

    def release_for_sweep(self) -> None:
        """Discard columns needed for catalog creation but not interval sweeps."""
        self.partition_ids.clear()
        self.haplotype_ids.clear()
        for holder in self.iter_holders():
            holder.contig_id = array("I")
            holder.sequence_index = array("I")

    @staticmethod
    def locator(holder_id: int, local_index: int) -> int:
        if not 0 <= holder_id < 2**32 or not 0 <= local_index < 2**32:
            raise OverflowError("assembly/local header index exceeds uint32")
        return (holder_id << 32) | local_index

    def decode_locator(self, locator: int) -> Tuple[AssemblyHeaderStore, int]:
        holder_id, local_index = divmod(locator, 2**32)
        try:
            holder = self.holders[holder_id]
        except IndexError as error:
            raise IndexError(f"invalid header locator {locator}") from error
        if holder is None or local_index >= len(holder):
            raise IndexError(f"invalid header locator {locator}")
        return holder, local_index

    def add_holder(self, holder: AssemblyHeaderStore) -> None:
        holder_id = holder.assembly.order
        if holder_id >= len(self.holders):
            raise ValueError(
                f"invalid assembly order {holder_id} for {holder.assembly.name}"
            )
        if self.holders[holder_id] is not None:
            raise ValueError(
                f"duplicate reconstructed assembly {holder.assembly.name!r}"
            )
        if len(holder) >= 2**32:
            raise OverflowError(
                f"{holder.assembly.name}: header count exceeds uint32"
            )
        for contig in holder.contigs:
            prior = self.contig_owner.get(contig)
            if prior is not None:
                raise ValueError(
                    f"source contig {contig!r} occurs in both "
                    f"{self.haplotypes[prior]!r} and "
                    f"{holder.assembly.name!r}; contig names must be "
                    "globally unique across assemblies"
                )
            self.contig_owner[contig] = holder_id
        self.holders[holder_id] = holder
        self._header_count += len(holder)

    def iter_holders(self) -> Iterator[AssemblyHeaderStore]:
        for holder in self.holders:
            if holder is None:
                raise RuntimeError("incomplete assembly header store")
            yield holder

    def iter_contigs(
        self,
    ) -> Iterator[Tuple[Tuple[str, str], int, array]]:
        for holder_id, holder in enumerate(self.iter_holders()):
            for contig_id, members in holder.by_contig.items():
                yield (
                    (holder.assembly.name, holder.contigs[contig_id]),
                    holder_id,
                    members,
                )

    def iter_partition_entries(
        self,
        partition: str,
    ) -> Iterator[Tuple[int, AssemblyHeaderStore, int]]:
        try:
            partition_id = self.partition_ids[partition]
        except KeyError as error:
            raise KeyError(f"unknown partition {partition!r}") from error
        for holder_id, holder in enumerate(self.iter_holders()):
            start, end = holder.partition_range(partition_id)
            for local_index in range(start, end):
                yield holder_id, holder, local_index

    def iter_partition_members(self, partition: str) -> Iterator[int]:
        for holder_id, _holder, local_index in self.iter_partition_entries(
            partition,
        ):
            yield self.locator(holder_id, local_index)

    def find_record(
        self,
        partition: str,
        record_id: str,
        haplotype: str,
    ) -> int:
        try:
            holder_id = self.haplotype_ids[haplotype]
            partition_id = self.partition_ids[partition]
        except KeyError as error:
            raise KeyError(record_id) from error
        prefix = f"g{partition}_{haplotype}_"
        if not record_id.startswith(prefix):
            raise KeyError(record_id)
        try:
            sequence_index = int(record_id[len(prefix):])
        except ValueError as error:
            raise KeyError(record_id) from error
        holder = self.holders[holder_id]
        if holder is None or sequence_index <= 0:
            raise KeyError(record_id)
        start, end = holder.partition_range(partition_id)
        local_index = start + sequence_index - 1
        if (
            local_index >= end
            or holder.sequence_index[local_index] != sequence_index
        ):
            raise KeyError(record_id)
        return self.locator(holder_id, local_index)

    def partition_name(self, locator: int) -> str:
        holder, local_index = self.decode_locator(locator)
        return self.partitions[holder.partition_id[local_index]]

    def haplotype(self, locator: int) -> str:
        holder, _local_index = self.decode_locator(locator)
        return holder.assembly.name

    def contig(self, locator: int) -> str:
        holder, local_index = self.decode_locator(locator)
        return holder.contigs[holder.contig_id[local_index]]

    def start(self, locator: int) -> int:
        holder, local_index = self.decode_locator(locator)
        return int(holder.source_start[local_index])

    def end(self, locator: int) -> int:
        holder, local_index = self.decode_locator(locator)
        return int(holder.source_end[local_index])

    def record_id(self, locator: int) -> str:
        holder, local_index = self.decode_locator(locator)
        return (
            f"g{self.partitions[holder.partition_id[local_index]]}_"
            f"{holder.assembly.name}_"
            f"{holder.sequence_index[local_index]}"
        )

    def record_length(self, locator: int) -> int:
        holder, local_index = self.decode_locator(locator)
        return int(
            holder.source_end[local_index] - holder.source_start[local_index]
        )


def read_query_sources(
    path: str,
    hotspot_dir: str,
) -> List[AssemblySource]:
    sources: List[AssemblySource] = []
    seen = set()
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) != 2:
                raise ValueError(
                    f"{path}:{line_number}: expected exactly NAME FASTA; index must be FASTA.fai"
                )
            name, fasta = fields[0], os.path.abspath(os.path.expanduser(fields[1]))
            fai = fasta + ".fai"
            hotspot = os.path.abspath(
                os.path.join(hotspot_dir, f"{name}_hotspot.txt")
            )
            if name in seen:
                raise ValueError(f"{path}:{line_number}: duplicate {name!r}")
            for required in (fai, hotspot):
                if not os.path.isfile(required):
                    raise FileNotFoundError(required)
            seen.add(name)
            sources.append(AssemblySource(
                len(sources), name, fasta, fai, hotspot,
            ))
    if not sources:
        raise ValueError(f"{path}: no assemblies")
    return sources


def read_blocks(path: str) -> Dict[int, BlockSource]:
    blocks: Dict[int, BlockSource] = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                fields = raw.split()
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected BED4+")
            start, end, index = int(fields[1]), int(fields[2]), int(fields[3])
            if start < 0 or end <= start or index <= 0 or index in blocks:
                raise ValueError(
                    f"{path}:{line_number}: invalid/duplicate block {index}"
                )
            blocks[index] = BlockSource(
                index, str(index), fields[0], start, end,
            )
    if not blocks:
        raise ValueError(f"{path}: no blocks")
    return blocks


def partition_name_from_target_path(path: str) -> str:
    """Recover the partition name from one target-list FASTA path."""
    normalized = path.rstrip("/")
    filename = os.path.basename(normalized)
    lowered = filename.lower()
    name = filename
    for suffix in (
        ".fasta_kmer.list", ".fna_kmer.list", ".fa_kmer.list",
        "_kmer.list", ".fasta", ".fna", ".fa", ".list",
    ):
        if lowered.endswith(suffix):
            name = filename[:-len(suffix)]
            break
    parent = os.path.basename(os.path.dirname(normalized))
    return parent if parent and parent == name else name


def remap_blocks_to_target_list(
    blocks: Dict[int, BlockSource],
    target_list: str,
) -> Dict[int, BlockSource]:
    """Map KmerSearcher target-list row IDs to original block coordinates."""
    remapped: Dict[int, BlockSource] = {}
    seen_partitions = set()
    with open(target_list, "rt") as handle:
        for target_index, raw in enumerate(handle, 1):
            target_path = raw.strip()
            if not target_path:
                raise ValueError(
                    f"{target_list}:{target_index}: empty target-list row"
                )
            partition = partition_name_from_target_path(target_path)
            try:
                original_index = int(partition)
            except ValueError as error:
                raise ValueError(
                    f"{target_list}:{target_index}: partition {partition!r} "
                    "is not a positive integer block name"
                ) from error
            if original_index <= 0 or original_index not in blocks:
                raise ValueError(
                    f"{target_list}:{target_index}: partition {partition!r} "
                    "is absent from the original block BED"
                )
            if partition in seen_partitions:
                raise ValueError(
                    f"{target_list}:{target_index}: duplicate partition "
                    f"{partition!r}"
                )
            original = blocks[original_index]
            remapped[target_index] = BlockSource(
                target_index,
                partition,
                original.contig,
                original.start,
                original.end,
            )
            seen_partitions.add(partition)
    if not remapped:
        raise ValueError(f"{target_list}: no target FASTA paths")
    return remapped


def read_fai_lengths(path: str) -> Dict[str, int]:
    lengths: Dict[str, int] = {}
    with open(path, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip():
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected FAI2+")
            lengths[fields[0]] = int(fields[1])
    if not lengths:
        raise ValueError(f"{path}: empty FAI")
    return lengths


def merge_and_anchor(
    intervals: Iterable[Tuple[int, int]],
    contig_length: int,
    merge_distance: int,
    anchor: int,
) -> List[Tuple[int, int]]:
    merged: List[List[int]] = []
    for start, end in sorted(intervals):
        if not merged or start - merged[-1][1] >= merge_distance:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [
        (max(0, start - anchor), min(contig_length, end + anchor))
        for start, end in merged
    ]


def extend_reference_overlaps(
    intervals: Iterable[Tuple[int, int]],
    assembly: str,
    contig: str,
    contig_length: int,
    block: Optional[BlockSource],
    extension: int,
) -> List[Tuple[int, int]]:
    translated: Optional[Tuple[int, int]] = None
    if block is not None:
        if block.contig == contig:
            translated = block.start, block.end
        else:
            prefix = f"{assembly}_{contig}_"
            if block.contig.startswith(prefix):
                coordinate_fields = block.contig[len(prefix):].split("_")
                if len(coordinate_fields) == 2:
                    try:
                        source_start, source_end = map(
                            int, coordinate_fields,
                        )
                    except ValueError:
                        pass
                    else:
                        source_length = source_end - source_start
                        if (
                            source_start >= 0
                            and source_length > 0
                            and block.start >= 0
                            and block.end > block.start
                            and block.end <= source_length
                        ):
                            translated = (
                                source_start + block.start,
                                source_start + block.end,
                            )
    output: List[Tuple[int, int]] = []
    for start, end in intervals:
        if (
            translated is not None
            and start < translated[1]
            and translated[0] < end
        ):
            start = max(0, start - extension)
            end = min(contig_length, end + extension)
        if output and start <= output[-1][1]:
            output[-1] = output[-1][0], max(output[-1][1], end)
        else:
            output.append((start, end))
    return output


def reconstruct_one_assembly(
    source: AssemblySource,
    blocks: Dict[int, BlockSource],
    valid_partitions: frozenset,
    partition_ids: Dict[str, int],
    kmer_length: int,
    merge_distance: int,
    anchor: int,
    reference_extension: int,
) -> AssemblyHeaderStore:
    lengths = read_fai_lengths(source.fai)
    raw_regions = defaultdict(lambda: defaultdict(list))
    hotspot_count = 0
    with open(source.hotspot, "rt") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) < 4:
                raise ValueError(
                    f"{source.hotspot}:{line_number}: expected "
                    "CONTIG BLOCK START END"
                )
            contig = fields[0]
            block_index = int(fields[1])
            hotspot_start, hotspot_end = int(fields[2]), int(fields[3])
            if (
                block_index <= 0
                or hotspot_start <= 0
                or hotspot_end < hotspot_start
            ):
                raise ValueError(
                    f"{source.hotspot}:{line_number}: invalid hotspot"
                )
            block = blocks.get(block_index)
            if block is None or block.name not in valid_partitions:
                continue
            start = max(0, hotspot_start - kmer_length)
            end = hotspot_end
            if start < end:
                raw_regions[block_index][contig].append((start, end))
                hotspot_count += 1

    contigs: List[str] = []
    contig_ids: Dict[str, int] = {}
    partition_column = array("I")
    contig_column = array("I")
    source_start_column = array("Q")
    source_end_column = array("Q")
    sequence_column = array("I")
    by_contig: Dict[int, array] = {}

    def intern_contig(contig: str) -> int:
        contig_id = contig_ids.get(contig)
        if contig_id is None:
            contig_id = len(contigs)
            contig_ids[contig] = contig_id
            contigs.append(contig)
            by_contig[contig_id] = array("I")
        return contig_id

    for block_index in sorted(raw_regions):
        block = blocks[block_index]
        partition_id = partition_ids[block.name]
        sequence_index = 0
        for contig in sorted(raw_regions[block_index]):
            try:
                contig_length = lengths[contig]
            except KeyError as error:
                raise KeyError(
                    f"{contig!r} from {source.hotspot} is absent from "
                    f"{source.fai}"
                ) from error
            intervals = merge_and_anchor(
                raw_regions[block_index][contig],
                contig_length,
                merge_distance,
                anchor,
            )
            intervals = extend_reference_overlaps(
                intervals,
                source.name,
                contig,
                contig_length,
                block,
                reference_extension,
            )
            contig_id = intern_contig(contig)
            for start, end in intervals:
                sequence_index += 1
                local_index = len(partition_column)
                partition_column.append(partition_id)
                contig_column.append(contig_id)
                source_start_column.append(start)
                source_end_column.append(end)
                sequence_column.append(sequence_index)
                by_contig[contig_id].append(local_index)
    return AssemblyHeaderStore(
        source,
        contigs,
        partition_column,
        contig_column,
        source_start_column,
        source_end_column,
        sequence_column,
        by_contig,
        hotspot_count,
    )


_WORKER_BLOCKS: Optional[Dict[int, BlockSource]] = None
_WORKER_VALID_PARTITIONS = frozenset()
_WORKER_PARTITION_IDS: Optional[Dict[str, int]] = None
_WORKER_SETTINGS: Optional[Tuple[int, int, int, int]] = None


def _initialize_reconstruction_worker(
    blocks: Dict[int, BlockSource],
    valid_partitions: frozenset,
    partition_ids: Dict[str, int],
    settings: Tuple[int, int, int, int],
) -> None:
    global _WORKER_BLOCKS
    global _WORKER_VALID_PARTITIONS
    global _WORKER_PARTITION_IDS
    global _WORKER_SETTINGS
    _WORKER_BLOCKS = blocks
    _WORKER_VALID_PARTITIONS = valid_partitions
    _WORKER_PARTITION_IDS = partition_ids
    _WORKER_SETTINGS = settings


def _reconstruct_worker(source: AssemblySource) -> AssemblyHeaderStore:
    if (
        _WORKER_BLOCKS is None
        or _WORKER_PARTITION_IDS is None
        or _WORKER_SETTINGS is None
    ):
        raise RuntimeError("hotspot reconstruction worker is not initialized")
    return reconstruct_one_assembly(
        source,
        _WORKER_BLOCKS,
        _WORKER_VALID_PARTITIONS,
        _WORKER_PARTITION_IDS,
        *_WORKER_SETTINGS,
    )


def reconstruct_headers(
    sources: Sequence[AssemblySource],
    blocks: Dict[int, BlockSource],
    partition_names: Sequence[str],
    jobs: int,
    kmer_length: int = 31,
    merge_distance: int = 30_000,
    anchor: int = 15_000,
    reference_extension: int = 5_000,
    progress=None,
) -> Tuple[CompactHeaderStore, int]:
    valid_partitions = frozenset(partition_names)
    ordered_partitions = [
        blocks[block_index].name
        for block_index in sorted(blocks)
        if blocks[block_index].name in valid_partitions
    ]
    if len(ordered_partitions) != len(valid_partitions):
        missing = valid_partitions - set(ordered_partitions)
        raise ValueError(
            "local-block partitions absent from block definitions: "
            + ",".join(sorted(missing))
        )
    store = CompactHeaderStore(ordered_partitions, sources)
    settings = (
        kmer_length, merge_distance, anchor, reference_extension,
    )
    total_hotspots = 0
    with ProcessPoolExecutor(
        max_workers=jobs,
        mp_context=mp_context(),
        initializer=_initialize_reconstruction_worker,
        initargs=(
            blocks,
            valid_partitions,
            store.partition_ids,
            settings,
        ),
    ) as executor:
        source_iterator = iter(sources)
        pending = set()
        for _index in range(jobs * 2):
            try:
                source = next(source_iterator)
            except StopIteration:
                break
            pending.add(executor.submit(_reconstruct_worker, source))
        finished = 0
        while pending:
            completed, pending = wait(
                pending, return_when=FIRST_COMPLETED,
            )
            for future in completed:
                holder = future.result()
                store.add_holder(holder)
                total_hotspots += holder.hotspot_count
                finished += 1
                if progress is not None:
                    progress(
                        finished, len(sources), len(store), total_hotspots,
                    )
                try:
                    source = next(source_iterator)
                except StopIteration:
                    continue
                pending.add(executor.submit(_reconstruct_worker, source))
    # Detect a missing process result before downstream code starts.
    list(store.iter_holders())
    return store, total_hotspots


def source_fingerprints(
    query_paths: str,
    blocks_bed: str,
    sources: Sequence[AssemblySource],
    additional_paths: Sequence[str] = (),
) -> Iterator[Tuple[str, int, int]]:
    for path in (query_paths, blocks_bed, *additional_paths):
        file_stat = os.stat(path)
        yield os.path.abspath(path), file_stat.st_size, file_stat.st_mtime_ns
    for source in sources:
        for path in (source.fai, source.hotspot):
            file_stat = os.stat(path)
            yield path, file_stat.st_size, file_stat.st_mtime_ns
