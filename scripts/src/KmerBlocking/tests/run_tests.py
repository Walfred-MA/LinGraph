#!/usr/bin/env python3

import pathlib
import struct
import subprocess
import sys
import tempfile


UINT32_MAX = (1 << 32) - 1
SEARCHER_HASH1_SIZE = 0x40000000
SEARCHER_VALUE_BUCKETS = 1 << 16


def write_fai(path, contigs):
    with path.open("w", encoding="ascii") as output:
        for name, length in contigs:
            output.write(f"{name}\t{length}\t0\t80\t81\n")


def write_index(path, contigs, groups):
    offsets = {}
    offset = 0
    for name, length in contigs:
        offsets[name] = offset
        offset += length
    with path.open("wb") as output:
        for key in sorted(groups):
            locations = sorted(offsets[name] + position
                               for name, position in groups[key])
            output.write(struct.pack("<Q", key))
            if len(locations) == 1:
                output.write(struct.pack("<I", locations[0]))
            else:
                output.write(struct.pack("<I", UINT32_MAX - len(locations)))
                output.write(struct.pack(f"<{len(locations)}I", *locations))


def run(binary, *arguments, ok=True):
    result = subprocess.run(
        [binary, *map(str, arguments)], text=True, capture_output=True
    )
    if ok and result.returncode != 0:
        raise AssertionError(result.stderr)
    if not ok and result.returncode == 0:
        raise AssertionError("command unexpectedly succeeded")
    return result


def read_bed(path):
    records = []
    for line in path.read_text(encoding="ascii").splitlines():
        name, start, end = line.split("\t")
        records.append((name, int(start), int(end)))
    return records


def hotspot_groups(query_contig, target_contig, target_positions,
                   first_key=1):
    groups = {}
    for index, target in enumerate(target_positions):
        key = first_key + index
        groups[key] = [(query_contig, 1_000 + index),
                       (target_contig, target)]
    return groups


def rotate_right_32(value, amount):
    value &= UINT32_MAX
    return ((value >> amount) |
            (value << ((32 - amount) & 31))) & UINT32_MAX


def searcher_hash(kmer):
    mask30 = (1 << 30) - 1
    hash2 = (kmer >> 30) & UINT32_MAX
    mixed = (rotate_right_32(hash2, 7) ^
             rotate_right_32(hash2, 13) ^
             rotate_right_32(hash2, 19)) & mask30
    return (kmer & mask30) ^ mixed, hash2


def verify_searcher_cache(path, expected):
    with path.open("rb") as cache:
        size, hash1_size, total_kmers, total_targets, collisions = (
            struct.unpack("<QQQQI", cache.read(36))
        )
        if (size, hash1_size, total_kmers, total_targets, collisions) != (
            len(expected) + 11,
            SEARCHER_HASH1_SIZE,
            len(expected) + 1,
            max(expected.values()) + 1,
            0,
        ):
            raise AssertionError("KmerSearcher cache header is incorrect")

        hash1_offset = 36
        hash2_offset = hash1_offset + 4 * hash1_size
        hash2_size_offset = hash2_offset + 4 * size
        values_offset = hash2_size_offset + hash1_size
        multi_sizes_offset = values_offset + 2 * size
        multi_values_offset = multi_sizes_offset + 4 * SEARCHER_VALUE_BUCKETS
        expected_size = multi_values_offset + 8 * SEARCHER_VALUE_BUCKETS
        if path.stat().st_size != expected_size:
            raise AssertionError("KmerSearcher cache file size is incorrect")

        for kmer, block_index in expected.items():
            hash1, hash2 = searcher_hash(kmer)
            cache.seek(hash1_offset + 4 * hash1)
            redirect = struct.unpack("<I", cache.read(4))[0]
            cache.seek(hash2_size_offset + hash1)
            bucket_size = cache.read(1)[0]
            found = False
            for index in range(redirect, redirect + bucket_size):
                cache.seek(hash2_offset + 4 * index)
                candidate = struct.unpack("<I", cache.read(4))[0]
                if candidate == hash2:
                    cache.seek(values_offset + 2 * index)
                    value = struct.unpack("<H", cache.read(2))[0]
                    if value != block_index:
                        raise AssertionError(
                            f"k-mer {kmer} has block {value}, "
                            f"expected {block_index}"
                        )
                    found = True
                    break
            if not found:
                raise AssertionError(f"k-mer {kmer} is absent from cache")

        cache.seek(multi_sizes_offset)
        if cache.read(4 * SEARCHER_VALUE_BUCKETS) != (
            b"\0" * (4 * SEARCHER_VALUE_BUCKETS)
        ):
            raise AssertionError("single-value cache has multi-value counts")
        cache.seek(multi_values_offset)
        if cache.read(8 * SEARCHER_VALUE_BUCKETS) != (
            b"\0" * (8 * SEARCHER_VALUE_BUCKETS)
        ):
            raise AssertionError("single-value cache has multi-value records")


def main():
    binary = str(pathlib.Path(sys.argv[1]).resolve())
    with tempfile.TemporaryDirectory() as temporary:
        directory = pathlib.Path(temporary)
        contigs = [("shift", 120_000), ("quiet", 75_000)]
        fai = directory / "assembly.fa.fai"
        index = directory / "assembly.idx"
        write_fai(fai, contigs)

        # The first 100 target hits fit within 500 bp. Two more hits remain
        # connected by gaps below 150 bp, so the entire component ends at
        # 52,749 + 31 = 52,780.
        targets = [52_000 + 5 * index for index in range(100)]
        targets.extend([52_600, 52_749])
        groups = hotspot_groups("shift", "shift", targets)
        write_index(index, contigs, groups)

        one = directory / "one.bed"
        four = directory / "four.bed"
        common = ("-i", index, "-f", fai, "--block-size", 50_000,
                  "--search-distance", 10_000)
        run(binary, *common, "-o", one, "-t", 1)
        run(binary, *common, "-o", four, "-t", 4)
        if one.read_bytes() != four.read_bytes():
            raise AssertionError("thread counts produced different BED output")
        expected = [
            ("shift", 0, 52_780),
            ("shift", 52_780, 100_000),
            ("shift", 100_000, 120_000),
            ("quiet", 0, 50_000),
            ("quiet", 50_000, 75_000),
        ]
        if read_bed(one) != expected:
            raise AssertionError(f"wrong hotspot adjustment: {read_bed(one)}")

        # Exactly 99 dense hits are not sufficient.
        no_hotspot_index = directory / "no_hotspot.idx"
        no_hotspot = hotspot_groups(
            "shift", "shift", [52_000 + 5 * i for i in range(99)])
        write_index(no_hotspot_index, contigs, no_hotspot)
        no_hotspot_bed = directory / "no_hotspot.bed"
        run(binary, "--index", no_hotspot_index, "--fai", fai,
            "--output", no_hotspot_bed, "--threads", 2)
        if read_bed(no_hotspot_bed)[0:3] != [
            ("shift", 0, 50_000),
            ("shift", 50_000, 100_000),
            ("shift", 100_000, 120_000),
        ]:
            raise AssertionError("99 hits incorrectly formed a hotspot")

        # A hotspot that would leave under 10 kb absorbs the entire neighbor.
        absorb_index = directory / "absorb.idx"
        absorb = hotspot_groups(
            "shift", "shift", [91_000 + i for i in range(100)])
        write_index(absorb_index, contigs, absorb)
        absorb_bed = directory / "absorb.bed"
        run(binary, "--index", absorb_index, "--fai", fai,
            "--output", absorb_bed, "--search-distance", 50_000)
        if read_bed(absorb_bed)[0:2] != [
            ("shift", 0, 100_000),
            ("shift", 100_000, 120_000),
        ]:
            raise AssertionError("small neighboring remainder was not absorbed")

        bad = run(binary, "--index", index, "--fai", fai,
                  "--output", directory / "bad.bed", "--block-size", 9_999,
                  ok=False)
        if "10000" not in bad.stderr:
            raise AssertionError("block-size validation message is missing")

        # Independent block annotation/hash mode. BED input is intentionally
        # unsorted; count>1 k-mers must not be loaded into this mode.
        annotation_bed = directory / "annotation-input.bed"
        annotation_bed.write_text(
            "quiet\t0\t1000\n"
            "shift\t500\t800\n"
            "shift\t100\t500\n",
            encoding="ascii",
        )
        annotation_index = directory / "annotation.idx"
        annotation_groups = {
            100: [("shift", 50)],       # outside every block
            101: [("shift", 100)],      # block 1 start
            102: [("shift", 499)],      # block 1 end - 1
            103: [("shift", 500)],      # block 2 start
            104: [("quiet", 10)],       # block 3
            200: [("shift", 200), ("quiet", 20)],  # count 2: ignored
            0x3FFFFFFF: [("shift", 101)],  # maximum 30-bit hash bucket
        }
        write_index(annotation_index, contigs, annotation_groups)
        annotated = directory / "annotated.bed"
        block_hash = directory / "blocks.hash.bin"
        run(binary, "-b", annotation_bed, "-k", annotation_index,
            "-f", fai, "-ob", annotated, "-ok", block_hash, "-t", 4)
        if annotated.read_text(encoding="ascii").splitlines() != [
            "shift\t100\t500\t1\t3",
            "shift\t500\t800\t2\t1",
            "quiet\t0\t1000\t3\t1",
        ]:
            raise AssertionError("annotated BED is incorrect")

        verify_searcher_cache(
            block_hash,
            {101: 1, 102: 1, 103: 2, 104: 3, 0x3FFFFFFF: 1},
        )

    print("all tests passed")


if __name__ == "__main__":
    main()
