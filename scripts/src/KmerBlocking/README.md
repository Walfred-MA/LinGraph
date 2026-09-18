# KmerBlocking

`kmerblocking` converts a KmerIndex binary index into adjusted, non-overlapping
assembly blocks in BED3 format. It does not modify the FASTA.

For every contig, the program creates a `uint64` array with one entry per base,
initialized to zero, and loads indexed canonical 31-mer keys at their start
locations. Contigs are processed independently and in parallel.

## Build and run

```sh
make
./kmerblocking \
  -i assembly.k31.idx \
  -f assembly.fa.fai \
  -o assembly.blocks.bed \
  -t 8
```

The initial block size is 50,000 bases (`--block-size` or `-S`). For each
current block, its indexed keys are queried against the external upstream and
downstream flanks only. The flank length is 10,000 bases by default
(`--search-distance` or `-D`). Each indexed flank occurrence is one shared hit.

Hotspot detection walks shared-hit starts in coordinate order:

1. Consecutive hits belong to the same connected component when their distance
   is strictly less than 150 bases.
2. A component is a hotspot when any 100 consecutive hits span at most 500
   bases.
3. Once qualified, the entire connected component is retained. Its BED-like
   interval runs from the first hit start through the final 31-mer.

The current block absorbs every qualifying flank hotspot and the gap between
the hotspot and block. This moves the adjacent breakpoint. If that leaves less
than 10,000 bases in the neighboring block, the whole neighboring block is
absorbed. Blocks are processed from left to right using all earlier boundary
updates. An initial terminal remainder shorter than 10,000 bases is also merged
into the preceding block.

Output has no header and uses standard zero-based, half-open BED3 coordinates:

```text
contig_name<TAB>start<TAB>end
```

The `.fai` must be the same file/order used by KmerIndex. Run `make test` for
integration tests.

## Independent block annotation and unique-k-mer hash mode

This mode is independent of hotspot-based block creation. It reads an existing
BED and KmerIndex file, loads only k-mers with occurrence count 1, and assigns
each unique k-mer start to a block:

```sh
./kmerblocking \
  -b assembly.blocks.bed \
  -k assembly.k31.idx \
  -f assembly.fa.fai \
  -ob assembly.blocks.scored.bed \
  -ok assembly.blocks.hash.bin \
  -t 64
```

`-ob` and `-ok` are independently optional, but at least one is required.
Blocks are sorted in `.fai` coordinate order and receive one-based indexes.
Block-start boundary records use `UINT64_MAX` in the combined 12-byte
`(uint64 kmer, uint32 position)` event vector, so a k-mer starting exactly at a
boundary is assigned to the new block.

The annotated BED is BED5:

```text
contig  start  end  block_index  unique_kmer_count
```

The score is the raw count and is not capped at 1000.

The optional `-ok` output is a native KmerSearcher cache. It uses the same
30-bit mixed `hashfunc` and the exact layout written by
`Kmer_hash::savehash`: header, `Hash1`, `Hash2`, `Hash2_size`, `values`,
`values_mul_bucksize`, and the 65,536 multi-value buckets. Because each
count-1 k-mer belongs to exactly one block, its one-based block index is stored
directly in the `uint16` `values` array and all multi-value buckets are empty.
Both tools use `0x40000000` first-level entries so every value from the 30-bit
hash, including `0x3fffffff`, is addressable. Caches made with the older
`0x3fffffff` table-size definition must be regenerated.

KmerSearcher reserves `UINT16_MAX` as its multi-value marker, so this output
supports at most 65,534 BED blocks. The file has an approximately 5 GiB fixed
logical size plus six bytes per assigned k-mer. KmerBlocking writes empty
ranges sparsely, although loading it in KmerSearcher still requires the full
KmerSearcher table memory.

For `kmer_searcher -T targets.txt`, KmerSearcher looks for
`targets.txt.bin`. Either write `-ok targets.txt.bin` or rename the completed
cache to that exact path. `targets.txt` itself remains a normal text target
list; do not pass the binary cache as the `-T` argument.
