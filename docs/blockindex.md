# Build a graph for a new species from one reference

Use `blockindex` to generate grouped blocks, then `graph --construct-only`
to build the reusable local graph package. Neither step requires a cohort.

Start with an uncompressed reference FASTA and its adjacent `samtools faidx`
index. Use the same soft masking as in your established blocking procedure:
lowercase sequence is excluded from the merger's unmasked k-mers, and
`--unmask` measures uppercase A/C/G/T bases. The input FASTA and its contig
names are preserved. `Species_h1` is a sample label, not a restriction on
species or ploidy.

```bash
# Prerequisite, if reference.fa.fai does not already exist:
samtools faidx /assemblies/reference.fa

# 1. Generate reference blocks and the partition search index.
python scripts/LinGraph.py blockindex \
  -i /assemblies/reference.fa -s Species_h1 \
  -o species_blocks -t 64

# 2. Build/cache the single-reference graph, keeping its templates fixed.
python scripts/LinGraph.py graph \
  -I species_blocks/query_paths.txt -r Species_h1 \
  -b species_blocks/blocks.bed --bed-grouped \
  -G species_graph -t 64 \
  --construct-only --static-block
```

`--construct-only` stops after the existing graph builder has generated its
graphs, reference alignments, breakpoint caches, and compact summary. It does
not run cohort variant calling, merging, or merged-VCF GFA export. It accepts
one or more genomes. For a single genome it also accepts native reference
contig names, such as `chr1`, without adding a sample prefix.

`--static-block` is independent: it skips additional sequence discovery and
local template refinement. Omit it when normal graph refinement is wanted.
`--construct-only` does not implicitly enable it. Existing graph runs without
`--construct-only` retain their build/call/merge behavior and still require
`-O`. Construction-only runs can omit `-O`; their CLI logs and metadata then
live under `-G/lingraph/`. If supplied, `-O` can equal `-G` or be a separate,
non-nested directory. Calling-only options and `--make-graph` are incompatible
with construction-only mode.

## Block generation

The wrapper runs these existing algorithms in dependency order:

1. First-genome self-cleaning from `minsetref.py`, exposed by
   `self_clean_reference.py`; retain source coordinates in FASTA headers.
2. Index the cleaned FASTA; run `kmerindex` and `kmerblocking`.
3. Assign `g1`, `g2`, etc.; extract and index the named block FASTA.
4. Map the **original reference** back to those templates using
   `map_assemblies_to_block_reference.py`.
5. Filter mappings with `--minimum-block-size 20000`.
6. Merge partitions with `kmerpartitionmerger -m 10 --only-kmer`.
7. Materialize local blocks using the filtered mappings, merged BED, and
   merger groups TSV; retain groups with at least 1,000 uppercase bases.
8. Export the retained source intervals and index their partition FASTAs
   using the existing KmerSearcher dummy-query method (`-c 500`).

The default initial block size is 50,000 bp. Adjusted boundaries and merged
groups need not be 50 kb. `--block-size`, `--minimum-block-size`, and
`--unmask` are configurable. `--unmask` is a retention threshold, not a
sequence-unmasking operation. The original reference's mapped repeat copies
and retained template origins contribute to the final BED, including template
origins without a surviving mapping.

## Outputs and resume

| File or directory | Content |
| --- | --- |
| `blocks.bed` | BED4: original contig, zero-based start, exclusive end, final group |
| `query_paths.txt` | One `Species_h1 FASTA` row, with an absolute reference path |
| `localblocks/GROUP/` | Retained template FASTA and detailed source-coordinate BED |
| `partition_targets.list` | Ordered absolute paths to the retained template FASTAs |
| `partition_targets.list.bin` | Matching KmerSearcher cache |
| `partition_targets.list.bin.fasta_target_unique_v2` | KmerSearcher cache-version marker |
| `uncovered.bed` | Reference spans absent from the final retained block union |
| `blocks.stats.json` | Partition/interval counts and reference coverage totals |
| `blockindex.json` | Reference file identity and blocking parameters |
| `work/`, `checkpoints/` | Intermediate artifacts and completed-stage inventories |

Repeated BED group names indicate related loci and must be passed to graph
construction with `--bed-grouped`, even when every group has only one row.
The graph builder generates its own final caches from this BED; the block
index targets are the initial templates, not the final graph paths.

Repeat the same command to resume; `--resume` is accepted but not required.
Completed stages are reused only when their recorded inputs and outputs match.
Failed stages have no completion marker. Missing artifacts are rebuilt before
their dependents. A changed reference or blocking policy requires a new output
directory, so it cannot silently reuse a different genome's cache. Use a
dedicated directory and do not move it while using the absolute-path target
list. A lock prevents simultaneous writers to the same run.

`--dry-run` validates input paths/index geometry and prints the stages without
creating files or running tools. It does not test dependency availability or
execute alignments. The self-cleaning wrapper uses this checkout's existing
first-genome algorithm; historical results from other versions should be
compared before treating them as identical.

Dependencies are the existing native k-mer programs under `scripts/`, plus
`samtools`, `minimap2`, `winnowmap`, `meryl`, `blastn`, `makeblastdb`, and `sort`.
Graph construction additionally requires the existing pinned Snakemake
environment and graph-alignment dependencies described in
[`graph_build_snakemake/README.md`](../scripts/graph_build_snakemake/README.md).
The existing KmerIndex 32-bit cumulative-coordinate limit still applies to
the cleaned reference; this wrapper does not extend the native index format.
