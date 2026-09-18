# KmerPartitionMerger

Merge partitions using the reference mappings in the mapper's filtered
16-column BED and shared canonical unmasked 31-mers. Only windows consisting
entirely of uppercase `A/C/G/T` contribute; lowercase soft-masked bases and
invalid bases exclude every 31-mer window containing them.

For each column-4 partition, the k-mer input is the coordinate union of:

1. every lifted interval in BED columns 1-3; and
2. the partition's own source interval decoded from column 4.

Column-4 names must use
`gN_SAMPLE_HAP_CONTIG_START_END`. The parser removes the first three
underscore-separated fields (`gN`, `SAMPLE`, `HAP`), treats the final two as
zero-based half-open coordinates, and rejoins every intervening field with `_`
to recover the true contig. For example,
`g60845_HG38_h1_chr6_GL000256v2_alt_3855379_3906047` decodes to
`chr6_GL000256v2_alt:3855379-3906047`. The target FASTA record with that name
must be exactly 50,668 bp, matching the decoded zero-based half-open span.

```bash
make -C KmerPartitionMerger -j

./KmerPartitionMerger/kmerpartitionmerger \
  -f cohort_minsetref_v3/full_ref_withnovels.fa \
  -T cohort_minsetref_v3/main_chroms_withnovel_blocks.fasta \
  -b cohort_minsetref_v3/main_chroms_withnovel_blocks.fasta_aligned.filtered.bed \
  -o cohort_minsetref_v3/main_chroms_withnovel_blocks.merged.bed \
  -m 20 \
  -t 64
```

`-f` supplies sequences for lifted BED columns 1-3. `-T` supplies the exact
sequence of every column-4 source interval and enumerates all targets, including
targets with no filtered BED alignment at all. No `query_paths.txt`, original
assembly FASTA, or `.fai` is needed.

The first rule merges partitions whose unioned query-coordinate mappings share
strictly more than 1,000 uppercase `A/C/G/T` bases. The resulting groups are
then used for k-mer association. The first rule intentionally uses only lifted
BED intervals; adding the column-4 source interval does not manufacture a
reference-liftover overlap. K-mer extraction uses the union described above.
A distinct canonical unmasked k-mer contributes once per group, is discarded when
associated with more than `-m` groups, and pair counters saturate at
`UINT16_MAX`. The second rule merges group pairs sharing strictly more than
1,000 retained k-mers. Both rules are transitively closed by a minimum-root
disjoint set.

To apply only the aligned-interval rule and completely skip unmasked-k-mer
extraction, sorting, pair counting, and the second merge rule, add:

```bash
--alignment-only
```

The `-p` unmasked-k-mer pair report is unavailable in this mode. The output BED4
and groups TSV formats are unchanged.

To skip aligned-reference-base merging and merge the original BED partitions
only through shared unmasked 31-mers, add:

```bash
--only-kmer
```

`-p` remains available in this mode and reports shared unmasked k-mer counts
between the original partitions. `--only-kmer` and `--alignment-only` are
mutually exclusive.

Regions on the same FASTA record are merged across gaps strictly below 300 bp.
The BED output remains BED4 in the lifted coordinate system; decoded column-4
source intervals expand only the k-mer input and are not emitted as artificial
lifted intervals. The groups TSV maps each original column-4 target to its final
name, such as `merged_g1g2g5`. It also contains fully unlifted targets even when
their final group consequently has no row in the output BED.
