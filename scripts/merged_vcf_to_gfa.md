# Cohort SNP/indel/SV VCFs to rGFA

The default `--gfa-mode rgfa` exports a reference backbone with insertion,
SNP, substitution, and deletion branches and stable coordinates. It includes
all sizes by default (`--all`). Add `--svonly` to include only non-SNP SVs
at or above `--svcutoff` (default 20 bp). This includes qualifying insertions
and deletions; `--insertion-only` remains a separate filter. `--gfa-mode query` retains the previous local
query-flank GFA behavior; that mode is not suitable as an anchored rGFA.

```bash
python merged_vcf_to_gfa.py \
  -v cohort_calls/cohort.all.vcf \
  -q query_paths.txt \
  --graph-folder graph \
  --local-reference-templates cohort_calls/checkpoints/local_reference_templates.fa \
  --local-path-fasta graph/summary/alternatives.fasta \
  -o cohort_calls/cohort.gfa \
  --anchor 50 --max-node-length 1024 --processes 8
```

`-v/--vcf` accepts multiple files and may be repeated. Plain and gzip-compressed
VCFs may be mixed. Shared contig headers are combined in input order; conflicting
lengths and duplicate variant IDs are rejected. Parent insertion definitions
can live in another input file, regardless of input order.

Use `--insertion-only [SIZE]` to retain the previous insertion-only behavior.
A bare `--insertion-only` uses 50 bp; `--insertion-only 20` uses 20 bp. This flag
is also available in `LinGraph.py ... --make-graph`. The older `--size-cutoff`
option remains available as a general minimum variant size. Use
`--insertion-only` with the legacy `--gfa-mode query` for mixed variant inputs.

SNP positions are 1-based base coordinates: a SNP at POS 5 replaces reference
interval `[4,5)`. Symbolic INS/SUB records use the caller's existing POS/END
breakpoint convention. SNP and substitution bases are read from their sample
assembly intervals and verified against ALT or INFO/SEQ. Deletions introduce
no sequence nodes: their branch skips `abs(SVLEN)` parent bases starting at
the normalized POS. A merged DEL's `END` is the cluster-wide envelope, which
can differ from the representative deletion's endpoint. Only legacy DEL rows
without `SVLEN` fall back to the POS/END interval.
A deletion's named path contains these flanks (at least one adjacent base on
each available side), and the reference path retains the deleted sequence.
SNP and substitution paths contain their alternate bases, like insertion cores.

Use the completed, correct merged VCFs and the actual reference/assembly FASTA
paths. A `.partial.vcf` or an input with unresolved query coordinates is not a
verified replacement for the completed input. The query list is `NAME FASTA
[FAI]` per line; relative paths are resolved against the list's directory.

`-r` is optional: the default is
`GRAPH/inputs/reference_alternatives_novels.fa`, with `GRAPH=graph` by default.
This includes the reference, imported alternatives, and discovered novels
without passing the novel catalog again. An explicit `-r selected.fa` replaces
that default with the chosen haplotype's backbone. Use
`--reference-haplotype SAMPLE_h1` with a plain assembly FASTA; otherwise the
sample is inferred from catalog `source=` metadata where available.

The cohort caller permits any reference haplotype in its assembly list and
defaults to the first assembly when its `-r` is omitted. `LinGraph.py -r NAME`
uses that sample's FASTA for calling and GFA export. Without a LinGraph `-r`,
export uses the combined default catalog. Changing a backbone requires calling
against it; a VCF from another backbone cannot be changed by relabeling coordinates.

For every local block with no accepted alignment for the chosen haplotype,
`local_reference_templates.py` selects its reference-class template paths.
Imported alternatives are fixed independent sequences. Their original bytes,
orientation, and bounds are preserved through graph construction, refinement,
packaging, and calling. `sequence_role=imported_alternative` records this policy;
source sample/contig fields remain provenance and do not authorize fetching or
expanding the template from a source assembly. Newly discovered novel templates
remain eligible for refinement.

`lift_local_templates.py` aligns a fixed imported template's actual sequence to
the chosen backbone with minimap2. It does not fetch source-assembly flanks for
these records. For refinable novel templates, it verifies source sequence and
adds 10 kb of source flanks for placement. Lifting changes placement metadata,
not template bytes. Old lift annotations are replaced before the templates are
loaded by the calling stages; public IDs and local coordinates remain unchanged.

Only uniquely supported placements reaching the actual template boundary
produce backbone edges. Two supported ends attach both sides; one supported
end attaches only that side. Tied, noncollinear, disagreeing, or missing placements
retain the full template as an unplaced path. Missing source assemblies retain
refinable templates unplaced; a source/sequence mismatch is an input error.
Fixed imported templates need no source assembly. In GFA these
paths have `LS:Z:unmapped` or `LS:Z:both` and `UP:Z:unplaced`.

`--local-reference-templates` consumes this freshly lifted catalog. Templates
already contained in included source paths reuse their sequence nodes; other
selected templates add their own sequence. The converter rejects lifts labeled
for another backbone. `OUTPUT.anchors.bed.template-lifts.tsv` records the lift
status, placements, and note for each template.

For **HG38_h1**, only `chr1`–`chr22`, `chrX`, and `chrY` are eligible. The same
rule applies to backbone records, fallback sources, lift targets, cohort routing,
GenomeLift, and insertion representative selection. `chrM`, alternate, decoy,
and unplaced scaffolds are excluded. Prepared names such as `HG38#1#chr1` are
recognized.

Do not pass `graph/summary/alternatives.fasta` to `-a`: it contains duplicated
local paths, and `-a` deliberately includes every nonempty record. Instead,
`--local-path-fasta` uses that catalog only to translate missing local VCF
CHROM/mapping targets onto an included root. It never exports catalog records
as additional paths. It reads the adjacent `.fai`, selected headers, and only
the sequence intervals needed for identity checks. Index it with `samtools
faidx` first if the `.fai` is absent.

Translation requires matching source sample/contig metadata, a covering source
interval, and case-insensitive sequence identity after strand normalization.
Both `source=SAMPLE:CONTIG:START-END[+-]` and the local catalog's second-token
source format are supported. `reference=` lift annotations are not treated as
base-level alignments. Ambiguous, uncovered, or sequence-mismatched local
targets stop conversion; no local record is silently included or dropped.
For example, the 24,136-base local target
`alternative993cc4211631_NA12879#2#pat-0000632_20487305_20511441` maps to
the included root `NA12879_h2_NA12879#2#pat-0000632_20484104_20527047`
at `[3201,27337)` on the plus strand. The root keeps its full original path.

All eligible nonempty records from the reference catalog and explicit `-a` become
complete backbone/source paths. Every target used by a reachable insertion
must have an insertion definition, an included FASTA record, or a verified
source translation onto an included record. Missing roots,
invalid breakpoints, dependency cycles, and ambiguous names fail before graph
output. Dependencies are emitted before their children, including necessary
ancestor insertions below the cutoff or with a non-PASS filter.

## Coordinates and paths

Each segment has `LN`, `SN`, `SO`, and `SR`. `SN` names its sequence of origin,
`SO` is a nonnegative, zero-based offset on that sequence, and `SR=0` identifies
the selected backbone's records. Source-tagged records from other haplotypes
in the combined catalog, explicit `-a` records, and fallback templates have
rank 1. Explicitly fixed imported alternatives also have rank 1 when their
provenance sample happens to match the backbone sample. New insertion
sources receive increasing ranks after all their dependencies.

The `P` path named for an insertion contains its core sequence, beginning at
offset zero; query flanks are not prepended. An aligned `=` run reuses its
target's segments, including reverse-oriented and nested mappings. Literal
`I`, `S`, and `X` runs use the insertion's own coordinates. Retained unaligned
pieces still get `#uniqueN` local paths with graph anchors.

An insertion's incoming and outgoing `L` edges attach to its parent at the
event breakpoints. Flanks extending past a parent insertion's beginning or
end continue through that parent's attachment to its own parent. The normal
parent/reference walk remains present. `--anchor 0` removes requested local
flanks but still creates attachment edges using the immediately adjacent
parent bases where available. Links have `0M` overlap and rank tags.

This reader follows the **minsetref cohort writer's breakpoint convention**:

| Parent coordinate system | Replacement interval |
| --- | --- |
| Reference or supplied alternative/novel FASTA | `[POS, END)` |
| Insertion parent, default | Start `POS - 1`; pure insertion ends at the same start; nonempty replacement ends at `END` |

For DEL records with `SVLEN`, the start follows the same POS convention, but
the end is always `start + abs(SVLEN)`. Size cutoffs also use `abs(SVLEN)`,
not the cluster envelope. For example, `POS=248386058;SVLEN=-659` deletes
`[248386058,248386717)` on a root even if the merged `END` is `248387849`.
The representative interval must still fit inside its parent; invalid intervals
are rejected rather than clipped to the FASTA length.

The caller now also rejects deletion operations outside a reference/template
before its coordinate helpers can clip them. Previously, a `108D` operation
starting at a path's length could produce `start == end == length` while
retaining a deletion size of 108; the merger then extended its END from that
size. This is separate from insertion-alignment splitting. The merged record
alone cannot establish a corrected position for such a call.

Audit all root deletions in an existing cohort before another graph export:

```bash
python3 minsetref/tools/audit_deletion_bounds.py \
  -v cohort_calls/cohort.sv.vcf \
  --fai reference.fa.fai --fai alternatives.fa.fai \
  --report cohort_calls/invalid_deletions.tsv
```

The audit writes a report and does not modify or exclude any VCF records.
It uses representative SVLEN, so an out-of-bounds cluster END alone is not
flagged when the representative fits. It checks only DEL records whose parents
appear in the supplied FASTA indexes; it does not validate nested insertion
parents, other variant types, sequences, or graph topology. Fixing the caller
does not repair already-written VCF records or require regenerating every
unaffected chromosome.

Older cohort merges could split internal `I`/`D` operations in an insertion's
template alignment into a chromosome deletion, including records named
`DEL_INS_..._DEL`. The corrected merger keeps an `INS` with `END == POS` as
an insertion even when its alignment contains deletions relative to a template.
For `INS_NC_060926.1_242696729_114963`, the original HG01109_h1 call has
`POS=END=242696729` and `SVLEN=46291`. Its reverse mapping contains 957 bases
of internal `D` operations; those are part of the insertion's sequence mapping
and must not become a chromosome deletion starting at the insertion point.
An existing split DEL row does not retain enough alignment context to repair
its parent coordinates safely. Regenerate affected cohorts from the original
per-sample VCFs with the updated `graphvcfmerge.py`, starting from a fresh scan
and merge workspace; old scan shards and refined groups already contain the
incorrect components. Changing `--nested-pos-base` or clipping to a chromosome
end does not recover that information.

The nested merger emits `POS = insertion_offset + 1`. Use
`--nested-pos-base 0` only for input that instead uses zero-based nested
breakpoints. These rules are specific to this cohort pipeline; this is not a
general converter for arbitrary VCF allele encodings. Graph mapping offsets
and query FASTA intervals are zero-based. In rGFA mode the mapping TSV's
`ref_start/ref_end` columns are also normalized to zero-based intervals.
Local source aliases are translated to included-root coordinates in both
the parent and mapped-target columns; the appended `ref_strand` column records
parent orientation. `OUTPUT.anchors.bed.local-paths.tsv` audits local name to
included-root interval/strand translations. Local flanks may extend past a
trimmed catalog record into its full included root.

`LinGraph.py --make-graph` passes the selected/default backbone, lifted template
catalog, query assembly list, and `--processes` in graph and individual modes.
With `-G graph/summary`, the default catalog is located in the parent graph
directory. Template selection/lift checkpoints record backbone identity; the
cohort DAG also tracks lifted templates as inputs to the calling stages.

## Unverifiable variants

Every literal allele must match a carrier's assembly bases at its `QUERYCOORD`
(on a `-` strand, a single SNP position is the traversal boundary, so its base
is at `QUERYCOORD - 1`). A SNP called on a reverse-traversed graph path, or
projected through a reverse mapping, has a complemented ALT while QUERYCOORD
keeps the query strand, and the VCF does not record that direction: both
orientations are tried at that position, and only a carrier base equal to ALT
is accepted. By default one mismatch stops the conversion. With
`--drop-unverified`, variants that no carrier verifies, and variants nested in
or aligned onto them, are left out of the graph instead; they are listed in
`OUTPUT.anchors.bed.unresolved.tsv` (dependents with `depends_on_dropped`) and
marked not exported (path `.`) in `OUTPUT.variants.tsv`. Bases that fail the
check are never exported.

## Alternative paths from the graph catalog

With `--alternative-catalog GRAPH/summary/alternatives.fasta` (the default when
the graph folder has it and neither `--local-path-fasta` nor
`--local-reference-templates` is given), every VCF target outside `-r` is found in
the catalog by source interval, never by name. Its source interval comes from the
VCF header's `##alternativeLocus` line, or from its own catalog header.
`--reference-haplotype` is required: records from that haplotype map straight
onto the backbone.

Graph construction can cut one alternative sequence into several local loci. On
each source contig, all other catalog records whose source intervals overlap or
touch are stitched into one forward path, `alt_SAMPLE_CONTIG_START_END`, covering
their union. Overlapping bases must agree; otherwise the pieces stay separate.
Each VCF target maps onto its merged path at its source offset.

A path built from one record uses that record's catalog placement (header fields
3-4, a graph CIGAR on a backbone contig) without realignment. Any other path is
aligned once, as a whole, with `lift_local_templates.lift_templates` (minimap2).
Only ends aligned exactly at the path boundary link into the backbone; other
paths stay in the graph with `UP:Z:unplaced`. The stitched sequences are written
to `OUTPUT.alternative_paths.fa`, and `OUTPUT.alternative_paths.tsv` lists each
path's records, placement, and VCF names. Targets that no merged path covers stop
the conversion.

## Full-locus duplications

An insertion with `INFO/PACLASS=fulllocusdup` copies a reference PA. All such
copies whose aligned source intervals (from `INFO/ALTERNATIVECIGAR`, pieces
joined by `&`) overlap on one reference path share a single new path,
`dup_<path>_<start>_<end>` with `TP:Z:duplication`, copying the union of those
intervals. Each copy's own P path walks the shared path where it aligns:
mismatches and inserted bases are literal bubbles, and deleted or uncovered PA
bases are skipped (deletions). Insertion bases outside every aligned piece stay
literal pieces. Copies of the same PA inserted at different places reuse the
same shared nodes.

A piece is used only if its bases, rebuilt from the PA sequence and the CIGAR
payloads, occur in the literal `INFO/SEQ` after the previous piece; otherwise it
stays literal. The log line `Full-locus duplications: ...` counts shared paths,
placed pieces, and pieces left literal. `graphvcfmerge.py` drops these
insertions unless it runs with `--keep-full-locus-dup-insertions`.

## Variant index sidecar

rGFA output also writes `OUTPUT.variants.tsv`, one row per exported VCF record:
`vcf_index` (data lines numbered from 0 across all `-v` inputs, in order),
`variant_id`, `path` (the allele's P line), `kind`, and `parent`,
`parent_start`, `parent_end`, `parent_strand`: the breakpoints on the parent P
path. A deletion's P line includes flanks; its allele is the parent interval.

## Per-sample GAF

`gfa_sample_gaf.py` writes `SAMPLE.gaf` for each sample VCF from its
`##pseudoLinearMapping` header and the merged-VCF alleles it carries. Pass the
same merged VCFs in the same order as to the converter:

```bash
python gfa_sample_gaf.py -g cohort.gfa -v cohort.sv.vcf cohort.indel.vcf cohort.snp.vcf \
  -s samples/*/*.vcf -q query_paths.txt -o gaf -t 8
```

It streams the merged VCFs once into per-sample shards (flushed every
`--shard-records`, default 100000, closing files after each flush), then
processes samples in parallel. Each contig's intervals are walked in query
order along the reference P path, splicing every carried variant's allele
(nested variants inside their parent insertion) at its breakpoints. One GAF
record covers a contig until the walk cannot continue: an unmapped interval,
an insertion without a graph allele, a missing link, or a mid-node jump.
Columns 10/11 assume the walked sequence equals the query; no CIGAR is written.
`nv:i:` counts spliced variants. `gaf_stats.tsv` summarizes breaks and unplaced
variants, and `variant_offsets.npy` holds each merged-VCF record's byte offset
(uint64, cumulative across files). Requires numpy.

For many samples, build the graph tables once and memory-map them in every job:

```bash
python gfa_gaf_index.py -g cohort.gfa -o gaf_index
python gfa_sample_gaf.py -g cohort.gfa --index gaf_index -v ... -s ... -o gaf
```

To stream the merged VCFs once for all samples, write shards first and give each
batch `--shards` (it then reads only its samples' shards; `variant_offsets.npy`
stays in the shard folder). Output matches one run over all samples; batches that
each stream the VCFs number contigs per batch, which can reorder `unplaced.tsv` rows.

```bash
python gfa_sample_gaf.py -g cohort.gfa --index gaf_index -v ... -s ALL... -o shards --write-shards
python gfa_sample_gaf.py -g cohort.gfa --index gaf_index --shards shards -v ... -s BATCH... -o gaf_1
```

The shard pass can also run chromosome by chromosome: `--plan-shards` cuts the
merged VCFs near CHROM changes (located by bisection; the sorted merge writes one
block per CHROM, nested INS_/DEL_/SUB_/DUP_ rows as one tail block), bundles small
pieces, splits pieces above `--shard-chunk-bytes`, and counts each chunk's lines
with `-t` processes, so every chunk knows its first VCF record and byte offset.
Each `--write-shards --chunk K` then runs on its own (in any order), and
`--finish-shards` numbers contigs over the chunks in order. Shards and GAFs are the
same as one pass.

```bash
python gfa_sample_gaf.py ... -o shards --plan-shards --shard-chunk-bytes 34359738368 -t 16
python gfa_sample_gaf.py ... -o shards --write-shards --chunk K     # each K in plan.json
python gfa_sample_gaf.py ... -o shards --finish-shards
```

The index keeps only what the walks read: segment lengths, P-line steps (with
each step's uint32 base offset within its path, so workers cache nothing per path) and
names (with a sorted name-hash table for lookups), canonical link keys, and the
variant index columns with paths resolved to numbers. Sequences, other tags and
variant IDs are dropped. Walks keep each GAF record as uint32 (path, start, end)
step ranges within a P line (reversed when start > end), so a worker's memory
follows the sample's variants, not the graph's node count; the path column is
streamed. Output is the same as without `--index`; a job refuses
an index built from a GFA or sidecar with another size or mtime.

## Sequence checks and large files

The converter reads the VCF body sequentially and retains coordinate/run
metadata and sequence digests, not the sample matrix or sequence payloads.
Query extraction uses scaffold batches and temporary files. Sequence copying
uses bounded chunks. Nearby short query intervals share a 64 KiB read window,
and their sequence checks happen in memory before buffered output, avoiding a
temporary-file write/read cycle for every SNP. Query workers handle up to 32
batches before recycling, rather than starting a new interpreter for each batch.
FASTA interval reads include wrapped lines in one contiguous read.

Before writing the GFA, the converter loads the reference and alternative
FASTA records supplying graph segments into memory once. Forked writer workers
share these immutable sequences; duplication paths reuse their backing record.
Lookup-only catalog records and query assemblies are not loaded wholesale.
Allow roughly one additional byte per loaded base, plus temporary loading
buffers (about 3 GB of retained sequence for a human reference, plus alternatives).
Use `--no-preload-reference` to keep bounded FASTA read windows instead. The
preload option applies to both normal and partitioned output and preserves the
same GFA and sidecars.

Memory still scales with the number of events, graph
segments, and edges; emitting the full reference and chopping nodes increases
both graph size and memory relative to the legacy local graph.

Topology (`gfa_topology.py`) is built from integer arrays rather than Python
sets and tuples: each cut point is one int64 key (source ID in segment order,
then coordinate), segments are pairs of consecutive keys, and links are
deduplicated with their minimum rank by sorting. The interval planner
(`gfa_interval_nodes.py`) prepares shared source intervals before topology:

1. Index `(start, end, variant ID, type)` by chromosome or parent insertion
   and sort each list. SNPs use `[POS-1, POS)`; insertion anchors have equal
   start and end. This index also supplies adjacency checks for carrier links.
2. Resolve composite variant paths in parent/copy dependency order. Each
   path refers to verified query intervals or reference/alternative intervals;
   nested `=` runs reuse already resolved sources. Complete forward copies
   share one interval table. Simple literal variants need no additional table.
3. Slice these shared tables with binary searches when building paths. Cached
   parent boundary flanks avoid repeated walks through nested ancestors. Each
   path's core is resolved once for both its P line and its link walk.
   Node chopping is computed once per source/start through the furthest used
   end, rather than generating the same cuts for every copy of that interval.

This retains run boundaries, node cuts, and the existing segment numbering;
it does not create duplicate sequence strings for copied alleles. Workers
inherit the interval tables through fork, and the parent concatenates S/P
blocks in their established order. Normal and partitioned output use the same
planner. The `Intervals: ...` log reports dependency depth and preparation time.
Log lines `Topology: ...` report the
time of each stage. Progress is written to stderr (the Slurm `.err` file),
including catalog/dependency resolution, metadata writing, and carrier junctions.

## Very large cohorts: `--partition-records`

Without it, memory grows with the total number of variants (roughly 3-4 KB
per variant). `--partition-records N` bounds it by N instead, trading extra
passes for memory (`gfa_partitioned.py`):

1. A light scan of the VCFs (no sample columns) records each record's byte
   offset, which variants are nested in or copy from which, the target names
   reachable variants need, and the full-locus duplication alignments.
2. Templates, catalog paths and shared duplication paths are resolved once,
   globally, exactly as in a one-pass run.
3. Whole variant trees (a variant plus everything nested in it or built on
   it) are packed in VCF order into partitions of about N records. Each
   partition is read by byte offset (streamed and skipped for `.vcf.gz`) and
   converted in its own forked process, so its memory is released afterwards.
4. Root segments are numbered first, as in one pass; variant ranks and
   segment IDs follow the one-pass dependency order (the partitions' orders
   merged by VCF line). Each partition then writes its S/P blocks and links;
   the blocks are copied in global rank order and links merged in ranges.

The GFA, `variants.tsv`, BED, mapping, local-path and template-lift outputs
are byte-identical to a one-pass run; only the unresolved report is ordered
by partition. A partition uses about 4 KB per record (e.g. N = 10,000,000
needs about 40 GB); the parent keeps roughly 40 bytes per variant plus about
100 bytes per reference-path segment. Temporary files beside the BED need
about the GFA's size plus ~30 GB per 100 million variants. Requires
`--gfa-mode rgfa` and a GFA output. Without `--drop-unverified`, the first
partition with an unverifiable variant stops the run (its variants are in
the unresolved report); later partitions are not checked.

The sample selector collects observations with both `SIZE == SVLEN` and query
span `== SVLEN`. In rGFA mode it first selects full-length exact representative
alignments: FORMAT `EXTENDGRAPHCIGAR=>N=` (also `>0HN=`), where `N == SVLEN`,
with `TEMPLATEOFFSET=0` when that field is supplied. If such observations exist,
only those are eligible; a failed FASTA check will not cause a non-exact member
to be substituted. Without an exact representative observation, the remaining
SIZE/span-matched candidates are eligible. VCF sample/observation order breaks
ties. The QUERYCOORD strand determines whether FASTA bases are reverse-complemented.

For `INS_NC_060938.1_48411877_572715`, this selects `GW00002_h1` directly and
extracts `GW00002#1#GWHBKCQ00000019:[51043879,51044006)` on the minus strand.
The BJ observation is not an eligible alternative to this exact representative.

In rGFA mode every declared literal block is checked against its
oriented query sequence using a case-insensitive SHA-256 digest. A mismatch,
missing scaffold, or out-of-bounds interval causes the next eligible candidate
to be tried. An alignment CIGAR alone is not accepted as proof of identity.
For plain INFO/SEQ this checks the entire insertion allele; for encoded
INFO/SEQ it checks the literal blocks, while `=` blocks reuse target nodes.

Alternative candidate coordinates are spooled to a temporary disk file beside
the output BED, retaining offsets/cursors in RAM. Retries remain grouped by
sample/scaffold and preserve successful extractions from earlier rounds. If
all candidates fail, the converter stops and writes the reason and candidate
count to the unresolved TSV. It never accepts a sequence mismatch merely to
complete the conversion. Missing literal payloads and ambiguous `M` runs also
stop conversion; supply literal sequence or explicit `=`/`X` mappings.

Sequences are written uppercase. Bases outside ACGTN are rejected in rGFA mode
because VG would convert them to N. Root extraction lengths are checked.
Graph output is written to a temporary file and replaces the requested output
only after successful completion. BED/mapping sidecars are written after query
selection and extraction succeed, using the verified candidate's coordinates;
they may still precede a later graph-writing error. Exhausted query candidates
leave an existing graph and BED/mapping files unchanged.

`--validate-only` extracts/checks sequence and resolves topology without writing
the graph. `--bed-only` writes coordinate metadata and does not perform the
sequence or VG import checks.

## VG import

Output uses numeric segment IDs, explicit complete `P` paths, forward
definitions, `0M` links, and a default maximum node length of 1024 bases.
`--max-node-length` changes that limit. Node chopping is useful for GCSA2
indexing; successful import alone does not guarantee every indexing workflow
will complete on a large or complex graph.

Import the explicit P paths without synthesizing duplicate paths from SN tags:

```bash
vg convert -g -r -1 cohort_calls/cohort.gfa > cohort_calls/cohort.vg
vg validate cohort_calls/cohort.vg
vg paths -x cohort_calls/cohort.vg -L
```

Here `-r -1` disables rGFA-derived path synthesis only; it preserves the explicit
reference and insertion P paths. VG's default rank-0 synthesis can otherwise
report the explicit reference P line as a duplicate. VG conversion does not
promise to retain every original optional GFA tag; keep the rGFA as the stable
coordinate source.

Run the regression suite with native tools installed:

```bash
VG=/path/to/vg GFATOOLS=/path/to/gfatools \
  python -m unittest discover -s tests -p 'test_gfa*.py' -v
```

With VG available, successful rGFA topology fixtures are imported, checked
with `vg validate`, and compared against the exact expected path sequences.
Coverage includes nested insertions, reverse walks, aliases, and chromosome
endpoints. The dedicated native-tool tests are explicitly skipped if their
tools are absent.

On 2026-09-13, all 40 converter tests passed without skips using native VG 1.77.0 and
gfatools 0.5-r296 on this Mac. VG is installed in the `vg-gfa-test` conda
environment. Another 23 LinGraph, 13 cohort, and 9 template/lift tests passed,
including a real minimap2 alignment. Snakemake 6.15.1 dry-runs verified the
template selection/lift DAG for default and switched backbones. Tests cover every root/local strand combination, reused
mapped runs, flanks beyond trimmed local paths, empty novel catalogs, and
failure on ambiguous or mismatched source mappings. This validates the fixtures, not the full mounted cohort or
all downstream indexing workflows.

The exporter represents retained **insertion alleles**, including their
replacement spans. It does not export standalone DEL/SNP records or reconstruct
phased, whole-sample haplotype walks. Producing a complete Minigraph-Cactus
pangenome or haplotype indexes for Giraffe requires those additional steps.

Format references: [rGFA specification](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md),
[VG conversion options](https://github.com/vgteam/vg/wiki/vg-manpage),
[Minigraph-Cactus pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).
