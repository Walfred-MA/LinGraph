# Cohort SNP/indel/SV VCFs to rGFA

The default `--gfa-mode rgfa` exports a reference backbone with insertion,
SNP, substitution, and deletion branches and stable coordinates. It includes
all sizes by default (`--all`). Add `--svonly` to include only non-SNP SVs
at or above `--svcutoff` (default 20 bp). This includes qualifying insertions
and deletions; `--insertion-only` remains a separate filter. `--gfa-mode query` retains the previous local
query-flank GFA behavior; that mode is not suitable as an anchored rGFA.

```bash
python scripts/merged_vcf_to_gfa.py \
  -v cohort_calls/cohort.sv.vcf cohort_calls/cohort.indel.vcf cohort_calls/cohort.snp.vcf \
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
is also available in `LinGraph.py ... --MC-graph`. The older `--size-cutoff`
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

`LinGraph.py --MC-graph` passes the selected/default backbone, lifted template
catalog, query assembly list, and `--processes` in graph and singular modes.
With `-G graph/summary`, the default catalog is located in the parent graph
directory. Template selection/lift checkpoints record backbone identity; the
cohort DAG also tracks lifted templates as inputs to the calling stages.

## Sequence checks and large files

The converter reads the VCF body sequentially and retains coordinate/run
metadata and sequence digests, not the sample matrix or sequence payloads.
Query extraction uses scaffold batches and temporary files. Sequence copying
uses bounded chunks. Memory still scales with the number of events, graph
segments, and edges; emitting the full reference and chopping nodes increases
both graph size and memory relative to the legacy local graph.

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

The exporter includes SNP, insertion, deletion, and replacement branches by
default. Use `--insertion-only` to restrict it to insertion alleles. It does not
reconstruct phased, whole-sample haplotype walks. Building haplotype indexes for
downstream tools requires those additional paths and indexing steps.

Format references: [rGFA specification](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md),
[VG conversion options](https://github.com/vgteam/vg/wiki/vg-manpage),
[Minigraph-Cactus pipeline](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md).
