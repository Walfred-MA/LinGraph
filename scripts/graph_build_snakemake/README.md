# Resumable graph-building workflow

This directory contains a Snakemake **6.15.1** workflow for building local
minimum-set graphs from:

- an assembly list (`-q`);
- either a cohort reference supplied as a FASTA path/list name (`-r`), or a
  template-defined grouped BED (`--bed-grouped`) that needs no global
  reference;
- an optional BED4+ block file (`-b`);
- an intermediate/checkpoint directory (`-G`); and
- the same graph directory for final graph-build artifacts (`-O` is deprecated).

The launcher checks the Snakemake version exactly and always adds
`--rerun-incomplete`. Every long-running downstream program that supports
`--resume` is invoked with it. Commands are run synchronously: a Snakemake job
is not marked complete while a background process is still running.

## Environment and executables

Create the pinned workflow environment:

```bash
conda env create -f scripts/graph_build_snakemake/envs/snakemake-6.15.1.yaml
conda activate minsetref-snakemake-6.15.1
```

The following k-mer executables are expected in `--script-folder` by default:

```text
kmerindex
kmerblocking
kmerpartitionmerger
kmer_searcher
build_local_graphs_parallel_hotspots
```

`kmerblocking` is not required when `-b` is supplied.
`kmerpartitionmerger` is not required when a supplied BED is already grouped.
Paths can be overridden in `config/config.yaml` under `executables`.

The normal minsetref dependencies must also be available: `samtools`,
`minimap2`, `winnowmap`, and, unless `--skip-blastn` is used, `makeblastdb` and
`blastn`. The graph-alignment stage additionally uses the tools resolved by
`local_graph_whole_cigar.py` (`KmerStrd`, `runmakeblastdb`, `runblastn`, and
`graphcigarlight.py`), either beside the minsetref scripts or in their existing
fallback locations/on `PATH`.

## Assembly list and validation

The list has two whitespace-separated columns. Every input FASTA must already
have a valid `samtools faidx` index at `FASTA.fai`; the workflow checks and
reuses it but does not create or replace it. Exactly two columns are required;
a third column specifying an index path is rejected.

```text
CHM13_h1 /assemblies/chm13.fa
HG002_h1 /assemblies/HG002.h1.fa
HG002_h2 /assemblies/HG002.h2.fa
```

Create any missing indexes before launching the workflow:

```bash
samtools faidx /assemblies/chm13.fa
samtools faidx /assemblies/HG002.h1.fa
samtools faidx /assemblies/HG002.h2.fa
```

Relative FASTA paths are resolved from the working directory where the
launcher is invoked, including paths in a separate discovery list. FASTAs must be
uncompressed, and their paths cannot contain whitespace because the downstream
minsetref query-list readers use whitespace-delimited fields.

The workflow enforces:

- names of the form `SAMPLE_hHAPLOTYPE`, such as `HG002_h1`;
- nonempty, valid `samtools faidx` indexes;
- uncompressed FASTA inputs; and
- assembly and FASTA paths that are safe for the downstream command-line tools.

Whether the first FAI contig begins with `SAMPLE#HAPLOTYPE` (normally
`HG002#1#contig-name`) is recorded in the input metadata for diagnostics, but
it is not required. Later contig names are ignored by this diagnostic. Hotspot
files carry their source assembly name separately, so multiple assemblies with
legacy or unprefixed contig names are accepted.

Validation failure stops before k-mer construction and tells the user to run
the assembly-preparation workflow in the configured `preparation_folder`.
The standalone preparation command is documented in
[assembly preparation](../../tools/assembly_preparation.md).

If `-r` is a list name, that entry is moved to the first row of the normalized
list. If it is a FASTA already present in the list, its entry is moved first.
If it is a new FASTA, it is inserted first; a name is derived from a
`SAMPLE#HAPLOTYPE` header or, for a legacy reference, from the filename with
`_h1` appended.

When `-r` is omitted, `--bed-grouped` is required. The normalized assembly
list preserves its input order, no assembly is labeled as a cohort reference,
and `inputs/assemblies.json` records `reference_mode=per_block_template`.

## BED behavior

The input BED is interpreted as zero-based, half-open reference coordinates.
Columns 1-4 and, when present, the last column are retained. Underscores are
removed from column 4. Other punctuation in annotation-derived names (for
example, semicolon-separated gene names) is converted to a readable, bounded
filename-safe name with a short deterministic hash. Repeated normalized
column-4 names mean the BED is already grouped.

Pass `--bed-grouped BED` when the supplied BED is known to enumerate all loci in
each repeat group, including when every group has only one row. In this mode,
the extractor writes exact template-to-block mappings directly and skips the
reference-to-block Minimap2 alignment and its mapping filter.

With `-r`, grouped BED rows retain the historical reference-coordinate BED4+
behavior. Without `-r`, the file must be BED4 or BED5:

```text
template_contig  start  end  group  [template_assembly]
```

BED column 5, when present, must exactly name an entry in `query_paths.txt`.
For BED4, the contig must occur in exactly one input assembly FAI, from which
the template assembly is inferred. If the same contig name occurs in multiple
assemblies, column 5 is required. Each partition's eligible block-template
assemblies are passed to `minsetref_light.py`; its deterministic priority
selection chooses that partition's reference, and all coordinate refinement
and graph construction for that partition use the selected assembly FASTA.
No arbitrary first sample is copied into a fixed-reference catalog.

For example:

```text
NC_060925.1  26099  59613  group_6355
```

becomes group `group6355`. If `group6355` occurs on several rows, those rows
are written to one merged partition and `kmerpartitionmerger` is skipped. A
BED with no repeated group names follows the mapping/filtering/k-mer partition
merging path. With no BED, `kmerblocking` generates blocks and they are named
`g1`, `g2`, and so on.

When `--alternative alternative.fa` is supplied, the workflow streams the
reference and every alternative record into one fixed catalog. Both FASTAs
must have adjacent `.fai` files, and duplicate record names are rejected. BED
rows may refer to reference or alternative contigs. Each alternative record
with no BED row is added as one whole-record block; a record with BED rows uses
only those defined regions. A BED source haplotype supplied by the alternative
FASTA does not have to occur in `-q`; coordinate analysis keeps that locus as
an external template with its original BED coordinates.

## Optional novel-locus discovery

The public `LinGraph.py graph --static-block` option fixes each partition to
its initial block templates. It skips per-partition minset discovery and
nearby-reference refinement while retaining sample alignment and downstream
variant calling. Original template sequence, coordinates, and orientation
are preserved, including supplied alternatives. The backend launcher accepts
the same option (`static_block: true` in direct Snakemake configuration).
It conflicts with `--find-novel-loci` and `--annotate-novels`.
Without this option, existing behavior is unchanged.

Static partitions use ordinary adjusted manifests containing only the original
templates and their materialized reference paths. No novel-query paths or
synthetic minset outputs are produced. A policy file under `GRAPH/inputs/`
tracks switching this option in either direction, without using script
timestamps or invalidating legacy dynamic runs merely because the option was
added. Completed graphs remain reusable through the public CLI.

Use `--find-novel-loci` to run cohort-wide discovery before block construction.
With no value it reuses `-q`; an optional value selects another indexed
assembly list:

```bash
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt \
  --bed-grouped groupblocks.bed \
  --alternative alternative.fa \
  --find-novel-loci \
  -G graph_work -j 64

python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt \
  --find-novel-loci discovery_query_paths.txt \
  -G graph_work -j 64
```

The fixed `reference + alternative` catalog is used by
`KmerRefCleaner --exclude-all`, then by per-assembly Minimap2 and convergent
Winnowmap+BLASTN cleaning. The residual cohort is self-cleaned by
`minsetref.py`; its finalized catalog is `GRAPH/novel_loci.fa`, with restartable
work under `GRAPH/novelloci/`. Exact k-mer residuals are retained there, while
alignment uses compact residual neighborhoods with 30-kb original-source
flanks so final loci can still be lifted. `--novel-locus-jobs INT` divides the
total `-j` CPU budget across concurrent discovery assemblies (default 4).
The finalized loci receive one last exact 31-mer comparison to the combined
reference+alternative catalog. Covered k-mer spans are unioned and
lowercase-masked without changing record coordinates; records with fewer than
1,000 remaining uppercase `A/C/G/T` bases are discarded. The raw BED3 spans,
including spans from discarded records, are written to
`GRAPH/novel_loci.reference_kmer_covered.bed`.
In `--exclude-all` mode, `KmerRefCleaner` uses KmerSearcher's compact bucket
and flat-suffix index rather than the node-heavy count map used by normal
masking-aware mode. Construction is a streaming two-pass read and preserves
rare bucket overflows in an exact side table.

Reference, alternative, and finalized novel records are then concatenated for
k-mer indexing and block construction. Novel records of at most 1,000,000 bp
become one block; only records longer than 1,000,000 bp run through
`kmerindex`/`kmerblocking`. Supplied reference/alternative BED coordinates are
not recomputed. Partition merging still runs only when grouping was not
provided, or when no BED was supplied. `--find-novel-loci` and
`--large-novels` are mutually exclusive.

## Run

Reference by assembly-list name, automatic blocking:

```bash
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 \
  -q query_paths.txt \
  -G cohort_graph_work \
  -j 64
```

Reference by FASTA and user-supplied blocks:

```bash
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r /assemblies/reference.fa \
  -q query_paths.txt \
  -b blocks.bed \
  --bed-grouped \
  -G cohort_graph_work \
  -j 64
```

Template-defined grouped blocks, with no cohort-wide reference:

```bash
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -q query_paths.txt \
  --bed-grouped groupblocks.bed \
  -G cohort_graph_work \
  -j 64
```

Reference-free template mode intentionally rejects alternative/novel reference
stages and novel annotation because those interfaces still require one
cohort-wide reference. The graph summary remains supported by embedding every
partition template.

Inspect the DAG without running commands:

```bash
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt \
  -G cohort_graph_work -j 64 --dry-run
```

To continue directly into cohort variant calling after graph construction,
give a separate calling output directory with `--continue-cohort-call`:

```bash
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt \
  -G cohort_graph_work -j 64 \
  --continue-cohort-call cohort_calls \
  --cohort-call-args "--SVonly 50 --merge-locus-processes 4"
```

The handoff occurs only after the graph Snakemake process exits successfully.
It reuses the graph folder, assembly list, reference name, CPU budget,
Snakemake executable, script folder, and Slurm settings. The cohort workflow
is independently checkpointed, so rerunning the same command resumes both
stages. In graph `--dry-run` mode, the launcher prints the deferred cohort
command without executing it. This continuation requires a named reference
assembly (or a reference FASTA matching exactly one assembly-list row), so it
is unavailable in reference-free template mode.

Extra Snakemake arguments can follow `--`, for example `-- --keep-going`.
To remove stale locks, append `--unlock` (or `-- --unlock`) to the same
launcher command. With `--continue-cohort-call`, this unlocks both workflows
and exits without building graphs or calling variants; graph completion is
not required. The cohort unlock reuses `DIR/inputs/cohort_call.run.json` from
the existing run. Only unlock when no other instance of either workflow is
running. Remove the unlock option to resume normal processing.

Restart with the same command after a failure. Do not delete `-G`; it contains
the Snakemake checkpoints and the native resume work directories.

When every partition batch is already complete but restored files have newer
timestamps than their Snakemake outputs, first verify the restored files and
then use Snakemake's explicit `--touch` mode on `minset_batches_complete`.
Run the normal pipeline command afterward. Breakpoint consistency builds the
cache and blocked BED for each eligible partition, then regenerates
`Graphs.list` from only completed partitions. No partition folder
is deleted.

Cohort local/global/small novel annotation is disabled by default. Add
`--annotate-novels` to run the optional post-graph mapping pass and write
`-G/cohort_annotations/` sidecars. This pass does not rewrite graph,
alignment, or linear coordinates.

The pass is driven by `run_partition_novelty_pipeline.sh --annotations-only`,
which reuses the workflow's `partition_coordinate_analysis/` and
`unified_novel_refinement/` results and runs only the three mapping stages
plus the sidecar writer. Its outputs are all under `-G/cohort_annotations/`:

```text
cohort_annotations/
  partition_novelty/        # map_partition_local_novel.py
  global_novel_mapping/     # map_global_novel.py (coordinates on the reference)
  small_novel_mapping/      # map_small_novel_minset.py (smallnovel.fa minset)
  PARTITION/PARTITION.annotations.tsv
  residual_novel_annotations.tsv
  _ANNOTATIONS_SUCCESS.json
```

`--skip-blastn` is forwarded to the three mapping stages. The small-novel
priority tiers default to the cohort reference haplotype followed by
`HG38_h1,HG002`; pass `--priority` to the driver directly to change them.

Imported alternatives are annotated against their stored catalog records,
using record coordinates and the `alternative` target namespace. Their source
assembly need not be listed in `-q`, and source provenance never substitutes
different assembly sequence for an imported record.

Annotation completion records the size and modification time of every output
table, including per-partition sidecars. The launcher schedules reconstruction
when that inventory is missing or changed; the final refinement/annotation
stage also verifies its outputs on `--resume`.

Each external graph-alignment command has a six-hour (21,600-second) cohort
timeout by default. Override it with `--alignment-timeout SECONDS`; this is
independent of the Slurm job time limit. A partition that reaches this timeout
is skipped without failing its batch, emits a `[WARNING]`, and is upserted in
`GRAPH/Graphs/fail.list`. Its completed bulk BLASTN and Winnowmap batches stay
in the stable
`GRAPH/Graphs/PARTITION/.PARTITION_align.txt.direct.work/` checkpoint
directory. A later forced rerun resumes those completed batches. Non-timeout
command, parsing, and data errors remain fatal.

Bulk fast alignment is now the default. `--fast` and the historical
`--fast-mode` spelling explicitly select the same behavior. BLASTN sees the
complete partition query, while only incompletely aligned repeat-heavy
queries are emitted to Winnowmap as complete-record batches capped at 256 MiB.
Each combined SAM is read once to index byte offsets by query; graph-CIGAR
conversion is then multiprocessed per query, with workers seeking only their
indexed SAM rows. `_linear.txt` is not generated in fast mode unless
`--linear` is explicitly supplied. `--slow-rigorous` runs both BLASTN and
Winnowmap for every sequence. `--max-rigorous` additionally disables the
k-mer redundancy prefilter, starts cleaning directly with Minimap2, and then
uses the same alignment policy as `--slow-rigorous`:

```bash
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt --bed-grouped groupblocks.bed \
  -G graph_work -j 64

# Also generate _linear.txt:
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt --bed-grouped groupblocks.bed \
  -G graph_work -j 64 --fast-mode --linear

# Run both graph aligners for every sequence:
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt --bed-grouped groupblocks.bed \
  -G graph_work -j 64 --slow-rigorous

# Also omit the k-mer redundancy prefilter:
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt --bed-grouped groupblocks.bed \
  -G graph_work -j 64 --max-rigorous
```

## Selective Slurm submission

Pass `--slurm` to submit optional novel-locus discovery, per-assembly k-mer
searches, and per-batch minimum-set jobs through `sbatch`. The wrapper waits
for each submitted job, so downstream Snakemake rules do not start before its
output is complete.
`--slurm-account` and `--slurm-partition` are required when using the dedicated
options; `--slurm-time` defaults to 24 hours. Alternatively, pass arbitrary
shell-quoted `sbatch` options together with `--slurm-args`. Values supplied in
`--slurm-args` override the corresponding dedicated account, time, or partition
option. `--mem`, `--mem-per-cpu`, and `--cpus-per-task` likewise replace the
managed memory and CPU requests rather than being emitted twice. These extra
arguments apply to every Slurm-submitted workflow stage. The text is parsed
into arguments and is not evaluated as shell code.
Defaults are:

- `discover_novel_loci.py`: 32 CPUs and 128G RAM when discovery is enabled
  (adjust with `--novel-slurm-cpus` and `--novel-slurm-memory`);
- `kmer_searcher`: 16 CPUs and 32G RAM per assembly;
- integrated partition processing: 32 CPUs and 64G RAM per batch, running
  four partitions concurrently with eight threads each; and
- at most 20 simultaneously submitted jobs (`--slurm-jobs`).

All other rules remain local. Slurm stdout/stderr files are written under
`GRAPH_FOLDER/slurm_logs/`.

The C++ hotspot FASTA builder writes each partition atomically through a
temporary file. The workflow persistently stores only
`Graphs/PARTITION/PARTITION_samples.header`. Initial block FASTA/BED files,
persistent headers, refined graph files, and alignment outputs all share this
single per-partition directory. The shared `Graphs/` tree is not owned as a
Snakemake `directory()` output, so rerunning the initial block rule cannot
delete downstream graph artifacts. Normal mode defaults to four
partition workers per batch job; fast mode defaults to one. Threads per worker
are calculated as the job's available CPUs divided by its partition workers.
A worker materializes its current partition exactly once
from the indexed assemblies, then keeps that same temporary FASTA through
k-mer cleaning, `.unique.bed` creation, local refinement, final `.FA`
construction, `_align.txt`, and `_linear.txt`. Only after all three final graph
files validate does the worker delete the temporary `_samples.fasta` and its
successful cycle-alignment work. `$SLURM_TMPDIR` is used when available.
Failed cycle work remains available for diagnosis. Sequence reads are streamed
in 16-MiB chunks, so the full cohort of partition sample sequences is never
held in RAM or retained on disk.

```bash
python3 scripts/graph_build_snakemake/run_graph_pipeline.py \
  -r CHM13_h1 -q query_paths.txt \
  --bed-grouped grouped_blocks.bed \
  -G graph_work -j 64 \
  --partition-batches 20 --slurm --slurm-jobs 20 \
  --slurm-args "--time=200:00:00 --account=mchaisso_100 --partition=qcb"
```

With the default 32-CPU minset Slurm allocation, normal mode therefore runs
four partitions with eight threads each, while fast mode runs one partition
with all 32 threads. `--slurm-args --cpus-per-task=N` replaces 32 in this
calculation. `--minset-partition-threads` is intentionally not an input.

Within each partition, that calculated thread budget is propagated to hotspot
FASTA materialization, `KmerRefCleaner`, `KmerStrd`, Minimap2, BLASTN,
Winnowmap, alignment parsing, and graph-CIGAR conversion. `KmerStrd` keeps its
k-mer orientation-learning pass ordered because later records depend on k-mers
learned from earlier records; its independent reverse-complement and header
finalization pass uses `-t/--threads`. Iterative cleaning cycles likewise remain
ordered, while the aligner and per-query work inside each cycle are parallel.

## Stages

The workflow performs the requested stages in this order:

1. validate/reorder inputs and require every input FAI;
2. concatenate reference and alternatives, tagging imported alternatives with
   `sequence_role=imported_alternative`; optionally k-mer-filter,
   alignment-clean, and cohort-self-clean novel loci;
3. concatenate fixed and novel catalogs and build their k-mer index;
4. automatically block or normalize the supplied BED, adding undefined
   alternatives and size-dependent novel blocks;
5. extract/index the block-reference FASTA;
6. map the catalog to blocks, filter mappings, and merge partitions when needed;
7. build local-block folders and the authoritative partition target list;
8. index partition targets and run `kmer_searcher` once per assembly;
9. reduce all hotspot reports into persistent per-partition `.header`
   coordinate manifests without reading assembly sequence bases;
10. divide the active headers into `min(partition_batches, target_count)`
    reusable batches and run four partitions concurrently per batch with eight
    threads per partition. For each partition, in one uninterrupted lifecycle:
    materialize `_samples.fasta` once; run the exact 31-mer `KmerRefCleaner`
    pre-clean and `minsetref_light.py`; validate `.unique.bed`; run local
    coordinate analysis and nearby-reference refinement; build the final
    locally refined `.FA`; align the same `_samples.fasta` to produce the
    payload-bearing `_align.txt`; build `_linear.txt`; validate every output;
    then delete `_samples.fasta` and successful cycle-alignment work. An
    imported alternative always remains its original independent template,
    including when its source haplotype occurs in the cohort or equals the
    primary reference. Its source fields are provenance only: assembly flanks
    cannot expand or replace its sequence. Newly discovered novel templates
    remain eligible for refinement, including in mixed partitions. An imported
    source haplotype absent from `-q` remains usable without that assembly;
11. after every graph batch finishes, run breakpoint consistency over the same
    active partition list. Partitions lacking `_align.txt` are skipped. The
    default uses every sample alignment and every source BED haplotype. A first
    pass learns one cohort-wide Assemblysmall control from all rows; a second
    pass replays that frozen control on every sample and the reference. It
    writes `PARTITION_breakpoint_consistency.initial.tsv`, the corresponding
    blocked BED6, and finally commits `PARTITIONcache.json` containing the
    learned control;
12. atomically regenerate `Graphs.list` without deleting folders.
    Only partitions with a nonempty `.FA`, `_align.txt`, current all-sample
    cache, and completed blocking output enter the final list; exclusions are
    reported in
    `checkpoints/Graphs.excluded.tsv`;
13. run `build_local_graphs.py` and write the self-contained package under
    `-G/summary/`. A present partition cache is packed only after its
    all-sample scope, alignment-row count, blocked-row count, and dependency
    timestamps validate. Missing caches are allowed;
14. with `--annotate-novels`, only after all final local graphs exist, run
    cohort local/global/small novel mapping and write annotation sidecars.
    This optional pass never rewrites `.FA`, `_align.txt`, `_linear.txt`, or
    their coordinates.

`match_partition_blocks.py` is deliberately not part of this graph-build DAG.
Reference/query block matching runs in the downstream
`run_sample_pipeline.py` stage, where `--reference-sample`, its alignment, its
blocked BED, and its FASTA are explicitly defined before GenomeLift.

When neither `--find-novel-loci` nor `--large-novels` is used, the workflow
supplies an empty established large-novel catalog. This keeps the three
required inputs sufficient and means no large novel loci participate in
annotation.

## Outputs

Important intermediates under `-G` include:

```text
inputs/query_paths.normalized.txt
inputs/assemblies.json
inputs/reference_alternatives.fa
inputs/reference_alternatives_novels.fa
novel_loci.fa
novelloci/
main_chroms_withnovel_blocks.named.bed
<REFERENCE>_withnovel_blocks.fasta
<REFERENCE>_withnovel_blocks.fasta_aligned.filtered.bed
<REFERENCE>_withnovel_blocks.merged.bed
Graphs/PARTITION/
kmermaps/
Graphs.template.list      # initial PARTITION.fasta template paths
Graphs.template.list.bin  # KmerSearcher cache for those templates
Graphs.list               # final PARTITION.FA graph paths
checkpoints/Graphs.active.list
partition_batches/
partition_local_adjusted_paths/
local_adjusted_partition_paths/
partition_coordinate_analysis/
unified_novel_refinement/
cohort_annotations/  # only with --annotate-novels
summary/
```

Final per-partition files are written under `-G/Graphs/PARTITION/`:

```text
PARTITION.FA
PARTITION_align.txt
PARTITION_linear.txt
PARTITION_breakpoint_consistency.initial.tsv
PARTITION_breakpoint_consistency.bed
PARTITIONcache.json
```

The alignment and linear files are generated from the same already-refined
`.FA`, so their path coordinates are compatible. Oriented query bases are
attached to every `I` and `X` operation in `PARTITION_align.txt` by default;
pass `--no-alignment-payload` to omit them. The file remains a 10-column table
and linearization continues to consume the operation lengths.
Breakpoint consistency uses all samples by default. The standalone summary,
single-partition, cache-builder, and graph-builder interfaces accept
`--onlyref` to restrict breakpoint derivation and blocking to the named
reference haplotype(s).
The downstream sample pipeline performs RefMatch only after a reference is
explicitly supplied there. Every query is compared with those reference
blocks; graph construction itself does not require or emit a RefMatch table.

Every completed run writes `-G/summary/` containing:

```text
local_graphs.tsv
alternatives.fasta
partition_caches.jsonl
summary.complete
```

`partition_caches.jsonl` is always present and may be empty. Each nonempty
JSON-line record names its partition and stores the exact text of that
partition's validated `PARTITIONcache.json`. `--compressgraph` remains accepted
as a deprecated compatibility flag but is no longer needed.

When a packaged graph path is byte-for-byte identical to a contiguous slice of
the designated reference (including strand and letter case), its bases are not
duplicated in `alternatives.fasta`. `local_graphs.tsv` records the verified
reference coordinate, SHA-256, and exact FASTA header instead. Reconstruction
fetches the slice from the supplied reference and rejects a reference whose
bases do not match the recorded digest. Paths that do not match exactly remain
fully materialized in `alternatives.fasta`.
