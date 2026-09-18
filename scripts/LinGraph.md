# LinGraph

Use `scripts/LinGraph.py` as the public command-line entry point for assembly
preparation, graph construction, variant calling, merging, and GFA export.
Run the examples from the repository root.
Local execution is the default. Top-level `-h` summarizes the commands, default
alignment modes, and common run options. Put run options after `graph` or
`singular`; use that command's `--help-all` for the complete calling, merging,
construction, SLURM, and GFA tuning options.
Users do not need to invoke the backend scripts directly.

```bash
python3 scripts/LinGraph.py --help
python3 scripts/LinGraph.py graph --help
python3 scripts/LinGraph.py graph --help-all
python3 scripts/LinGraph.py singular --help-all
python3 scripts/LinGraph.py prepare --help
```

## Public controls

| Task | Options |
| --- | --- |
| Variant selection | `--all`, `--svonly`, `--svindel`, `--snp`, `--svcutoff` |
| Local graph alignment | `--fast` (graph only), `--rigorous` |
| Partial runs | `--recall-only`, `--merge-only`, `--gfa-only`, `--vcf-list` |
| GFA export | `--mc-graph`, `--insertion-only [SIZE]`, `--gfa-mode`, `--max-node-length`, anchor and size filters |
| Workflow inspection | `--dry-run`, `--dag-dry-run`, `--snakemake-args`, standalone `--unlock` |
| SLURM | `--slurm`, `--slurm-jobs`, `--slurm-args`, account, partition, time, and memory options |
| Calling and merging | `--samples`, worker counts, per-stage memory, realignment, provenance tags, template fallback, nested insertion and merge thresholds |
| Construction | `--static-block`, partition batches, novelty discovery, alignment timeout, payloads, and per-partition workers |

Use `graph --help-all` or `singular --help-all` for the complete options for that
mode. `--alternative` is a construction FASTA; `--alternative-bed` overrides
the cohort caller's original-template intervals. Settings that cannot apply to
a partial run are rejected. `--sequence-similarity` is retained for compatibility;
the SV backend currently uses its fixed KmerMatch insertion-similarity rule.

Graph mode defaults to `--fast`: bulk BLASTN first, then Winnowmap for qualifying
incomplete, repeat-heavy queries. `--rigorous` / `--slow-rigorous` runs both
aligners for every query. Singular mode always uses rigorous alignment for
both the query and reference; it does not accept `--fast` or `--fast-mode`.
Its added query anchors can otherwise trigger fast-mode fallback repeatedly.
An explicit alignment mode makes a normal graph run check the construction
workflow even when a graph already exists; changing the mode can invalidate
alignment results.

Add `--static-block` to a graph run to use the initial block templates without
calling additional novel loci. New builds still create the initial blocks from
the reference, supplied BED, and any supplied alternative templates. Within
each block, this option skips minset discovery and local template refinement;
the original template coordinates, orientation, and sequence are retained.
Sample alignment, variant calling, and merging still run, with `--all`,
`--svonly`, `--svindel`, or `--snp` controlling the VCF contents as usual.
Graph mode retains both `--fast` and `--rigorous` with static blocks too.

```bash
python scripts/LinGraph.py graph --graph graph -I NA19240trios.txt \
  -O HG38SVs -r HG38_h1 --all --static-block -t 16
```

`--static-block` cannot be combined with `--find-novel-loci` or
`--annotate-novels`. Omitting it preserves the existing discovery/refinement
behavior. A completed graph is reused when no construction is requested;
this flag alone does not rebuild it or remove previously discovered loci.
When construction runs, changes between static and dynamic blocks invalidate
the affected partition results, while repeated runs with the same policy
reuse completed partitions.

In graph mode, `--mc-graph` controls export of `cohort.gfa`. Without it, normal calling and
merging produce VCFs. The local graphs needed for calling remain under `-G`.
`--gfa-only` explicitly exports existing merged VCFs and implies `--mc-graph`.
GFA tuning options require export to be enabled.

Alignment intermediates (`*.align.txt` and `*_align.txt`) now omit `I`/`X`
sequence payloads by default. Graph-CIGAR comparison restores those bases
from the oriented local query and reference sequences for each comparison,
using the same payload writer as hotspot alignment. Existing payload-bearing
files remain readable. Use `--alignment-payload` to retain the previous
storage format. Graph-to-reference outputs still contain payloads; gap
filling fetches FASTA sequences for its new rows and retains that output format.

## Prepare inputs separately

Both modes require prepared, uncompressed assemblies with adjacent indexes
at `FASTA.fai`. The input list has exactly two whitespace-separated columns:

```text
CHM13_h1 /assemblies/prepared/chm13.fa
HG002_h1 /assemblies/prepared/HG002.h1.fa
HG002_h2 /assemblies/prepared/HG002.h2.fa
```

Relative FASTA paths are relative to the list's directory. Names must have
the form `SAMPLE_h1`, `SAMPLE_h2`, etc. FASTA paths cannot contain whitespace.
Do not add an FAI column.

Run preparation explicitly with the standalone tool:

```bash
python3 tools/prepare_assemblies.py \
  -q raw_assemblies.list -O prepared_assemblies -j 4 --contignamefix
```

Use `prepared_assemblies/query_paths.prepared.txt` as LinGraph's `-I` input.
Add `--remask` to the preparation command when a fresh mask is needed.
`LinGraph.py prepare` remains an explicit shortcut to this tool. For name
fixing alone, use `tools/namecontigsfix.py`. See
[assembly preparation](../tools/assembly_preparation.md) for examples.

Before starting either mode, LinGraph checks each assembly's **first sequence
only**: its identifier must be `SAMPLE#HAPLOTYPE` or start with
`SAMPLE#HAPLOTYPE#`, it must contain lowercase `a/c/g/t`, and the first entry
of `FASTA.fai` must match the header and sequence offset. This includes the
reference assembly, except that **`CHM13_h1` and `HG38_h1` skip the contig-prefix
and masking checks**. These two references can retain native names such as
`NC_060925.1` and `chr1`; their FASTAs and matching adjacent indexes are still
required. Later sequences are not scanned by this preparation check.
Preparation checks apply to the selected samples and chosen reference;
unselected assemblies in the saved graph cohort are not checked for masking or
sample-prefixed contig names. The main pipeline never masks, renames, or creates an assembly index.
A failed check stops before computation and prints a command using
`tools/prepare_assemblies.py`. The standalone `run_sample_pipeline.py` uses
the same read-only checks and also requires inputs to be prepared separately.

## Graph mode

```bash
python3 scripts/LinGraph.py graph \
  -I prepared_assemblies/query_paths.prepared.txt \
  -G cohort_graph -O cohort_calls -t 16
```

The first assembly is the reference by default; select another with
`-r CHM13_h1`. A reference FASTA path is also accepted. Its sample name is
inferred from its first sequence prefix, or supplied with `--reference-name`.

LinGraph builds/resumes an empty or incomplete `-G`, calls every cohort
assembly, and always merges their VCFs. For a completed graph with a saved
cohort list, `-I` can be omitted. Supplying construction inputs invokes the
existing Snakemake workflow so it can check its pending work. When resuming a
build, repeat the original construction options, including BED and novel-locus
settings. The graph-building backend requires at least four threads.

Optional construction settings:

```bash
python3 scripts/LinGraph.py graph -I prepared_assemblies/query_paths.prepared.txt \
  -G cohort_graph -O cohort_calls -r CHM13_h1 \
  -b grouped_blocks.bed --bed-grouped --find-novel-loci --MC-graph
```

`-L partitions.list` restricts calling to those partitions. Omit it to call
all active graph partitions. Construction still completes the graph workflow.
The default merge is `--svonly`, writing `cohort.sv.vcf`. With `--mc-graph`
(or `--MC-graph`), the default becomes `--all`, writing `cohort.sv.vcf`,
`cohort.indel.vcf`, and `cohort.snp.vcf` and converting them to `cohort.gfa`. Explicit modes override these defaults:
`--snp` writes only `cohort.snp.vcf`; `--svindel` writes `cohort.sv.vcf` and
`cohort.indel.vcf`, splitting the same merged SV calls at `--svcutoff` (default
20 bp). `--all` includes SNPs and SVs of every size in these three separate VCFs.

Individual calling extracts SNPs from encoded insertion alignments by default,
including short insertions; `--svonly` suppresses SNPs. Existing VCFs must be
recalled to add these SNPs. Use `--recall-only --all` to require existing
alignments and regenerate only the VCFs and merged output (see below).

The GFA converter receives `--svonly` when the main run uses `--svonly`, and
`--all` otherwise. To merge existing calls without graph building, calling,
or GFA export, use:

```bash
python scripts/LinGraph.py graph -G graph -O cohort_calls --merge-only --all -t 16
```

### Recall individual VCFs and merge, using existing alignments

```bash
python scripts/LinGraph.py graph -G graph -O HG38SVs -r HG38_h1 \
  --recall-only --all -t 16
```

This reads `HG38SVs/inputs/cohort_call.run.json` and directly runs the individual
VCF caller on each saved `samples/NAME/NAME.graphcigartoreffix.tsv`, then merges
the recalled VCFs into `cohort.sv.vcf`, `cohort.indel.vcf`, and `cohort.snp.vcf`
under `HG38SVs/`. No Snakemake DAG is launched.
Template discovery, graph alignment, gap filling, graph completeness checks,
and preparation scans are skipped. A missing required input stops the run
before any individual VCF is changed.

The saved reference, FASTA paths, templates, coordinate maps, and biological
calling settings are reused. `-r`, if supplied, must match that reference.
`-I samples.list` can select a subset of saved samples using the same FASTA
paths; the saved reference is always included. Omit `-I` for all samples in
the saved calling run. `--svonly`, `--snp`, and `--svindel` also work.
`--samples` can select saved sample names. Calling overrides such as
`--no-realignment`, `--edge-blackregion`, `--no-pa-tag`, and
`--separate-adjacent-indels`, plus merge thresholds, can be supplied explicitly;
omitted values retain the saved settings.

The first recall replaces individual VCFs only after successful calling.
Repeating the command resumes completed samples using checkpoints under
`OUTPUT/checkpoints/vcf_recall/`; logs are in `OUTPUT/logs/vcf_recall/`.
`--dry-run` prints the commands without writing outputs. The original workflow
configs are not modified, and `--merge-only` and `--recall-only` are mutually
exclusive. Run only one workflow against an output directory at a time.

Add `--mc-graph` to a merge-only or recall-only run to export its merged result.
To export an existing merged VCF by itself:

```bash
python3 scripts/LinGraph.py graph -G graph -O HG38SVs -r HG38_h1 --gfa-only --all -t 16
```

`--merge-only --vcf-list vcfs.list` selects an explicit list of individual VCFs.

The cohort calling workflow ignores script modification times when deciding
whether to reuse completed stages, including template discovery and lifting.
Switching `--svonly` to `--all` changes VCF calling and merging settings; it does
not invalidate templates or alignments. Changed data or stage settings, missing
outputs, and incomplete jobs can still schedule upstream work. To apply a script
update to completed outputs, explicitly force the affected rule or change its
stage protocol. `--forcerun graphreftovcf_sample` does not restrict the upstream
DAG.
Without `-r`, export defaults to `GRAPH/inputs/reference_alternatives_novels.fa`.
With `-r NAME` or `-r FASTA`, that selected haplotype is the calling and GFA
backbone. Blocks without an accepted alignment for that haplotype supply their
reference-class templates as alternatives. Imported alternatives retain their
original sequence and bounds and bypass refinement; source sample names are
provenance only. Newly discovered novel templates remain refinable. Before
loading, imported templates are placed using their own sequence, while novel
templates use source flanks. Unresolved placements retain their sequence as
unplaced paths. HG38_h1 is restricted to chr1–22/X/Y, including for alternative sources.
`GRAPH/summary/alternatives.fasta` is passed as a lookup-only catalog to resolve
local target coordinates; its duplicated records are not exported as extra
paths. Coordinate translation requires matching source metadata and sequence.
LinGraph writes the query assembly list and passes the requested worker count.
See [merged_vcf_to_gfa.md](merged_vcf_to_gfa.md) for coordinates and VG checks.

The finalized valid partition list is `GRAPH/Graphs.list`. Cohort calling does
not require a final KmerSearcher cache. Singular calling searches the original
`PARTITION.fasta` templates through `GRAPH/Graphs.template.list` and its
`GRAPH/Graphs.template.list.bin` index; graph alignment still uses the refined
`PARTITION.FA` files. Refined graph sequences are never substituted for the
original hotspot templates. Reconstruction creates the templates by default, with
`--original_fasta` available when their chromosome source differs from the
graph-path reference. When a calling subset/order differs from an existing
build-time template list, a separate `.selected.template.list[.bin]` preserves
that original list and keeps hotspot indices aligned with the calling list.
Legacy `Graphs.list.bin` files are not used for hotspot discovery. Singular
startup uses the supplied graph and template index directly: it does not scan
partition `.FA`/`.fasta` files or validate/rebuild the k-mer index. Index creation
belongs to graph reconstruction or an explicit `graph_list_cache.py` run.
For a custom subset passed as `-L selected.list`, supply the matching
`selected.template.list` and `selected.template.list.bin` in the same partition
order. No separate dummy-query index-building stage runs during singular startup.
KmerSearcher retains its normal missing/legacy-cache handling when it actually
searches for hotspots.

Resume checks track the graph list, template list, index and `graphs.complete`
marker, rather than individual partition files. After changing partition data,
refresh the graph package/index before calling again. Existing sample
checkpoints are upgraded without rerunning completed stages when their
settings, stage inputs and outputs still match and the graph manifests have
not changed since completion. This upgrade does not inspect partition files.

Reference-cache completeness is determined only by the presence of nonempty
`NAME.align.txt` and `NAME.align.txt_blocks.bed` files in the selected reference
cache directory and a matching `signature.fai_lengths_sha256` in
`reference_cache.json`. This hash compares sorted reference contig names and
lengths, ignoring FAI row order and offsets. Completed caches are reused without
checking other signatures, the completed-stage list, or hotspot files. Missing
or empty outputs use the existing preparation/resume logic. If any nonempty
output exists but the FAI-length signature is missing or mismatched, LinGraph
stops without overwriting the cache. Use a separate reference name for a
different reference version, or restore the matching metadata.
`--reference-caches DIR` continues to select supplied reference results without
validation.
The graph alignment list is also copied to
`GRAPH/summary/Graphs.list`; LinGraph refreshes this package copy when
reusing a completed graph. Entries may be bare names (`group1`) or paths
(`group1/group1.FA`, including absolute paths). Readers take the final path
component and strip a recognized suffix; dots inside names are preserved.
For an existing `Graphs/` directory, singular mode can use the summary copy if
the graph-root list is absent. It does not discover partitions by scanning their
FASTA files.

## Singular mode

**Singular mode calls samples independently and accepts one or many samples
in one run.** Each haplotype assembly gets its own VCF and coverage report,
using the same existing graph and reference alignments. `-i` supplies one
haplotype FASTA; `-I` supplies a list of haplotypes from one or multiple samples.

One assembly:

```bash
python3 scripts/LinGraph.py singular \
  -i /assemblies/prepared/HG003.h1.fa --sample HG003_h1 \
  -G cohort_graph -O query_calls
```

Several assemblies, called independently:

```bash
python3 scripts/LinGraph.py singular -I prepared_queries.list \
  -G cohort_graph -O query_calls
```

Singular mode requires an existing graph. It aligns the reference once and
reuses that alignment for the queries. It produces one VCF and coverage report
per input assembly. With only one sample (from `-i` or a one-entry `-I` list),
calling skips cohort merging, including when `--all` or `--merge` is supplied;
the result is `OUTPUT/samples/NAME/NAME.vcf`. With multiple samples, add `--merge`
to also merge the independent calls. Explicit variant-selection options
(`--all`, `--svonly`, `--svindel`, or `--snp`) also enable merging with multiple
inputs. Without those options, only separate calls are produced. Singular mode does not
export GFA, regardless of the number of samples; `--mc-graph`/`--MC-graph`,
`--gfa-only`, and GFA export settings are available only in graph mode.
Explicit `--merge-only` runs still merge existing VCFs.
The coverage report audits the selected reference assembly's FAI contigs.
Independent loci declared by `##alternativeLocus` are excluded from reference
coverage when absent from that FAI; their query alignments still count toward
query coverage. Their source coordinates describe provenance and are not used
as placements on the selected reference.

A compact graph means the pipeline's **summary directory**, containing
`local_graphs.tsv` and `alternatives.fasta` (and, when available,
`partition_caches.jsonl`). It is not an arbitrary compressed archive. Pass
either a graph root containing `summary/` or the summary directory itself.
Summary-only packages without a graph list/layout are reconstructed into the
output directory. Existing graphs are used as supplied, without checking for or
repairing missing individual FASTAs:

```bash
python3 scripts/LinGraph.py singular -I prepared_queries.list \
  -G graph_summary -r /assemblies/prepared/chm13.fa \
  -O query_calls --all
```

The current exporter stores all path bases, including main-reference paths,
in `alternatives.fasta`. It reuses final partition `.FA` files and original
`.fasta` templates when available, following the finalized `Graphs.list`.
Only partitions with missing FASTAs use producer-side sequence extraction.
Older reference-slice packages remain readable.

Reconstruction uses the summary package and the supplied graph-build reference
FASTA; it never reads cohort assembly sequences. Non-reference original
templates are packaged alongside graph paths. Older packages can recover an
original interval from stored paths covering the same haplotype and contig,
including reverse strands and overlapping/adjacent pieces. A package that
does not contain all required bases is reported as incomplete; the exporter
must include them before distributing it. `--reference-haplotype` identifies
the supplied reference when the package's omitted reference rows are ambiguous.
Omit `-L` to process all partitions; a supplied list may contain partition
names or paths, one per row.

For default-reference `--MC-graph`, retain the original graph's
`inputs/reference_alternatives_novels.fa`. Passing `-G graph/summary` finds it
in `graph/`. An explicit `-r` uses the selected assembly instead. Supply original
source assemblies for template lifting; templates with unavailable sources are
retained as unplaced paths.

Coverage reports describe the reference/query bases represented by the call
and the uncovered regions; they are not sequencing-depth reports. An optional
`--satellite-bed satellites.bed` adds satellite annotations. Without that BED,
the satellite category is unannotated (reported as zero), with those regions
included in the remaining coverage categories.

For GFA conversion in graph mode, the summary FASTA is detected as a coordinate lookup only.
The calling run's lifted template catalog supplies fallback alternatives.
`--alternative` is a graph construction input available only in graph mode.
GFA conversion includes SNPs, indels, and SVs by default. Add
`--insertion-only` to retain insertion-only export with a 50 bp minimum, or
`--insertion-only SIZE` to choose another minimum.

## Outputs and restarting

```text
OUTPUT/
  samples/NAME/NAME.vcf
  samples/NAME/NAME.coverage.summary.tsv        # singular mode
  samples/NAME/NAME.coverage.missing_reference.bed
  samples/NAME/NAME.coverage.missing_query.bed
  cohort.snp.vcf                             # --snp or --all
  cohort.indel.vcf                           # below-cutoff SVs with --svindel or --all
  cohort.sv.vcf                              # --svonly (default), --svindel, or --all
  cohort.gfa                                  # --MC-graph
  lingraph/run.json
  lingraph/logs/
  lingraph/checkpoints/
  lingraph/reconstructed/Graphs/               # when reconstruction is needed
```

Keep `-G` and `-O` separate and non-nested. Repeat the same command after a
failure. Graph/cohort workflows retain their Snakemake checkpoints. Singular
stages reuse successful outputs only when their recorded commands, input
sizes/timestamps, and output sizes/timestamps still match. Only VCFs from the
current input list enter the singular merge; old samples in `-O` are excluded.

After an interrupted run has stopped, clear stale workflow locks with:

```bash
python scripts/LinGraph.py --unlock
```

This clears `.snakemake/locks` in this repository's `graph_build_snakemake/`
and `cohort_call_snakemake/` directories, then exits. It requires no run mode,
`-G`, `-O`, assembly inputs, or Snakemake executable. Checkpoints, incomplete-job
records, and outputs are preserved. Use `--unlock --dry-run` to preview the
directories. Unlocking does not stop running jobs; use it after those jobs exit.

Use `--dry-run` to validate prepared inputs and print the planned commands.
It does not create files, launch genomics tools, or submit jobs. For an empty
graph, later commands are shown with the paths that construction will create;
this is an orchestration preview, not a Snakemake DAG validation.

For the actual Snakemake job plan on an existing completed graph:

```bash
python3 scripts/LinGraph.py graph -G graph -O HG38SVs -r HG38_h1 --all -t 16 --dag-dry-run
```

This refreshes the saved cohort list and stage configurations, then runs the
backend launchers with Snakemake's dry-run option. It runs no workflow jobs.
Because it updates configuration files, use it when that output directory has
no active workflow. Ordinary `--dry-run` remains the preview that writes no files.
Additional calling-workflow controls can be supplied through LinGraph, for example
`--snakemake-args '--forcerun graphreftovcf_sample'`. A forced rule still includes
its required upstream dependencies.

## SLURM and dependencies

```bash
python scripts/LinGraph.py graph --graph graph -O HG38SVs -r HG38_h1 \
  --all --rigorous -t 16 --slurm --slurm-jobs 40 \
  --slurm-args '--time=20:00:00 --account=mchaisso_100 --partition=qcb' \
  --dry-run
```

`-slurm` is an alias for `--slurm`. Normal graph runs submit the existing
workflow stages as separate SLURM jobs, with at most `--slurm-jobs` simultaneous
wrapper jobs (default 20). Snakemake's scheduling capacity is raised to allow
that limit independently of `-t`; the example allows up to 40 jobs while
retaining 16-thread calling workers. Build jobs retain their stage CPU defaults.
The controller stays running and waits for the workflows. Calling job logs are
under `OUTPUT/slurm_logs`; construction logs are under the graph directory.

`--slurm-args` accepts a shell-quoted list of sbatch options and forwards it to
each submitted job. It is parsed as arguments, not evaluated by a shell.
Explicit sbatch resource options override stage defaults. You can also use
`--slurm-account`, `--slurm-partition`, `--slurm-time`, and `--slurm-memory`.
Per-stage controls such as `--graphcigar-memory`, `--vcf-memory`, and
`--merge-memory` appear in `--help-all`. Optional GFA export gets its own SLURM job.

Singular and partial runs (`--merge-only`, `--recall-only`, `--gfa-only`) use one
allocation and also accept `--slurm-args`. `--slurm-jobs` applies to full graph
workflows and is rejected for partial runs. Their launch script and
`slurm-JOBID.log` are under `OUTPUT/lingraph`. Activate the required environment
before launching; binaries and input paths must be available on compute nodes.

Graph mode uses the existing Snakemake **6.15.1** environment and genomics
executables documented in [the graph workflow guide](graph_build_snakemake/README.md).
Singular mode uses the existing KmerSearcher, graph alignment, matching,
liftover, and VCF programs. The interface does not install dependencies.
