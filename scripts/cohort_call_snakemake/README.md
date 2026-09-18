# Merge modes and outputs

The default is `--svonly`. With `LinGraph.py --mc-graph`, the default is `--all`.
An explicit mode overrides that default:

| Mode | Output |
| --- | --- |
| `--all` | `cohort.sv.vcf`, `cohort.indel.vcf`, and `cohort.snp.vcf` |
| `--svonly` | `cohort.sv.vcf`: merged SV representatives >= cutoff |
| `--snp` | `cohort.snp.vcf`: SNPs only |
| `--svindel` | `cohort.sv.vcf` and `cohort.indel.vcf`: merged SVs split by size |

`--svcutoff` defaults to 20 bp (`--minsvsize` remains an alias in the cohort
launcher). Indels in `--svindel` are the SV merger's below-cutoff results,
including nested calls. They use the same merging algorithm as larger SVs.
`--all` saves the three categories in separate VCFs; it does not produce a combined VCF.
The historical `--SVonly [SIZE]` spelling is still accepted by the cohort launcher.

SNPs merge independently, chromosome by chromosome, using exact CHROM/POS/REF/ALT.
In `--all`, SV merging saves SNPs from the existing member-to-representative
insertion alignments, including short insertions and emitted recursive insertion
children. This SNP extraction does not require `--var-in-insert` for top-level
insertions. Original SNP chromosome jobs can run alongside SV merging. The final
SNP concat waits for the SV insertion manifest, then loads, sorts and appends one
insertion locus at a time in memory. These distinct contigs do not create
chromosome jobs or external sort runs. `--snp` skips SV merging and merges only
SNPs already present in the input VCFs. `--svonly` and `--svindel` do not create
these insertion SNP spools.

Production SNP scans write compact binary observations, separately from the SV
record shards. Each SNP observation occupies **21 bytes** (`align=False`):

| Field | NumPy storage |
| --- | --- |
| Reference/template position (VCF, 1-based) | `uint32` |
| REF and ALT (two IUPAC codes) | `uint8` |
| Assembly query coordinate (existing traversal-boundary convention) | `uint32` |
| Query dictionary index | `uint32` |
| `ALLELENAME` dictionary index | `uint32` |
| `LABEL_H` distance or joined-label dictionary index | `uint32` |

The query dictionary stores sample, contig, strand, SNP type, call state,
filter status and the `LABEL_H` encoding mode (missing, scalar or joined).
Allele names have a separate interned list. Scalar `LABEL_H` uses reversible
signed encoding inside `uint32`, preserving negative distances; missing values
remain missing. Cross-graph values such as `25980H&0H` use a separate interned
list, indexed by the same `uint32` field. Contributor order and missing components
are preserved, and these indices are remapped across shards before sorting.
Earlier scalar-only compact shards remain readable. Initial scanning uses `-t` workers.
Each worker reads its assigned VCFs once and writes a single `records.bin` for
SNPs, independently of the existing SV shards. Chromosome offsets/counts and
coverage intervals live in metadata, outside the 21-byte records. In `--all`,
the same scan writes both categories; SNP merging never rescans the input VCFs.

SV chromosome scheduling and resource allocations are unchanged. After the
shared scan, the SNP stage launches one job per original input chromosome.
Each Slurm SNP chromosome job requests **32 CPUs and 64G memory**, independently
of the SV merge settings. Local jobs request 32 threads, scaled down by
Snakemake's available cores. Each job reads all compact observations for its
chromosome into RAM, remaps the provenance dictionaries and sorts numerically
with Polars using the requested threads. Shards remain 21 bytes per observation;
no sorted binary runs or pair-merge rounds are created. Sorting also needs
workspace and dictionaries beyond the raw record size, so 64G is an allocation,
not a guarantee that every cohort size fits. VCF formatting currently uses one
writer; logs report loading, sort completion and VCF-writing progress.

Install the added dependency in the environment running the SNP merger with
`python -m pip install 'polars>=1.0'`. It is also listed by the repository's
installer, currently located at `notused/install.py`. Manual chromosome jobs
should pass `graphvcfmerge_snp.py --stage chrom ... --processes 32` and request
`sbatch --cpus-per-task=32 --mem=64G`. Changing the worker count does not alter
the original compact shards. Interrupted RAM sorts restart from those shards;
existing disk-sort checkpoints are retained but are no longer used.

Saved resume Snakefiles retain their old commands. After stopping the previous
controller and its unfinished SNP jobs, make a copy with just the updated SNP
chromosome rule:

```bash
python tools/upgrade_snp_merge_workflow.py "$WORK"
snakemake graphvcfmerge_small_cohort \
  --snakefile "$WORK/Snakefile.snp-ram" \
  --config "run_config=$WORK/run.json" \
  --cores 64 --resources slurm_jobs=40 \
  --rerun-incomplete --keep-going --notemp --printshellcmds \
  --allowed-rules graphvcfmerge_small_chrom graphvcfmerge_small_cohort
```

This command requires the existing SNP manifest/shards and any insertion
manifest required by the saved concat rule. It schedules unfinished SNP parts
and concat, preserving completed parts and source shards; it does not rerun
scan or SV jobs. The original `Snakefile` and `run.json` are retained. The new
copy refuses to overwrite an existing `Snakefile.snp-ram`.

For manual completion, `graphvcfmerge_snp.py --stage concat --shards-dir PATH
--insertion-snps SV_SHARDS/insertion_snps --snp-only -o OUTPUT_PREFIX` combines
the saved original chromosome parts and appends insertion SNPs to
`OUTPUT_PREFIX.snp.vcf`. It adds insertion contig definitions to the header and
preserves sample order, reference/missing calls, and contributor labels. It
retains all input shards and parts. Repeating concat replaces the output
atomically rather than appending duplicate rows; older manifests that already
merged insertion sources into chromosome parts remain readable without adding
those sources twice. `--snp-only` at concat leaves existing indel VCFs untouched.
Passing `--insertion-snps` at scan/prepare instead merely saves its path for
concat; the insertion directory need not exist yet.

SV concat accepts `graphvcfmerge.py --sv --stage concat --shards-dir PATH
-o cohort.sv.vcf --output-small --insertion-snps PATH/insertion_snps -t 8`.
For uncompressed output, workers copy chromosome parts directly to separate
byte ranges in a temporary output, retaining chromosome and sample order.
Each output is replaced atomically after its workers finish; source parts are
retained. Insertion IDs are collected during copying, so finalizing the SNP
manifest does not reread the merged VCFs or parse sample columns. This requires
the insertion SNP spools already saved by the SV chromosome jobs; concat cannot
recreate missing spools. Start with 4–8 workers on shared storage, where disk
throughput can limit scaling. Gzip output uses one writer. The main and newly
generated resume workflows allocate up to eight concat CPUs, bounded by
`merge_processes`; existing saved Snakefiles retain their prior allocation.

`tmp/resume_unfinished_graphcigar.py` uses these same staged workflow rules,
including separate SV and SNP jobs under Slurm. Existing final files are
preserved. On success it installs missing final VCFs and removes its temporary
work directory; failed jobs retain their working files. Repeating the same
resume command reuses a failed workspace when its inputs, settings, and scripts
are unchanged. It resumes incomplete jobs without repeating completed scans or
SV merges. Changed inputs or settings start a new workspace.

The SV `--insertion-snps DIR` option saves these temporary sources, and the SNP
merger accepts the same option to consume the completed source manifest.
Chromosome completion checks include the saved insertion SNP files, and the
normal merge cleanup removes them after successful cohort publication.

Individual VCF calling also extracts SNPs from mismatches in encoded insertion
alignments by default, including short insertions. These calls use the insertion
template's coordinates and retain sample assembly coordinates and provenance.
`--svonly` suppresses them along with other SNPs. This does not require
`--var-in-insert`, which controls recursive SV encoding/merging.
Different ALT alleles remain separate. Duplicate observations are removed,
while assembly coordinates and PA provenance are retained. `--snp` does not run
SV clustering or alignment. Missing samples get `0` only when reference coverage
is known, otherwise `.`; explicit missing calls remain missing.

## Merge existing calls only

From the directory containing `graph/` and `cohort_calls/`:

```bash
python /path/to/LinGraph/scripts/LinGraph.py graph -G graph -O cohort_calls --merge-only --all -t 16
```

This reads existing per-sample VCFs and writes separate `cohort.sv.vcf`, `cohort.indel.vcf`, and `cohort.snp.vcf` files in `cohort_calls/`.
It does not rebuild the graph, recall samples, or export GFA. An existing
`cohort_calls/inputs/mergevcf.list` is preferred; otherwise saved sample metadata
or `samples/NAME/NAME.vcf` supplies the inputs. No reference FASTA is needed for
this merge-only run. Missing listed inputs cause an error.

To obtain insertion SNPs missing from older individual VCFs using existing
alignments, use `--recall-only --all` as described below.

The direct equivalent is:

```bash
python scripts/cohort_vcf_merge.py run -O cohort_calls --all -t 16
```

Merge scratch directories are deleted automatically after all final cohort VCFs
are successfully written, for local and Slurm runs. Failed merges retain scratch
for resumption. Add `--keep-merge-tmpdir` to retain it after success as well.

The merge-only runner uses `OUTPUT/tmp/merge_only/INPUTS_AND_SETTINGS_SIGNATURE`.
Normal cohort merges use `OUTPUT/tmp/graphvcfmerge/SETTINGS_SIGNATURE`,
`OUTPUT/tmp/graphvcfmerge_small`, and `OUTPUT/tmp/cohort_merge`. Cleanup removes
the current merge's directories; other run directories and non-merge files
under `tmp/` are left alone. An explicit `--merge-tmpdir` is also removed after
success unless `--keep-merge-tmpdir` is supplied. The workflow keeps a small
chromosome manifest in `OUTPUT/checkpoints/` so cleanup does not trigger a
completed merge to rerun.

Script modification times do not invalidate completed cohort calling stages,
including template discovery and lifting. Copying or touching scripts therefore
does not restart templates and alignments. Switching `--svonly` to `--all`
updates VCF calling and merging settings only. Data and stage configuration
dependencies still apply, and missing outputs or incomplete jobs still rerun.
To apply a script update to completed results, explicitly force the affected
rule or update its stage protocol.

GFA export defaults to `--all`. If the main run selects `--svonly`, LinGraph passes
`--svonly --svcutoff SIZE` to the converter; other modes pass `--all` and their
selected VCFs. `--insertion-only [SIZE]` is an additional export filter, with a
50 bp default when the flag is supplied without a size.

## Recall existing alignments and merge only

```bash
python /path/to/LinGraph/scripts/LinGraph.py graph -G graph -O HG38SVs \
  -r HG38_h1 --recall-only --all -t 16
```

The cohort launcher also accepts this mode:

```bash
python scripts/cohort_call_snakemake/run_cohort_call_pipeline.py \
  -G graph -O HG38SVs -r HG38_h1 --recall-only --all -j 16
```

This launches the individual VCF caller directly on the saved gap-filled
graph CIGARs and merges the results. It never launches Snakemake, discovers or
lifts templates, aligns graphs, or fills gaps. Every required CIGAR, coordinate
map, FASTA, index, and template file must already exist; missing inputs fail
before any calls are replaced. Stop any other workflow using this output
directory first.

The saved `OUTPUT/inputs/cohort_call.run.json` supplies the reference, assembly
paths and selected samples; its VCF and merge stage configs supply the existing
calling settings. The mode, `--svcutoff`, and CPU budget can change. No upstream
configs are rewritten. A reference change or new sample requires a normal run.
Use LinGraph `-I samples.list` or the cohort launcher's `--samples samples.list`
to select saved samples; the reference is always included. The cohort launcher
uses the saved assembly paths, so `-q` and extra Snakemake arguments are not
accepted in recall-only mode.

Completed samples resume from `OUTPUT/checkpoints/vcf_recall/`; an unsuccessful
recall preserves that sample's previous VCF. Logs are in
`OUTPUT/logs/vcf_recall/`. Add `--dry-run` to print the plan without changing
files. This mode is distinct from `--merge-only`, which reads existing VCFs
without recalling them. Normal workflow runs, even with `--forcerun
graphreftovcf_sample`, may schedule upstream jobs when their dependencies change.

# Reference selection and fallback templates

`-r SAMPLE_h1` selects any haplotype in the cohort assembly list as the calling
backbone; omitting `-r` selects the first assembly. For local blocks with no
accepted alignment for that haplotype, the workflow selects reference-class
templates, then runs `lift_local_templates.py` before loading them into calling.
Imported alternatives marked `sequence_role=imported_alternative` are fixed
independent sequences: their source labels are provenance, and their original
bounds and bytes survive graph refinement unchanged. Newly discovered novel
templates can still be refined. The lift aligns the actual imported sequence
to the selected backbone; it uses 10 kb source flanks only for refinable
templates. Both modes preserve template sequence, orientation, IDs, and local
coordinates.

`checkpoints/local_reference_templates.unlifted.fa` stores the selected templates;
`checkpoints/local_reference_templates.fa` stores their fresh placement metadata.
Unresolved, ambiguous, and disagreeing placements retain the template as an
unplaced path. One exact supported boundary creates only that attachment. The
exporter uses these placements for graph links without relabeling novel bases
as reference bases. HG38_h1 contributes only chr1–22/X/Y; chrM and non-primary
scaffolds are excluded from backbone and alternative use.

`LinGraph.py --MC-graph` uses the selected haplotype FASTA when `-r` is supplied.
Otherwise its GFA reference defaults to
`GRAPH/inputs/reference_alternatives_novels.fa`. The local summary catalog is
lookup-only and is never exported wholesale. See `../merged_vcf_to_gfa.md` for
manual conversion and VG validation commands.

# Merge completion and resume

A completed graph rebuild invalidates cohort matching and extraction, then
liftover, template selection/lifting, variant calling, and merging. Old calls
against expanded imported templates must be recalled; changing VCF headers
alone does not restore the original template sequence or its missing calls.

To remove a stale cohort-workflow lock after confirming no other cohort
Snakemake controller is running, only the existing output directory is needed:

```bash
python scripts/cohort_call_snakemake/run_cohort_call_pipeline.py --unlock -O cohort_calls
```

This reads `cohort_calls/inputs/cohort_call.run.json`, asks Snakemake to unlock
the cohort workflow directory, and exits without preparing or running jobs.

The split Slurm merge now tracks each chromosome independently. Rerun the same
`run_cohort_call_pipeline.py ... --slurm ... -- --keep-going graphvcfmerge_cohort`
command with the same output directory and merge settings. The scheduler skips
checked chromosomes, submits unfinished ones, and enables final concatenation
only after every chromosome has checked outputs. It also rechecks chromosomes
when an older aggregate `chroms.complete.tsv` exists but individual markers are
missing. Earlier pipeline outputs continue to use their existing Snakemake
completion records.

Under the merge directory (default `cohort_calls/tmp/graphvcfmerge`):

| File | Meaning |
| --- | --- |
| `chroms/<hashed-chrom>/merge.done` | Completed final `.part` and `.promoted.tsv`; skip this chromosome. |
| `checkpoints/<hashed-chrom>/refined.ready` | Top-level grouping finished; reuse saved groups and member CIGARs. |
| `checkpoints/<hashed-chrom>/nested.ready` | Top-level rows and nested candidates saved; start directly at nested calling. |
| `chroms.complete.tsv` | The scheduler verified completion of all chromosomes; final concatenation can run. |

The merger writes these markers automatically, after the files they describe.
Restart files stay on shared storage, outside the old temporary spool directory.
CPU and memory allocations may change on restart. Saved scan metadata and merge
settings must still match. Changed outputs or incompatible checkpoints cause an
error instead of silently treating them as complete.

Nested resume restarts the **whole nested stage**, using an immutable copy of
the top-level rows to prevent duplicate appended rows. It does not resume the
exact insertion path that was running. In multi-core runs, paths coordinate
outside a shared alignment pool. Independent alignments from every active path
use that pool, including work from the last remaining path. Progress logs
include path completions, active/queued alignments, and pool utilization.

The pool has `--processes` worker slots (the chromosome job's CPU allocation).
Final member CIGARs, independent later-group comparisons, and KmerMatch timeout
fallback tasks share these slots. Online grouping decisions still follow their
original order; a single pairwise alignment is not divided across workers.
The queue admits at most twice the worker count in unfinished tasks and 64 MiB
of estimated input payload, with one oversized task admitted alone when needed.
This bounds submitted payloads, not total chromosome or native-aligner memory.
Small alignment batches and worker reuse limit sequence-copy and startup costs.
Worker request files use Python's temporary directory (`TMPDIR` when set);
durable checkpoints and ordered path-result spools stay in the merge directory.

The two-hour alignment-first cluster deadline includes queue waiting. Expiry
cancels that cluster's queued tasks and kills/reaps only its running worker
process groups before entering KmerMatch fallback; other paths keep running.
The affected slots are replaced as needed. Errors other than deadline expiry
still fail the merge. Shutdown also cancels native work before joining path
coordinators. Existing commands and checkpoint formats are unchanged.

## Register an existing run made by the older scripts

These are local code changes; copying them to the cluster does not update a
running Python process. Put the updated files in a separate copy of the existing
repository first. Include `graphvcfmerge.py`, `graphvcfmerge_kmer.py`,
`graphvcfmerge_nested.py`, `graphvcfmerge_scheduler.py`, `graphvcfmerge_checkpoints.py` and
`graphvcfmerge_resume.py`, plus `recover_graphvcfmerge.py` and
`diagnose_graphvcfmerge.py`. Also update `cohort_call_snakemake/`'s launcher,
`Snakefile`, and `workflow/scripts/graphvcfmerge_stages.py`. Keep their existing
dependencies and the compatible Linux KmerMatch executable.

From that updated repository copy, with the original Python environment:

```bash
MERGE_ROOT=/project2/mchaisso_100/walfred/projects/newphase/newBuild/annotation/test2/cohort_calls/tmp/graphvcfmerge

# Adopt already finished chromosome output pairs without rewriting them.
python scripts/graphvcfmerge_resume.py mark-done --shards-dir "$MERGE_ROOT" --all-existing

# Preserve chr1's completed top-level groups BEFORE stopping the old job.
python scripts/graphvcfmerge_resume.py mark-refined --shards-dir "$MERGE_ROOT" --chrom NC_060925.1

# Check what the next scheduler run will do.
python scripts/graphvcfmerge_resume.py status --shards-dir "$MERGE_ROOT"
```

The helper reads the original settings from
`cohort_calls/inputs/cohort_call.merge.json`. Use `--merge-config PATH` if the
configuration is elsewhere. If the original merge requested small-variant
output, also pass `--output-small` to each helper command.

`mark-done --all-existing` adopts only loci with both final `.part` and
`.promoted.tsv` files, plus `.small.part` when requested. An unfinished working
`.part` alone is insufficient. Existing checked markers are validated and kept.
For a specific chromosome, use `mark-done --chrom NAME` instead. This is an
explicit declaration that these older outputs are finished; the helper checks
file completeness and metadata, not biological correctness.

`mark-refined` discovers the saved `refined-groups.pkl` for chr1 and verifies that
all eligible top-level clusters were saved. If several spools match, supply
`--spool /absolute/path/to/refined-groups.pkl`. It copies the spool into the
durable checkpoint directory without modifying the original. An incomplete
spool is rejected. This registration can read a stable spool while the old job
is already in nested calling. Preserve it before stopping that job because old
temporary-directory cleanup may remove the source spool.

After registration succeeds, stop the old chromosome job and old launcher
before starting another run against the same output directory. Deploy the
updated scripts to the original repository path, then rerun your original
launcher command. For this migrated run, expect:

```text
[merge:status] NC_060925.1: refined-ready
[merge:resume] NC_060925.1: reusing completed groups; preparing nested calling
[merge:checkpoint] NC_060925.1: nested-ready saved in ...
[merge:nested] NC_060925.1: calling internal variants on ...
```

This first migration rebuilds top-level rows and nested candidates from the
saved groups, without repeating top-level alignment. Later restarts from
`nested.ready` skip that preparation too. The retained top-level groups still
come from the old run's algorithm; registration does not recompute them using
new grouping rules. Other unfinished chromosomes, if any, must also finish
before final concatenation.

## Manually adding a completion marker

`status` prints each chromosome's exact `merge.done` path, so you do not need to
calculate its hashed folder name. You may create an **empty** `merge.done` there
for a chromosome you know finished. The next run validates the output pair and
replaces that empty file with a marker recording settings and file stamps.
The `mark-done` command above performs those checks immediately.

Do not manually touch `refined.ready` or `nested.ready`: they require the saved
payload and associated data files. Use `mark-refined` for older runs. Keep the
checkpoint directory until recovery is no longer needed; it uses additional
disk space and is retained after success. Keep `merge.done` and final parts for
future scheduler checks. New workers hold a per-chromosome lock to prevent two
new-version jobs from writing the same chromosome simultaneously; older
workers do not honor that lock.

# Chromosome 1 merge resources

Slurm chromosome merging now gives chr1 a dedicated job with **twice the base
CPU and memory allocation** by default. The merger's `--processes` is set to
that job's allocated CPU count, so the additional CPUs are available to its
worker pools.

Nested variant calling uses one shared alignment pool within each chromosome
job. With 32 workers, a single remaining path can submit independent alignments
to all 32 slots, while multiple paths share the same 32-slot limit.

For example, with `--merge-processes 16` and the default memory policy:

| Merge job | CPUs / merger workers | Memory |
| --- | --- | --- |
| chr1 | 32 | 64G |
| Other chromosomes | 16 | 32G |

Control this using the cohort launcher's option:

```text
--merge-chr1-resource-multiplier 2
```

Use `1` to disable the extra allocation and restore the original round-robin
grouping. Larger positive integers are also supported. The setting is saved in
`inputs/cohort_call.run.json`. For an existing run invoked directly through
Snakemake, add `merge_chr1_resource_multiplier=2` to its `--config` values.
The stage helper exposes the corresponding `--chr1-resource-multiplier` option
under `run-chroms`.

Detection matches primary names `chr1`, `1`, and the versioned RefSeq accessions
`NC_060925` (CHM13) and `NC_000001` (GRCh38). It does not match chr10, chr11,
alternative contigs, or derived insertion-locus names containing these strings.

Explicit CPU and memory values in `--slurm-args` / `--sbatch-arg` are treated as
the base allocation before multiplying chr1's resources. With `--mem-per-cpu`,
the per-CPU amount stays the same; doubling CPUs doubles total memory.
Slurm's special `--mem=0` continues to request all node memory.

The job limit still bounds concurrent Slurm jobs. With a limit of one, the
dedicated chr1 job and the ordinary chromosome group run sequentially.
`tmp/graphvcfmerge/chrom_jobs/assignments.tsv` records each job's chromosome,
CPU count, and memory request (under the configured merge temporary directory).

Update `run_cohort_call_pipeline.py`, `Snakefile`, and
`workflow/scripts/graphvcfmerge_stages.py` together on the cluster. This applies
to newly submitted merge jobs; existing allocations are not resized. The saved
scan shards can be reused when rerunning the chromosome merge stage. Local runs,
scan/concatenation allocations, and wall-time limits are unchanged. More CPUs
help parallel work; they do not guarantee a twofold speedup for serial alignment.
## Original alternative calling intervals

Graph summarization writes `graph/summary/alternative_intervals.bed` from original
template intervals. Cohort calling detects this file automatically; override it
with `--alternative PATH`. The graph-CIGAR converter carries the projected calling
regions to the VCF caller, which excludes whole boundary-crossing variants. The
BED does not change path names or coordinate origins.
