# LinGraph

LinGraph is a variants-guided **pangenome graph builder** and an
assembly-based structural variant (SV) caller. Use `individual` for independent VCF calling against an existing
graph cache, or `graph` for a cohort run that calls variants together while
building a pangenome graph. LinGraph records lossless variant sequences and assembly coordinates
information in grVCF and can export cohort calls as a pangenome graph.

## Contents

1. [Workflow overview](#workflow-overview)
2. [Installation and requirements](#installation-and-requirements)
3. [Choose a run mode and prepare inputs](#choose-a-run-mode-and-prepare-inputs)
4. [For individual mode: download and reconstruct Win50KGraph](#for-individual-mode-download-and-reconstruct-win50kgraph)
5. [Call one or many samples independently](#call-one-or-many-samples-independently)
6. [Build a new pangenome graph from a fresh cohort run](#build-a-new-pangenome-graph-from-a-fresh-cohort-run)
7. [grVCF format and standard VCF conversion](#grvcf-format-and-standard-vcf-conversion)
8. [Build a pangenome GFA from grVCF](#build-a-pangenome-gfa-from-grvcf)
9. [Benchmark trio consistency](#benchmark-trio-consistency)
10. [Questions and answers](#questions-and-answers)

## Workflow overview

![LinGraph workflow: prepare and align local graphs, pair query and reference sequences, call and merge grVCF variants, and export a pangenome GFA.](figures/lingraph_workflow.png)

LinGraph partitions assemblies into loci, builds and aligns local graphs, and
uses those alignments to call variants. Cohort calls retain sequence and graph
relationships in grVCF and can be exported as a pangenome GFA. Novel-locus
discovery is optional. For independent VCF calls, `individual` reuses a supplied
graph cache; a cohort `graph` run builds or resumes the graph from the cohort
assemblies. The downloadable summary below supplies precomputed graph
information for individual calling.

[Download the workflow figure (PDF)](figures/lingraph_workflow.pdf).

## Installation and requirements

Run the examples from the repository root, which contains `readme.md`. Replace
example FASTA paths with your own files.

You need Python **3.9 or later**, Conda, `make`, and a C++17-capable compiler.
The workflows use Unix command-line tools and support local execution or SLURM.

```bash
cd /path/to/LinGraph
python3 install.py --conda-env lingraph -j 8
conda activate lingraph
export PATH="$PWD/scripts:$PATH"
python3 install.py --check-only
```

The installer creates or updates the named environment, builds the C++ tools
under `scripts/src/`, and copies the executables into `scripts/`. It installs:

| Dependency | Requirement or purpose |
| --- | --- |
| Snakemake **6.15.1**, tabulate **0.8.10** | Pinned workflow versions |
| HTSlib **≥1.19**, samtools | Indexed sequence access |
| BLAST+ **≥2.14.0**, including WindowMasker | Alignment and assembly soft masking |
| minimap2, Winnowmap2, meryl | Mapping and k-mer processing |
| bcftools | VCF normalization, compression, and indexing |
| numpy, polars ≥1.0, mappy, parasail, edlib, ssw-py | Python analysis and alignment libraries |

**Graph-alignment helpers** are included in `scripts/`, including
`graphcigarlight.py` and the BLAST/Winnowmap wrappers. The installer builds
`KmerStrd` and the other native tools. After installation, check helper discovery:

```bash
PYTHONPATH="$PWD/scripts" python3 -c 'from local_graph_whole_cigar import resolve_extension_dir; print(resolve_extension_dir())'
```

To use a separate helper installation, set `MINSETREF_GRAPH_TOOLS` to its directory.

`install.py --check-only` checks runtime dependencies without installing or
compiling. The helper check above verifies that the graph-alignment tools can
be located. Use `--no-install-deps` to build with an environment you manage yourself.

Graph mode requires at least **4 CPUs**; examples below use `-t 16`. For a full
cohort run, we recommend **32 CPUs and at least 64 GB RAM**. Larger cohorts or
high repeat complexity may need more memory and temporary storage. For cohorts
with thousands of haplotypes, use shared storage and the SLURM options described
in [the full command guide](scripts/LinGraph.md#slurm-and-dependencies).

## Choose a run mode and prepare inputs

Choose the mode based on the result you need:

| Mode | Use it for | Data needed and reuse |
| --- | --- | --- |
| `individual` | **Independent VCF calling**; one call per input haplotype, optionally merged afterward | A reconstructed graph cache and reference FASTA. Reuse a matching reference alignment cache when available. |
| `graph` (cohort mode) | **Call a cohort together while building or resuming its pangenome graph** | Cohort assembly list and reference. A fresh run builds the graph under `-G`; a completed graph can be resumed. |

`individual` can process multiple samples in one command, but calls each
haplotype independently and writes a separate VCF. The `graph` command is the
cohort workflow: it calls the cohort as part of building the graph.

Recommended windows: for individual calls, reuse a graph cache built with
balanced blocks (the supplied Win50KGraph cache uses these). For cohort graph
construction, use gene blocks; balanced blocks are an alternative when
prioritizing SV calling. The supplied gene and balanced BEDs use CHM13
coordinates, so use the matching CHM13 reference.

**Automatic blocks** are also available for cohort graph construction. Omit
`-b` and `--bed-grouped` and LinGraph chooses windows from the construction
inputs. Individual mode does not choose windows; it uses the ones in its graph
cache. Add `--make-graph` to a cohort run when you also want the pangenome GFA
file `cohort.gfa`.

### Prepare assembly FASTAs

Use one assembly per haplotype. A diploid sample normally has two entries,
such as `HG002_h1` and `HG002_h2`. The assembly-list label and the FASTA
contig prefix use different formats: the list uses an underscore before `h`,
while FASTA headers use `#` separators. For example, `HG002_h1` in the list
must have contig IDs beginning with `HG002#1#` (such as `HG002#1#chr1`). This
format change is required before graph construction or calling.

For example, preparation changes a raw contig header like `>chr1` to
`>HG002#1#chr1` for the `HG002_h1` assembly. The sample and haplotype prefix
identifies which assembly owns each contig; keep each contig under exactly one
haplotype prefix.

An input list has exactly two whitespace-separated columns, with no header:

```text
HG002_h1 /data/raw/HG002.h1.raw.fasta
HG002_h2 /data/raw/HG002.h2.raw.fasta
```

Save this as `raw_assemblies.list`, then prepare the assemblies:

```bash
python3 tools/prepare_assemblies.py \
  -q raw_assemblies.list -O prepared_assemblies \
  --contignamefix -j 4
```

Preparation soft-masks unmasked inputs, fixes contig prefixes, and writes
uncompressed FASTAs with adjacent `.fai` indexes. Use the generated
`prepared_assemblies/query_paths.prepared.txt` as the calling input list.
Add `--remask` when you want to regenerate an existing mask.

Input rules:

- The first column uses `SAMPLE_hN`, for example `HG002_h1`; each FASTA
  identifier must be `SAMPLE#N#CONTIG`, for example `HG002#1#chr1`. If the
  source FASTA has plain names such as `>chr1`, run the preparation command
  above with `--contignamefix`; it writes a prepared copy with the required
  prefix and does not change the source FASTA.
- For name fixing alone, this short command accepts either `HG002_h1` or
  `HG002#1` and writes `assembly.namefixed.fa` by default:

  ```bash
  python3 tools/namecontigsfix.py -i assembly.fa -n HG002_h1
  samtools faidx assembly.namefixed.fa
  ```

  Use `assembly.namefixed.fa` in the assembly list. Add `-o output.fa` to
  choose a different output path.
- Calling requires uncompressed, soft-masked FASTAs and matching `FASTA.fai`
  indexes. Preparation accepts compressed inputs.
- Relative FASTA paths are resolved from the input list's directory. Paths
  cannot contain whitespace; do not add an index column.
- **Every first-column haplotype name must be unique across the entire input
  list.** Do not assign the same name to different FASTAs, for example, two
  separate `HG002_h1` rows. Name collisions make inputs and per-haplotype
  outputs ambiguous and can cause one assembly's results to be mistaken for
  another's.
- Do not put one haplotype's prefix inside another's contig name, such as
  `HG002#1#HG002#2#chr1`. The name fixer rejects this collision instead of
  adding another prefix; correct the source FASTA header or assembly-list
  label so the contig belongs to exactly one haplotype.
- **`CHM13_h1` and `HG38_h1` can keep their native contig names and masking.**
  They still need matching indexes, created with `samtools faidx reference.fa`.
- Calling checks the first sequence and index entry. Run preparation explicitly
  before calling; otherwise calling stops, unless `--force-prepare` is given
  (see [the input format question](#why-is-a-specific-input-assembly-format-required-and-what-is-it)).

See [assembly preparation](tools/assembly_preparation.md) for more options.

## For individual mode: download and reconstruct Win50KGraph

**`Win50KGraph.tar.gz`** is the precomputed LinGraph graph-summary package for
the **balanced-block option**, recommended for individual sample SV calling.
It reconstructs all local graphs used in variant calling, including repetitive
regions. Download and extract it, reconstruct the graphs once, then reuse them
with `individual` mode for new haplotype assemblies.

**Figshare DOI:** [10.6084/m9.figshare.33930604](https://doi.org/10.6084/m9.figshare.33930604).
Download **`Win50KGraph.tar.gz`** (about 2.43 GB) from the Figshare record.
The download is available after the record is published and the DOI is active.

To download it from a terminal, run this command from the directory where you
want the archive saved. `-C -` resumes an interrupted download:

```bash
curl -L -C - https://ndownloader.figshare.com/files/69451449 -o Win50KGraph.tar.gz
```

Allow additional disk space for extraction and the reconstructed graphs, and
RAM for the uncompressed reference and packed sequences during reconstruction.

### 1. Extract the summary

From the directory containing the downloaded archive:

```bash
tar -xzf Win50KGraph.tar.gz
```

The package extracts into `cohort_minsetref_v3/summary/`. Its main files are:

| File or directory | Purpose |
| --- | --- |
| `local_graphs.tsv` | Local graph definitions, source intervals, and path metadata |
| `alternatives.fasta` and `.fai` | Packed alternative and novel sequences, with their index |
| `references/CHM13_h1_rig/` | Supplied CHM13 reference alignment and block coordinates |
| `references/HG38_h1_rig/` | Supplied GRCh38 reference alignment and block coordinates |

The archive contains the graph summary, not the reconstructed full graph. Run
LinGraph's `scripts/reconstruct_local_graph_folders.py` in Step 2 to build all
local graphs before using `Win50KGraph` for variant calling.

### 2. Reconstruct all local graphs

Use the matching **CHM13 FASTA used to build the package**, with its adjacent
`.fai` index. Reconstruction reads the summary and reference; the original
cohort assemblies are not required.

```bash
python3 scripts/reconstruct_local_graph_folders.py \
  -i cohort_minsetref_v3/summary \
  -r /data/references/chm13.fa --reference-haplotype CHM13_h1 \
  -o Win50KGraph -j 16 --resume
```

This writes the local graph FASTAs, template FASTAs, BED files, BLAST databases,
and graph lists under `Win50KGraph/`. It also links `summary/` and `references/`
to the extracted package, so **keep the extracted summary in place**. Repeat
the command with `--resume` to reuse completed reconstruction outputs.

### 3. Call variants in a new sample

After preparing your query assemblies, call against CHM13:

```bash
python3 scripts/LinGraph.py individual \
  -I prepared_assemblies/query_paths.prepared.txt \
  -G Win50KGraph \
  -r /data/references/chm13.fa --reference-name CHM13_h1 \
  --reference-caches Win50KGraph/references/CHM13_h1_rig \
  -O sample_calls_chm13 -t 16
```

The explicit cache option reuses the supplied reference alignment directly;
use it with the matching reference assembly. For GRCh38 calling, use its FASTA,
`--reference-name HG38_h1`, and
`--reference-caches Win50KGraph/references/HG38_h1_rig`, with a separate output
directory. Reconstruction still uses the package's CHM13 source reference.
For another calling reference, omit `--reference-caches` so LinGraph prepares
that reference's alignments.

## Call one or many samples independently

Use `individual` to call and annotate **one sample or a batch of samples
independently** against an **existing LinGraph graph**.
Pass the graph root or its `summary/` directory with `-G`. To obtain a graph,
[download and reconstruct Win50KGraph](#for-individual-mode-download-and-reconstruct-win50kgraph)
or build one with [cohort mode](#build-a-new-pangenome-graph-from-a-fresh-cohort-run).

The examples below use `cohort_graph` as that existing graph. Reference names
such as `CHM13_h1` must occur in its saved cohort list; otherwise, supply the
reference FASTA explicitly.

### Call multiple samples in one command

List the prepared haplotype assemblies in `samples.list`. For example, two
diploid samples have four entries:

```text
HG002_h1 /data/prepared/HG002.h1.prepared.fasta
HG002_h2 /data/prepared/HG002.h2.prepared.fasta
HG003_h1 /data/prepared/HG003.h1.prepared.fasta
HG003_h2 /data/prepared/HG003.h2.prepared.fasta
```

```bash
python3 scripts/LinGraph.py individual \
  -I samples.list -G cohort_graph -r CHM13_h1 \
  -O sample_calls_chm13 -t 16
```

This independently calls all four haplotypes and writes a VCF and coverage
reports for each under `sample_calls_chm13/samples/NAME/`. The existing graph
and reference alignments are reused across inputs. By default, this command
produces separate calls; add `--merge` to also merge them after calling.

### Use CHM13, GRCh38, or another assembly as the reference

Call against CHM13:

```bash
python3 scripts/LinGraph.py individual \
  -I prepared_assemblies/query_paths.prepared.txt \
  -G cohort_graph -r CHM13_h1 -O sample_calls_chm13 -t 16
```

Call the same assemblies against GRCh38:

```bash
python3 scripts/LinGraph.py individual \
  -I prepared_assemblies/query_paths.prepared.txt \
  -G cohort_graph -r HG38_h1 -O sample_calls_hg38 -t 16
```

Use an external reference FASTA, including a reference absent from the saved
cohort:

```bash
python3 scripts/LinGraph.py individual \
  -I prepared_assemblies/query_paths.prepared.txt \
  -G cohort_graph -r /data/references/custom.fa \
  --reference-name CUSTOM_h1 -O sample_calls_custom -t 16
```

The custom reference must follow the preparation rules for `CUSTOM_h1`. For
native CHM13 or GRCh38 FASTAs, use `--reference-name CHM13_h1` or
`--reference-name HG38_h1`. You can also select any other haplotype in the saved
cohort by name, for example `-r HG002_h1`.

Omitting `-r` uses the first saved cohort assembly. Keep a separate output
directory for each reference. Changing the reference requires a new calling
run. For `HG38_h1`, the pipeline uses chr1–22, chrX, and chrY; chrM and
non-primary scaffolds are excluded.

To supply one FASTA directly, replace `-I ...` with
`-i /data/prepared/HG002.h1.prepared.fasta --sample HG002_h1`.

### Add another reference: HG19 example

A graph can call against any reference once that reference's cache exists in
`GRAPH/references/NAME_rig`. Win50KGraph ships caches for `CHM13_h1` and
`HG38_h1`; this example builds one for HG19 (UCSC, chr1–22, X and Y).

Download HG19, keep the main chromosomes and prepare it. A reference other than
CHM13 or GRCh38 uses the `NAME#1#` contig prefix, so `--contignamefix` renames
`chr1` to `HG19#1#chr1`:

```bash
mkdir -p logs references && cd references
curl -L -O https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gzip -dc hg19.fa.gz > hg19.fa
samtools faidx hg19.fa
samtools faidx hg19.fa chr{1..22} chrX chrY > hg19_main.fa
python3 ../LinGraph/tools/prepare_assemblies.py \
  -i hg19_main.fa --name HG19_h1 -O hg19 \
  --contignamefix --remask --threads 16 > ../logs/prepare_hg19.log 2>&1
cd ..
```

Build the cache once. `--reference-only` aligns the reference to every local
graph (hours, like calling one haplotype), writes
`Win50KGraph/references/HG19_h1_rig` (or `HG19_h1_<fingerprint>_rig` if a
different HG19 already uses that name), and calls no samples. Rerun the same
command to resume:

```bash
python3 LinGraph/scripts/LinGraph.py individual \
  -G Win50KGraph -r references/hg19/hg19_main.fa \
  --reference-only -O hg19_cache -t 64 > logs/hg19_cache.log 2>&1
```

The name `HG19_h1` is read from the contig prefix; `--reference-name HG19_h1`
states it explicitly. To keep the cache elsewhere, add `--reference-caches DIR`:
the cache is built in `DIR`, and later calls use it with the same
`--reference-caches DIR`. Without this step, the first HG19 call builds the
cache in its own output folder. Then call against HG19 like any other reference; the
cache is found and reused automatically:

```bash
python3 LinGraph/scripts/LinGraph.py individual \
  -I prepared_assemblies/query_paths.prepared.txt \
  -G Win50KGraph -r references/hg19/hg19_main.fa \
  -O sample_calls_hg19 -t 64 > logs/calls_hg19.log 2>&1
```

VCF records use the prepared names, e.g. `HG19#1#chr1`. To restore plain
names, see [Reference contig names](#reference-contig-names).

### Reuse reference alignments

Calling needs a reference cache: the reference aligned to every local graph.
Caches live in folders such as `cohort_graph/references/CHM13_h1_rig/`, each
with a `reference_cache.json` recording the reference name, the alignment mode
and a **contig fingerprint**: the SHA-256 of sorted chromosome-name/length
pairs from the `.fai`. The fingerprint ignores FASTA contents, file names,
masking, wrapping, `.fai` row order, offsets and timestamps.

`individual` finds its cache automatically:

1. It checks `GRAPH/references/NAME_rig`, then
   `GRAPH/references/NAME_<fingerprint>_rig` (the first 12 characters of the
   reference's own fingerprint), and uses the first that is a complete cache
   with the same reference name, mode and fingerprint. It prints the folder and
   the fingerprint it matched.
2. Without a match, it builds the cache in `OUTPUT/references/NAME_rig` and
   uses it. Rerunning with the same `-O` resumes or reuses it.

To build a cache once for later calls, use `--reference-only` (see the HG19
example above). It writes to `GRAPH/references/NAME_rig`; when that folder
already holds a different reference of the same name, it writes to
`GRAPH/references/NAME_<fingerprint>_rig` instead, which step 1 also finds.

To use a specific cache folder instead, add:

```bash
--reference-caches /path/to/reference-caches
```

For `-r CHM13_h1` (or `--reference-name CHM13_h1`), that directory supplies
`CHM13_h1.align.txt` and `CHM13_h1.align.txt_blocks.bed`. LinGraph uses these
files directly without reference preparation or cache, hash, mode, or content
validation. Their compatibility is the user's responsibility. The reference
FASTA and its `.fai` are still required for variant calling.

When reconstructing a graph package, `reconstruct_local_graph_folders.py` also
links `OUTPUT/references` to `INPUT/references` if that input directory exists.
It retains a matching existing link and refuses to overwrite a different destination.

### Outputs and useful options

Each input haplotype produces:

```text
sample_calls_chm13/samples/HG002_h1/
  HG002_h1.vcf
  HG002_h1.coverage.summary.tsv
  HG002_h1.coverage.missing_reference.bed
  HG002_h1.coverage.missing_query.bed
```

The coverage files report represented and uncovered assembly bases.

| Option | Effect |
| --- | --- |
| `--all` | Include SNPs, indels, and SVs of all sizes; multiple inputs also produce separate merged VCFs. |
| `--merge` | Merge calls when more than one haplotype is supplied; by default, writes `cohort.sv.vcf`. |
| `--svcutoff 50` | Set the SV size boundary to 50 bp; the default is 20 bp. |
| `-L partitions.list` | Call only the listed graph partitions; omit to use all partitions. |
| `--satellite-bed satellites.bed` | Add satellite annotations to coverage reports. |
| `--dry-run` | Check inputs and show planned commands without running or writing files. |

The default selects independent SV calling with one VCF per input haplotype.
With multiple inputs, `--merge` or an explicit variant-selection option
(`--all`, `--svonly`, `--svindel`, or `--snp`) also merges the resulting VCFs.
The calling step remains independent for each input. A single input produces
only its individual VCF. The cutoff applies to **merged representatives**;
individual `.vcf` files can retain smaller candidate events for merging.

A compact summary contains `local_graphs.tsv` and `alternatives.fasta`.
Reconstruct it with the package's source reference before selecting a different
calling reference. See the [Win50KGraph instructions](#for-individual-mode-download-and-reconstruct-win50kgraph).

## Build a new pangenome graph from a fresh cohort run

**Cohort mode builds or resumes a pangenome graph from a cohort run.** The
workflow calls variants across the assemblies as part of graph construction.
Use this mode for a fresh graph build; use `individual` when you already have a
graph cache and need individual VCFs. It supports thousands of haplotypes.

The `graph` command builds/resumes local graphs and runs cohort variant calling
as part of that workflow. Use `--exact` for the merged cohort output; it is the
default and uses the indexed assemblies to realign merged SVs. Add `--make-graph`
to export the resulting pangenome GFA in the same run.

Create `cohort.list` with the reference and prepared haplotypes:

```text
CHM13_h1 /data/references/chm13.fa
HG002_h1 /data/prepared/HG002.h1.prepared.fasta
HG002_h2 /data/prepared/HG002.h2.prepared.fasta
HG003_h1 /data/prepared/HG003.h1.prepared.fasta
HG003_h2 /data/prepared/HG003.h2.prepared.fasta
```

### Recommended alternatives for CHM13

[data/alternatives.fa](data/alternatives.fa) is the recommended alternative
sequence catalog to use alongside CHM13. Its adjacent `.fai` index is included.
If you edit or replace the FASTA, regenerate the index:

```bash
samtools faidx data/alternatives.fa
```

Import the catalog when building a new cohort graph against CHM13:

```bash
python3 scripts/LinGraph.py graph \
  -I cohort.list -G cohort_graph -O cohort_calls \
  -b windowprofs/geneblocks.bed --bed-grouped \
  -r CHM13_h1 --alternative data/alternatives.fa --exact -t 16
```

The catalog complements the CHM13 reference; it does not replace its FASTA or
its entry in `cohort.list`. Imported alternatives become part of the saved
graph inputs. Use the same catalog when resuming a build. Reuse the resulting
graph and its saved templates for later calling and GFA export.

### Choose or change the reference/backbone

Follow these steps to choose which cohort haplotype is used as the reference
for variant calling and graph output:

1. **Check `cohort.list`.** Every row starts with a unique haplotype name,
   followed by the path to its prepared assembly. Use the exact name from this
   column.
2. **Choose any haplotype in the cohort.** For example, to use `HG002_h2`, add
   `-r HG002_h2`. If you leave out `-r`, LinGraph uses the first haplotype in
   `cohort.list`.
3. **Build the graph and calls.** For a fresh run, pass the cohort list,
   graph-cache directory, output directory, and your selected backbone:

   ```bash
   python3 scripts/LinGraph.py graph \
     -I cohort.list -G cohort_graph -O cohort_calls_HG002_h2 \
     -b windowprofs/geneblocks.bed --bed-grouped \
     -r HG002_h2 --exact --make-graph -t 16
   ```

4. **Create results with a different backbone later.** After the first run has
   completed, reuse the same `-G` cache, omit `-I` and graph-construction
   options, change `-r`, and give `-O` a new directory. For example:

   ```bash
   python3 scripts/LinGraph.py graph \
     -G cohort_graph -O cohort_calls_HG003_h1 \
     -r HG003_h1 --exact --make-graph -t 16
   ```

   LinGraph reuses the saved local graph alignments and reruns the
   reference-dependent matching, lifting, VCF, merge, and GFA steps. Use a
   normal `graph` run for this switch; `--recall-only` uses the reference from
   the previous run and cannot change backbones.

The selected `-O` directory receives per-assembly calls under
`samples/NAME/` and the merged cohort SV, indel, and SNP grVCFs. Keep these
grVCFs for downstream graph export. The `--all`, `--svonly`, `--svindel`, and
`--snp` options are primarily for selecting VCF contents in individual runs;
they do not control graph construction.

Keep `-G` and `-O` **separate and non-nested**. Repeat the same command to resume
completed work and retry unfinished stages. A completed graph retains its
cohort list, so `-I` can then be omitted.

Useful additions:

- `--static-block`: retain initial templates and skip additional novel-locus
  discovery and refinement. Calling and merging still run.
- `--slurm --slurm-jobs 20 --slurm-args '--account=myaccount --partition=compute'`:
  submit workflow jobs through SLURM. Activate the environment before launching.

For existing results, use `--merge-only` to merge saved individual VCFs, or
`--recall-only` to regenerate VCFs and merge from saved calling alignments. For
example, to regenerate exact merged cohort outputs:

```bash
python3 scripts/LinGraph.py graph -G cohort_graph -O cohort_calls \
  -r CHM13_h1 --recall-only --exact -t 16
```

Recall requires the saved alignments and calling inputs. See
[cohort calling](scripts/cohort_call_snakemake/README.md) for restart and merge details.

## grVCF format and standard VCF conversion

See the [Graph Recursive VCF (grVCF) guide](docs/LinGraph_Graph_Recursive_VCF.pdf)
for the CIGAR and recursive representations, nested variants, and validation.

### What grVCF stores

grVCF preserves information needed to describe variation within repetitive and
non-reference sequence. It uses the familiar VCF columns (`CHROM`, `POS`, `ID`,
`REF`, `ALT`, `QUAL`, `FILTER`, `INFO`, `FORMAT`, and sample columns), with extra
graph-aware annotations. LinGraph currently writes these files with a `.vcf`
extension, including individual and merged outputs.

A merged row describes one **representative allele**: the largest insertion of
its group, but the smallest deletion. Row IDs are unique: `I_`/`D_`/`S_` for
insertions, deletions and SNPs. A nested row's CHROM is its parent's ID and its
ID is `<I|D>_<parent>_<pos>_<size>`; insertions of different bases at one
position and size get `a`, `b`, ... appended. SNPs are `S_<chrom>_<pos>_<alt>`.
Each haplotype's FORMAT fields retain its supporting observations and how they
relate to that allele.
This preserves variation both on the chosen reference and within alternative,
novel, or inserted sequence, including variants nested inside other alleles.

Two downstream uses are supported:

| Goal | Input and tool |
| --- | --- |
| Use reference-based VCF software | Convert grVCF with `tools/grvcf_to_vcf.py` |
| Build a pangenome graph | Pass the original grVCFs to `merged_vcf_to_gfa.py` or use `LinGraph.py graph --gfa-only` |

### Main fields

| Field | Meaning |
| --- | --- |
| `ALT` | SNP bases or symbolic alleles such as `<INS>`, `<DEL>`, and `<SUB>` (replacement). |
| `INFO/SVTYPE`, `END`, `SVLEN` | Event type, endpoint, and representative signed length. |
| `INFO/SEQ` | Representative inserted/replacement sequence, which may be graph encoded. |
| `INFO/EXTENDGRAPHCIGAR` | Representative alignment encoding, including graph targets when present. |
| `FORMAT/GT` | Haplotype call: `1` supports the variant, `0` is callable reference, `.` is no call. |
| `FORMAT/TYPE`, `SIZE` | Each observation's event type and signed size. |
| `FORMAT/EXTENDGRAPHCIGAR` | Each observation's alignment to the representative. |
| `FORMAT/BASE` | Observed alternate base for a SNP. |
| `FORMAT/ASSEMBLYCONTIG`, `QUERYCOORD` | Original assembly contig and zero-based query coordinate/range with strand. |
| `FORMAT/TEMPLATEOFFSET` | Signed displacement from the row's `POS` to an observation's breakpoint. |
| `FORMAT/ALLELENAME`, `LABEL_H` | Source allele and local graph provenance. |

The new FORMAT layout exposes these annotations as separate named fields.
Multiple observations are comma-separated in matching order across the
per-observation fields; they do not represent extra genotype copies. The older
packed `HSV` and `HSNP` layouts remain supported by the converter.

For example, the default SV FORMAT is:

```text
GT:TYPE:SIZE:EXTENDGRAPHCIGAR:ASSEMBLYCONTIG:QUERYCOORD:TEMPLATEOFFSET:ALLELENAME:LABEL_H
```

`TEMPLATEOFFSET=5` means that observation's breakpoint lies five bases after
the row's displayed `POS`. Thus a merged site retains the individual
breakpoints as well as its representative allele. Optional provenance fields
may be omitted by the calling options; always read the record's FORMAT column.

**Coordinates need care.** SNP positions are one-based. Symbolic variants on
reference or supplied template sequences describe the interval `[POS, END)`;
an insertion has `END=POS`, and a boundary can be zero. Nested variants can use
an insertion allele ID as `CHROM`; their default `POS` is the insertion offset
plus one. Graph CIGARs can reference other alleles or named FASTA sequences.
These details are resolved by LinGraph's converters.

### Reference contig names

LinGraph names reference contigs so that chromosomes of different references
never share a name: `chr1` of GRCh38, CHM13, and HG19 are different sequences,
and one graph or cohort can use several references. CHM13 keeps its NCBI RefSeq
accessions (`NC_060925.1` ... `NC_060948.1`), GRCh38 keeps `chr1` ... `chrY`,
and every other reference or assembly gets the `NAME#N#` prefix
(`HG19#1#chr1`, `HG002#1#chr1`). VCF records, `##contig` lines, and coverage
headers use these names.

If you prefer the original chromosome names, rename the finished VCFs with
`tools/vcf_chromfix.py`. It also renames the reference fields of LinGraph's
headers and leaves nested rows pointing at their parents:

```bash
# CHM13: NC_060925.1 ... NC_060948.1 -> chr1 ... chr22, chrX, chrY
python3 tools/vcf_chromfix.py -i cohort.sv.vcf -o cohort.sv.chr.vcf --CHM13fix
# Other references: HG19#1#chr1 -> chr1
python3 tools/vcf_chromfix.py -i cohort.sv.vcf -o cohort.sv.chr.vcf --noprefix
```

Rename the indel and SNP files the same way. Keep the original files as input
to LinGraph's merge, GFA export, and converters: they match the prepared FASTA
names. See the
[renaming guide](tools/README.md#rename-reference-chromosomes-in-a-vcf) for
`--fixtable` and other options.

### Known incompatibility: records at POS 0

Some records have `POS=0`: an insertion before the first base of its sequence,
mostly a nested variant at the very start of its parent allele (`CHROM` = the
parent row ID). grVCF uses the interval convention above, where a boundary can
be zero, but standard VCF positions start at 1 (`POS=0` is reserved for
telomeric breakends), so some standard tools do not handle these records. With
bcftools/htslib 1.22, for example, `bcftools view` and `bcftools sort` accept
them, but `tabix` warns `Coordinate <= 0`, region queries (`-r`) skip them, and
`bcftools norm` reports a REF mismatch because there is no base 0. LinGraph's
own tools read them as intended. `tools/grvcf_to_vcf.py` removes them, together
with all nested rows, when it exports a standard VCF (see
`position_zero_removed` below).

### Convert to standard VCF

`tools/grvcf_to_vcf.py` exports a plain VCF in one streaming pass, without
reading a reference FASTA. Convert one file at a time:

```bash
python3 tools/grvcf_to_vcf.py \
  -i cohort_calls/cohort.sv.vcf \
  -o cohort_calls/cohort.sv.standard.vcf
```

Repeat for the indel and SNP files. Plain and gzip-compressed inputs are
accepted. The converter:

- Removes nested records on `I_`, `D_`, `SUB_`, or `DUP_` parents (`INS_`/`DEL_`
  in older runs), and
  records at `POS=0`, which lack a VCF anchor base.
- Converts `<INS>` with plain `INFO/SEQ` into explicit bases (`ALT=REF+SEQ`).
  Other symbolic alleles, including insertions with graph-encoded sequence,
  remain symbolic; it does not decode them during this export.
- Preserves record order, IDs, REF, QUAL, FILTER, sample names, and GT values.
  INFO retains only `SVTYPE`, `END`, `SVLEN`, and `NSUP`; FORMAT retains GT.
- Does not validate against a reference, normalize or sort alleles, calculate
  allele frequencies, or produce an unplaced-record sidecar. Non-nested
  alternative-locus coordinates remain in their original coordinate system.

This is a representative-allele export, **not a lossless assembly export**.
Keep the original grVCF and its dependencies for nested variation and graph
construction. Check the reported `written`, `nested_removed`, and
`position_zero_removed` counts before using the output downstream.

Use `.vcf.gz` output for BGZF compression (requires `bgzip`); the converter
does not create an index. `--force` replaces an existing output. The legacy
`-r` option is accepted but ignored; `-a` and `--normalize` are not supported.
See [the converter guide](tools/README.md#convert-grvcf-to-standard-vcf) or
`python3 tools/grvcf_to_vcf.py --help` for supported options.

## Build a pangenome GFA from grVCF

Graph export uses the **original merged grVCFs**, their graph annotations, and
the matching assembly sequences to create a pangenome `.gfa` file. Keep these
grVCFs when converting copies to standard VCF: standard conversion removes
information needed for graph reconstruction.

### Export through LinGraph

To export the merged cohort VCFs produced above without calling or merging again:

```bash
python3 scripts/LinGraph.py graph -G cohort_graph -O cohort_calls \
  -r CHM13_h1 --gfa-only -t 16
```

The output is **`cohort_calls/cohort.gfa`**. Use the same graph, reference, and
merged cohort VCFs as the calling run. Keep its original assembly FASTAs,
indexes, and saved template files available: graph export needs their sequences.

To build with gene blocks, call the cohort, and export the GFA in one run, add
`--make-graph`:

```bash
python3 scripts/LinGraph.py graph \
  -I cohort.list -G cohort_graph -O cohort_calls \
  -b windowprofs/geneblocks.bed --bed-grouped \
  -r CHM13_h1 --alternative data/alternatives.fa --exact --make-graph -t 16
```

The default export is **rGFA**, a GFA graph with stable sequence coordinates.
It includes the reference backbone and branches for the selected SNPs, indels,
and SVs. Its named paths describe source sequences and variant alleles.

- Cohort graph runs use `--exact` for merged call inputs. Use
  `--insertion-only 50` to export only insertions of at least 50 bp.
- `--max-node-length` sets the maximum segment length; the default is 1,024 bases.
- An explicit `-r` selects that haplotype as the backbone. Without `-r`, export
  uses `GRAPH/inputs/reference_alternatives_novels.fa`, the combined catalog of
  reference, alternative, and novel sequences.

### Convert grVCF files directly

To select the input files explicitly, use `merged_vcf_to_gfa.py`. An exact
cohort run produces SV, indel, and SNP files; pass all three so nested variants
can resolve their parent alleles across files. Here `--all` tells this
standalone converter which input categories were supplied; cohort runs use
`--exact` as shown above.

```bash
python3 scripts/merged_vcf_to_gfa.py \
  -v cohort_calls/cohort.sv.vcf \
     cohort_calls/cohort.indel.vcf \
     cohort_calls/cohort.snp.vcf \
  -q cohort.list \
  -r /data/references/chm13.fa --reference-haplotype CHM13_h1 \
  --graph-folder cohort_graph \
  --local-reference-templates cohort_calls/checkpoints/local_reference_templates.fa \
  --local-path-fasta cohort_graph/summary/alternatives.fasta \
  --all --processes 16 -o cohort_calls/cohort.gfa
```

`cohort.list` must map VCF sample names to their original prepared FASTAs.
The reference and template files must match the calling run. All sequence
FASTA inputs need adjacent `.fai` indexes. The local path catalog resolves graph
coordinates; the lifted template file supplies placements relative to the
selected reference. Use `--local-path-fasta` for the summary catalog, since
`-a` instead adds entire FASTA records as extra graph paths.

The converter writes GFA segments (`S`), links (`L`), and named paths (`P`),
with rGFA tags recording sequence origins and offsets. SNPs and inserted or
replacement sequences form branches; deletions connect across skipped sequence.
BED and mapping sidecars accompany the GFA for tracing coordinates. Use
`--validate-only` to check sequence and topology without writing a GFA.

See [the GFA export guide](scripts/merged_vcf_to_gfa.md) for coordinate conventions,
additional filters, and graph validation.

## Benchmark trio consistency

The [benchmark directory](benchmark/README.md) provides two tools for checking
whether child variants have a match in either parent:

- **[QuickTriocheck.py](benchmark/QuickTriocheck.README.md):** compare separate
  haplotype VCFs or one multi-sample VCF using position, variant type, and size.
  Requires only Python's standard library.
- **[truvari_trio.sh](benchmark/truvari_trio.README.md):** run eight pairwise
  Truvari benchmarks for two child and four parental haplotypes.

Each guide explains requirements, coverage filtering, commands, and outputs.
Use calls made against the same reference and retain their coverage headers.

## Questions and answers

### Why is a specific input assembly format required, and what is it?

Each haplotype is one uncompressed FASTA with an adjacent `.fai` index, and
every contig is named `SAMPLE#N#CONTIG` (for example `HG002#1#chr1` for
`HG002_h1`); see [Prepare assembly FASTAs](#prepare-assembly-fastas). The
prefix keeps every contig name unique across all assemblies of a run, so each
call can always be traced back to its contig in the original assembly; plain
names such as `chr1` would collide between assemblies. Assemblies must also be
soft-masked, which helps the aligners.

If inputs are not in this format, LinGraph stops and prints a preparation
command for each one. The lazy mode, `--force-prepare`, instead prepares
formatted copies in `OUTPUT/PreparedAssemblies/` and deletes them after the run
finishes. It needs extra temporary disk space for a full copy of every
assembly it prepares. Use it for a few samples; for a large cohort, we highly
recommend preparing the assemblies in the correct format beforehand with
`tools/prepare_assemblies.py`.

### What are local graphs, and how are they made?

Local graphs apply divide and conquer to make graph building and SV calling
faster and simpler. The genome is divided into blocks; each block's sequences
from all assemblies form one small local graph, which is built, aligned and
called on its own, and the results are joined back on the reference.

Two blocking profiles are provided in `windowprofs/` (both CHM13 coordinates,
used with `-b BED --bed-grouped`):

1. **Balanced blocks** (`balanceblocks.bed`): chosen by automatic code to
   optimize the local alignments and avoid breakpoints in repeats. This is
   better for SV calling (about 0.1-0.2% fewer errors, a 10-20% improvement).
   The supplied Win50KGraph cache uses these blocks.
2. **Gene blocks** (`geneblocks.bed`): adjusted to gene coordinates, which
   helps downstream gene-based annotation, for example by genotyping tools
   such as ctyper. **Highly recommended when building a graph.** SV calling
   is slightly less accurate because block breakpoints can fall in repeats,
   but the additional errors are generally in tandem repeats, not in gene
   bodies.

### How long does a run take?

About **1 hour per haplotype on 64 cores** for individual VCF calling. Merging
depends on the cohort size but usually adds 10-20% of the individual calling
time. With `--slurm`, haplotypes are called in parallel (up to `--slurm-jobs`,
default 20, at once).

### How much memory is used, and how do I increase it after an out-of-memory failure?

Locally, LinGraph uses the memory of the machine. With `--slurm`:

| Job | Default memory | Increase with |
| --- | --- | --- |
| `individual`: one job per haplotype (and the `--merge` job) | 64G | `--slurm-memory` |
| `graph`: per-sample matching, lifting, graph CIGAR, gap filling, VCF | 32G each | `--match-memory`, `--genomelift-memory`, `--graphcigar-memory`, `--gapfill-memory`, `--vcf-memory` |
| `graph`: cohort extraction and local templates | 64G | `--extraction-memory` |
| `graph`: run-once merge scans, concats, publish, GFA export | 64G below 100 sample VCFs, 128G from 100 on | `--slurm-memory` |
| `graph`: per-chromosome merge jobs | 2G per CPU below 100 samples, 64G from 100 on (chr1 doubled) | `--merge-memory` |
| `graph`: novel-locus discovery | 128G | `--novel-slurm-memory` |

To find the failed job, check `sacct -j JOBID` (state `OUT_OF_MEMORY`) and the
logs under `OUTPUT/lingraph/`. Then repeat the same command with the larger
memory option: finished stages are reused and only the failed ones run again.
See `--help-all` for every per-stage option.

### How do I output versions of the graph with different references as the backbone?

Run `graph` again with the same `-G`, a different `-r` (any haplotype in the
cohort list) and a new `-O`, with `--make-graph`. Saved local graph alignments
are reused; only the reference-dependent steps run again. See
[Choose or change the reference/backbone](#choose-or-change-the-referencebackbone).

### How are non-reference sequences represented in the VCF?

At three levels:

1. **Alternatives**: novel genes found in open-source assemblies, supplied
   with `--alternative` (for example `data/alternatives.fa`).
2. **Novel loci**: large insertions found in the cohort, usually novel
   immune-related genes.
3. **Insertions**: called dynamically in each sample. Variants inside an
   inserted sequence are nested rows whose `CHROM` is the insertion's row ID.

### Where do I find the alternative and novel sequences?

Which sequences are alternative or novel depends on the linear reference: a
sequence absent from CHM13 can be present in GRCh38. So each calling run
writes the sequences for its own `-r` reference:

| Run | Files |
| --- | --- |
| `graph` | `OUTPUT/checkpoints/local_reference_templates.fa` and `OUTPUT/checkpoints/alternative_loci.fa` |
| `individual` | `OUTPUT/samples/NAME/local_reference_templates.fa` and `OUTPUT/samples/NAME/alternative_loci.fa` |

These are the sequences that VCF rows not placed on the reference are reported
on. Each record header also records where it lifts onto the reference
(`reference=CONTIG:START-END:STRAND`), or `lift_status=unmapped` when it does
not. Runs from before `alternative_loci.fa` was added have only
`local_reference_templates.fa`.

The graph itself keeps every local graph path that is not an exact slice of
its construction reference in `GRAPH/summary/alternatives.fasta`. When a graph
is built, `GRAPH/inputs/reference_alternatives.fa` holds the reference with
the supplied alternatives, `GRAPH/novel_loci.fa` the discovered novel loci
(with `--find-novel-loci`), and `GRAPH/inputs/reference_alternatives_novels.fa`
all of them together. Supplied alternatives are tagged
`sequence_role=imported_alternative` in their headers.

### How are duplicated genes represented?

By a primary alignment, reported as an insertion, and an alternative alignment
to the closest reference gene.

### How do I convert grVCF to a normal VCF?

Use `tools/grvcf_to_vcf.py`; see [Convert to standard VCF](#convert-to-standard-vcf).

### How do I add my own samples to known VCFs?

1. Call your samples with `individual`, using the same graph and reference as
   the known VCFs.
2. Split the known merged grVCFs (exact or CIGAR version) into one grVCF per
   sample; each sample's own call is rebuilt from the CIGAR in its sample
   fields:

   ```bash
   python3 tools/convert_merged_grvcf.py split \
     -v known/cohort.sv.vcf known/cohort.indel.vcf known/cohort.snp.vcf \
     -r reference.fa -o known_individual
   ```

3. Merge the split grVCFs together with your samples' grVCFs into new merged
   grVCFs:

   ```bash
   python3 tools/merge_grvcfs.py -I all_vcfs.list -O merged_new -t 16
   ```

See [the tools guide](tools/README.md#split-merged-grvcfs-or-convert-between-exact-and-cigar)
for `--exact` merging and other options.

### Which regions are not included?

1. Regions with N gaps.
2. The edges of scaffolds.
3. Constitutive heterochromatin and other simple repeats longer than 1 Mb.

Check each haplotype's coverage report, `samples/NAME/NAME.coverage.summary.tsv`,
with the uncovered intervals in `NAME.coverage.missing_reference.bed` and
`NAME.coverage.missing_query.bed`.

## Repository layout

```text
readme.md
install.py               Root installer entry point
scripts/                 Main pipeline scripts and installer implementation
  graph_build_snakemake/ Graph-building workflow
  cohort_call_snakemake/ Cohort-calling workflow
  src/                   C++ sources
figures/                 Workflow figure (PNG and PDF)
data/                    Recommended CHM13 alternatives and FASTA index
docs/                    Format and usage guides
tools/                   Preparation and VCF conversion utilities
windowprofs/             Gene and balanced block-interval BED files
benchmark/               Trio benchmarks and a README for each script
tests/                   Regression tests
```

Run the trio-reporting regression tests with:

```bash
python3 -m unittest discover -s tests -v
```

Compile the native tools during installation; platform-specific binaries are
not included.

## More help

```bash
python3 scripts/LinGraph.py --help
python3 scripts/LinGraph.py individual --help-all
python3 scripts/LinGraph.py graph --help-all
```

For advanced options, compact graph packages, and SLURM execution, see
[the full LinGraph guide](scripts/LinGraph.md).
