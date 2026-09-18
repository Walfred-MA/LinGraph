# Assembly preparation

`prepare_assemblies.py` replaces the PATs masking Snakemake workflow with one
Python command. It processes the standard two-column assembly list:

```text
HG002_h1 /assemblies/HG002.h1.fa
HG002_h2 /assemblies/HG002.h2.fa
```

Run it with:

```bash
python3 tools/prepare_assemblies.py \
  -q query_paths.txt \
  -O prepared_assemblies \
  -j 4 \
  --contignamefix
```

To prepare one FASTA without creating an input assembly list, use `-i` and
provide its required sample/haplotype name:

```bash
python3 tools/prepare_assemblies.py \
  -i assembly.fa \
  -O prepared_assembly \
  --name 'HG002#1' \
  --contignamefix
```

`--name HG002_h1` is equivalent. Direct-input mode still creates
`prepared_assembly/query_paths.prepared.txt`, so its result can be passed
directly to the graph pipeline.

The command writes uncompressed prepared FASTAs, their `.fai` indexes, and
`prepared_assemblies/query_paths.prepared.txt`. Use that rewritten list as the
graph pipeline's `-q` input.

The tool writes separate copies and refuses to overwrite an input FASTA.
Use the generated list with `LinGraph.py graph -I` or `LinGraph.py singular -I`.
The calling pipelines only check input format; they never mask, rename, or
index the supplied assemblies automatically. `LinGraph.py prepare` remains
an explicit shortcut to this standalone tool.

Prepare native `CHM13_h1` and `HG38_h1` references separately without
`--contignamefix` to preserve names such as `NC_060925.1` and `chr1`.

For name fixing alone, without masking, use a separate output file:

```bash
python3 tools/namecontigsfix.py -i assembly.fa -n 'HG002#1' -o assembly.namefixed.fa
samtools faidx assembly.namefixed.fa
```

## Behavior

- A FASTA is considered soft-masked when at least one lowercase `a`, `c`, `g`,
  or `t` occurs in its sequence. An already masked FASTA skips WindowMasker.
- `--remask` runs WindowMasker regardless of the detected input state.
- `--contignamefix` checks every identifier against the list name. For
  `HG002_h1`, valid identifiers are `HG002#1` or start with `HG002#1#`. A
  prefix is added only when necessary, and duplicate prefixes are collapsed.
- Name repair aborts if it would create duplicate FASTA identifiers.
- `samtools faidx` is always run on every resulting FASTA.
- Inputs ending in gzip data are accepted and written as ordinary,
  faidx-compatible FASTA files.

WindowMasker is the executable used by the linked PATs workflow (the tool name
is `windowmasker`, despite the common `winnowmask` spelling). Install the two
runtime dependencies, for example, with:

```bash
conda install -c bioconda -c conda-forge blast samtools
```

To force a fresh mask:

```bash
python3 tools/prepare_assemblies.py \
  -q query_paths.txt \
  -O prepared_assemblies \
  -j 4 \
  --remask \
  --contignamefix
```

`-j` controls how many assemblies are processed concurrently. Each
WindowMasker process is single-threaded, so choose this based on available
memory and I/O bandwidth.

For a list containing exactly one assembly, `--name` can override the name in
the list. Both PATs list notation and FASTA-prefix notation are accepted:

```bash
python3 tools/prepare_assemblies.py \
  -q one_assembly.txt \
  -O prepared_assembly \
  --name 'HG002#1' \
  --contignamefix
```

This example fixes headers to start with `HG002#1` and writes `HG002_h1` in
`query_paths.prepared.txt`. Using `--name` with more than one list entry is an
error because one name cannot identify multiple assemblies.

Assembly lists use exactly two columns, `NAME FASTA`. Relative FASTA paths
are resolved from the working directory. Preparation writes adjacent
`FASTA.fai` indexes; a third index-path column is not accepted.
