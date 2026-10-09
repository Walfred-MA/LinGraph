# SmallVariantMerge

C++17 workers for the existing graph-VCF SNP/indel pipeline. Build with `make`
(zlib development files required); the repository installer also builds and
copies the executable into `scripts/`. The Python wrapper can compile a missing
binary into a source-hashed, locked temporary cache.

Use the established commands:

```sh
python scripts/graphvcfmerge.py --snp -i sample1.vcf sample2.vcf -o cohort.snp.vcf -t 32
python scripts/graphvcfmerge_snp.py -i sample1.vcf sample2.vcf -o cohort --svcutoff 20 -t 32
```

The second command emits exact SNP and small-indel unions. Existing staged
scan/chrom/concat commands and protocol 2–4 SNP shards remain supported. The
worker's private binary request protocol is written by `graphvcfmerge_native.py`.
SV clustering/alignment and SV-derived indels retain the existing implementation.

## Memory and work

- Coordinates and dictionary/site indices are `uint32_t`; byte offsets, record
  totals, memory accounting and coverage endpoints use 64 bits.
- Packed SNP shards remain 21 bytes; aligned working records are 24 bytes.
- SNP sorting is five stable radix passes over `(POS, REF/ALT)`, using the
  requested threads. Grouping preserves the Python observation sort/dedup order.
- Half the RAM budget bounds records and radix workspace. The remaining half
  allows for dictionaries, the Python parent, grouping and output. Larger inputs
  spill sorted runs with a maximum fan-in of 64. Small-indel grouping uses exact
  keys, bounded observation runs and a site metadata table; it rejects metadata
  exceeding its reserved quarter of the budget.
- Genotype defaults live in one reusable output buffer. Coverage changes patch
  the affected samples; variant observations and explicit calls replace only
  their own fields. Computation avoids a full cohort scan at each site. The
  dense VCF still requires O(sites × samples) output bytes.
- SNP emission reuses buffers for touched samples and formats directly from
  interned strings. Duplicate observations retain the Python tuple ordering.
  Oversized per-sample buffers are released after high-copy sites, so isolated
  large sites cannot accumulate large caches across the cohort.
- Indel emission consumes the sorted sample/allele groups directly, removing
  adjacent duplicates without rebuilding maps and sets. The site lookup table
  is released before output, and legacy label/offset parsing uses linear scans.
- Defaults: 64 GiB below 100 samples, 128 GiB through 2,000; larger cohorts scale
  in 64-GiB increments to 512 GiB at 10,000. Override with `--small-ram-gb` or
  `GRAPHVCFMERGE_RAM_GB`. Slurm allocation limits take precedence. Concurrent
  local small-indel chromosome workers share the configured budget.

The native parser preserves custom escaping, supported legacy FORMAT schemas,
parallel FORMAT vector cardinalities, joined/negative labels, exact indel keys,
stable IDs, filters, missing/reference calls and assembly provenance. Original
chromosome realignment additions and removals are applied in C++ before sorting.
Completed parts replace their predecessors only after a successful worker exit.

## Comparison

`GRAPHVCFMERGE_BACKEND=python` selects the retained reference implementation.
Run `pytest tests/test_native_small_merge.py` for differential checks and
`python benchmark/benchmark_small_merge.py` for a reproducible 2,000-sample
synthetic comparison. The benchmark checks output hashes and measures both
elapsed time and peak memory in separate worker processes. It is a chromosome
merge benchmark, not an end-to-end cohort or storage benchmark.
Add `--kind indel` to compare the exact small-indel merge using text shards.
Use `--native-binary /path/to/SmallVariantMerge` to compare a saved executable
against the current build on the same synthetic workload.
