# Adaptive-partitioning MI rank-cache benchmark data

This benchmark isolates the cost of pairwise adaptive-partitioning mutual
information (AP-MI). It is intended to measure the proposed change from sorting
the same gene once per candidate edge to ranking each gene once per bootstrap.
It is not a DPI benchmark.

## Storage policy

The generator, specification, and small tie-heavy correctness fixture are
tracked in Git. The large derived matrices are written to `generated/`, which is
ignored by Git. They should not be added to `tests/inputs`: they are performance
data, not unit-test inputs, and committing multiple projections of the same
source matrix would permanently inflate repository history.

The generator never modifies the source expression file. Its output includes a
manifest with SHA-256 checksums, dimensions, selection hashes, and provenance.

## Source data

The default specification was designed for the Conservative T&S all-gene
matrix from the aging-ovary analysis:

- 19,936 genes
- 1,000 cells
- 1,681 genes overlapping `SJARACNe/config/TF_list.txt`
- FOXP1 is forced to be the first hub
- CDKN1A is forced into every gene panel

The matrix is approximately 90% zero, including 2,523 all-zero genes and 104
all-zero eligible TFs. They are retained deliberately because ties and sparse
expression are part of the workload that the optimization must handle.

The source matrix remains outside this repository because it is a derived public
dataset and is already present in the analysis directory.

## Generate the panels

From Git Bash with the repository as the working directory:

```bash
python benchmarks/rank_cache/prepare_rank_benchmarks.py \
  --expression "E:/analysis/Aging/Ovarian/Wu_et_al_2024_Nature_Aging/results/regulation/02_sjaracne_foxp1_cdkn1a_gcs_ts_all_genes/inputs/Conservative_TS/Conservative_TS.exp.txt"
```

The command refuses to replace existing artifacts unless `--force` is supplied.
Large outputs are created under `benchmarks/rank_cache/generated/`.
For the audited source, the complete generated set is 64.95 MiB (seven
expression matrices plus small hub lists and metadata).

## Factor-at-a-time design

The generator writes only seven unique expression matrices, then reuses them
across 12 unique benchmark cases:

| Sweep | Fixed dimensions | Varied dimension |
|---|---|---|
| Hub count | G=5,000, N=1,000 | H=1, 10, 50, 100, 500, 1,681 |
| Observation count | G=5,000, H=100 | N=100, 250, 500, 1,000 |
| Gene count | N=500, H=100 | G=1,000, 5,000, 10,000, 19,936 |

Shared cases are listed once in `benchmark_cases.csv`. Gene panels and cell
panels are nested. Selection uses a fixed SHA-256 ordering, so regeneration from
an unchanged source/specification produces the same matrix and hub-list files
on every platform. The provenance manifest also records absolute local paths,
so the manifest itself is machine-specific.
Hash ordering also avoids taking the first N source columns, which would be
confounded because the source cells occur in contiguous donor blocks.

For the large timing runs, use a deliberately unreachable MI threshold and
disable DPI so output construction does not hide AP-MI timing. Build the native
Linux executable in WSL first; the executable is generated locally and is not
tracked by Git. From Git Bash:

```bash
MSYS_NO_PATHCONV=1 wsl.exe -e bash -lc 'cd /mnt/d/GitHub/SJARACNe && \
  make -C SJARACNe && \
  SJARACNe/bin/sjaracne.exe \
    -i benchmarks/rank_cache/generated/expression_g05000_n1000.exp \
    -s benchmarks/rank_cache/generated/hubs_h0100.txt \
    -S 17 -r 1 -t 100 -e 1 \
    -o /tmp/rank-cache-baseline.adj'
```

Run each case several times and report wall time, CPU time, and peak RSS. Do not
start with `hubs_h1681.txt`: that case evaluates 8,403,319 candidate pairs and is
intentionally expensive.

## Correctness fixture

`fixtures/tied_counts.exp` is a small synthetic matrix with many tied values.
Running it with bootstrap resampling also creates repeated observations. Use it
with `fixtures/tied_hubs.txt`, a fixed seed, `-t 0`, and `-e 1` to compare MI
values and edge sets before and after the rank-cache refactor.
`fixtures/tied_seed17_reference.tsv` records the current Linux baseline for
`-S 17 -r 1`. The regression test compares numerical values with an explicit
floating-point tolerance. A second reference covers conditioning followed by
bootstrap resampling, where confusing original cell IDs with resample positions
would change ranks. Output metadata contains paths and is not byte-compared
wholesale. Keep the build and C runtime fixed when comparing bootstrap seeds
because `std::rand()` is not portable across C runtimes.
