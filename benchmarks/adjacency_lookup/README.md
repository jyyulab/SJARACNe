# Accession lookup benchmark

This benchmark measures the native `-j` path that reloads an existing adjacency
matrix. It isolates the effect of replacing repeated linear accession searches with
one expression-matrix accession index.

## Complexity

Let:

- `G` be the number of expression-matrix genes;
- `R` be the number of adjacency source rows; and
- `E` be the number of retained adjacency edges whose target labels must be
  resolved.

The previous implementation scanned up to all `G` markers for every label, so
label translation cost `O((R + E) * G)` string comparisons in the worst case. The
new implementation builds the hash index once in `O(G)` expected time and uses
`O(1)` expected lookup per label, reducing expected translation cost to
`O(G + R + E)`. The index requires `O(G)` additional memory.

Duplicate accessions retain the legacy behavior: the first expression row wins.

## Workload

The three fixtures hold the retained-edge count near 100,000 while increasing the
expression-matrix gene count. Each source row connects to every other gene. The
generator writes fixtures outside the repository; the large adjacency files are not
tracked.

```bash
python3 benchmarks/adjacency_lookup/generate_dense_adjacency.py \
  --expression <expression.exp> --sources <source-count> \
  --output <fixture.adj>
```

The matched runner alternates baseline/candidate order, records three repetitions,
and verifies that all network data rows have the same SHA-256 digest:

```bash
python3 benchmarks/adjacency_lookup/run_benchmark.py \
  --baseline <baseline-sjaracne.exe> \
  --candidate <candidate-sjaracne.exe> \
  --expression-dir <rank-benchmark-generated-dir> \
  --adjacency-dir <generated-adjacency-dir> \
  --results-dir <results-dir> --repeats 3
```

Both builds used GCC 13.3.0 with `-O3` under WSL on an Intel Core i7-10700F. The
timed command used `-p 1 -e 1`, retaining all supplied edges and disabling DPI. The
measurement is end-to-end: expression loading, adjacency parsing/indexing, network
storage, and output serialization are all included.

## Results (2026-08-08)

| Genes | Source rows | Retained edges | Linear median (s) | Hash median (s) | Speedup |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 100 | 99,900 | 0.177 | 0.130 | 1.365x |
| 5,000 | 20 | 99,980 | 0.864 | 0.401 | 2.156x |
| 19,936 | 5 | 99,675 | 3.445 | 1.415 | 2.434x |

These numbers are specific to adjacency loading. They do not imply a comparable
speedup for MI calculation or DPI, which are separate execution paths. See
`raw_timings.csv` for every measurement and output-data digest.
