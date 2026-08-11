# Adaptive DPI neighbor-intersection benchmark

This benchmark measures the DPI traversal used after an adjacency matrix has
been created or loaded with `-j`. It targets sparse multi-hub and degree-skewed
all-gene networks, where most pairs of neighbors around a source gene do not
close a triangle.

## Optimization

For one source gene `A`, the legacy loop sorts `A`'s neighbors by decreasing MI.
For each candidate edge `A-B`, it then probes every sufficiently stronger
`A-C` edge until it finds a qualifying `B-C` edge. Sparse graphs therefore spend
most of their DPI time asking for edges that do not exist.

The optimized loop builds a reverse-neighbor index once, then intersects:

1. the prefix of `A`'s neighbors whose `MI(A,C)` exceeds the DPI threshold for
   `A-B`; and
2. the effective neighbors of `B`.

It scans whichever side is smaller. When the effective-neighbor side is used,
the original position of every `A` neighbor is retained so the selected
intermediate remains the first one in the legacy MI-descending order.

Imported directed adjacency files have an unusual legacy rule that is
deliberately preserved: a nonempty `B` row is authoritative, even if it lacks
`B-C`; reverse `C-B` lookup is used only when `B`'s row is absent or empty.
Strict DPI comparisons, TF protection, and the use of already-pruned edges as
support for later decisions are also unchanged.

## Complexity

Let `s_AB` be the number of stronger eligible neighbors of `A` for candidate
`A-B`, and let `Eff(B)` be the direct row of `B`, or its incoming row when the
direct row is empty. Index construction costs `O(G + E)` time and memory. After
the existing per-row MI sort, the optimized traversal examines approximately

```text
sum_A sum_{B in N(A)} min(s_AB, |Eff(B)|)
```

candidate entries. Prefix-side probes still use `std::map` lookups and retain
their logarithmic factor.

For a simple sparse symmetric graph, the full-neighborhood analogue is the
standard `O(E^(3/2))` smaller-side triangle-enumeration bound, compared with
work approaching `sum_A d_A^2` in the old wedge loop. For a directed hub network
with `H` computed source rows and `T` target-only genes, the intended sparse
case changes roughly from `O(H*T^2)` failed probes to `O(T*H^2)` intersection
work when `H << T`.

This is not `O(number of triangles)`: intersecting two sets can still inspect
nonmatching entries, and dense or tied-MI graphs may show no improvement.

## Generated fixtures

The large adjacency fixtures are deterministic but are not tracked. Generate
them from a 5,000-gene expression panel:

```bash
python3 benchmarks/dpi_intersection/generate_benchmarks.py \
  --expression benchmarks/rank_cache/generated/expression_g05000_n0100.exp \
  --output-dir /tmp/sjaracne-dpi-intersection-fixtures
```

The generator creates three selected-hub cases (`H=10, 50, 100`, each with
4,000 target-only genes), one degree-skewed all-gene case, and one dense tied-MI
control. Keeping these generated files outside Git avoids committing hundreds
of megabytes of derivative adjacency data.

Run matched legacy and candidate binaries with alternating order and three
repetitions:

```bash
python3 benchmarks/dpi_intersection/run_benchmark.py \
  --baseline /path/to/legacy/sjaracne.exe \
  --candidate /path/to/candidate/sjaracne.exe \
  --expression benchmarks/rank_cache/generated/expression_g05000_n0100.exp \
  --fixtures-dir /tmp/sjaracne-dpi-intersection-fixtures \
  --results-dir /tmp/sjaracne-dpi-intersection-results \
  --repeats 3
```

The runner uses `-p 1 -e 0`, verifies that DPI actually runs, and requires the
network data SHA-256 to agree across every baseline and candidate run.

For broader semantic validation, the differential runner creates random
symmetric and directed networks, including empty rows, MI ties, TF lists, hub
subsets, and several DPI tolerances:

```bash
python3 benchmarks/dpi_intersection/differential_fuzz.py \
  --baseline /path/to/legacy/sjaracne.exe \
  --candidate /path/to/candidate/sjaracne.exe \
  --cases 500 --seed 20260808
```

## Results (2026-08-08)

Both executables were compiled with GCC 13.3.0 and `-O3` under WSL on an Intel
Core i7-10700F. Times are end-to-end wall times, including expression loading,
adjacency parsing, DPI, and output writing. Values are medians of three runs.

| Case | Directed edges | Legacy (s) | Intersection (s) | Speedup |
|---|---:|---:|---:|---:|
| 10 hubs, 4,000 targets | 40,040 | 0.296 | 0.103 | 2.88x |
| 50 hubs, 4,000 targets | 200,200 | 1.135 | 0.241 | 4.71x |
| 100 hubs, 4,000 targets | 400,400 | 2.186 | 0.443 | 4.94x |
| Degree-skewed all-gene | 800,000 | 6.941 | 1.212 | 5.73x |
| Dense tied all-gene control | 249,500 | 0.233 | 0.236 | 0.99x |

All 30 measured outputs had identical network-data hashes before and after the
optimization. A separate deterministic differential run matched all 500 random
cases. The dense tied control was about 1.5% slower, which is within the cost of
building the index and confirms that the optimization should not be described
as a universal DPI speedup. See `raw_timings.csv` and `summary.csv` for exact
measurements and hashes.
