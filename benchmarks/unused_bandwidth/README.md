# Unused marker-bandwidth benchmark

SJARACNe previously called `computeMarkerBandwidth()` before every standard
run. That function recomputed every marker variance, copied and sorted every
marker's observations, calculated an interquartile range, and stored the result
in `Marker::bandwidth`. No MI, DPI, filtering, noise, or output code reads that
field. Without noise correction, the computed variance was also unread; with
noise correction, it duplicated the required explicit variance pass.

The optimization removes only the automatic call from `runStandard()`. The
legacy method and field remain available to avoid an unrelated source/API
change. When noise correction is requested, the existing explicit
`computeMarkerVariance()` call remains in place.

## Time bound

For `G` markers and `N` selected observations, the deleted pass performed a
linear variance scan and an `N`-value sort for every marker:

```text
Theta(G * N * log N) time
Theta(N) transient memory
```

The rest of the algorithm is unchanged. Because this setup pass ran once per
bootstrap, its relative end-to-end benefit decreases as the number of evaluated
hub-target pairs grows.

## Matched results

The baseline is commit `0bfaeb9`, immediately before the removal. Baseline and
candidate were built with GCC 13.3.0 at `-O3` and copied to the WSL-native
filesystem. Inputs were also staged on that filesystem. Runs were serial, with
one warm-up per implementation and three measured repetitions in alternating
`baseline/candidate`, `candidate/baseline`, `baseline/candidate` order.

Wall values below are measured with Python's high-resolution monotonic timer
and are medians of the three runs.

| Case | Baseline | Candidate | Speedup | Time saved |
|---|---:|---:|---:|---:|
| Preprocessing: G=1,000, N=500 | 0.068 s | 0.062 s | 1.10x | 0.006 s |
| Preprocessing: G=5,000, N=1,000 | 0.691 s | 0.606 s | 1.14x | 0.084 s |
| Preprocessing: G=19,936, N=500 | 1.343 s | 1.190 s | 1.13x | 0.153 s |
| AP-MI: H=1, G=5,000, N=1,000 | 1.060 s | 0.978 s | 1.08x | 0.082 s |
| AP-MI: H=100, G=5,000, N=1,000 | 10.530 s | 10.462 s | 1.01x | 0.067 s |
| BRCA100 adjacency replay | 2.555 s | 2.526 s | 1.01x | 0.029 s |

The 1.006x H=100 and 1.012x BRCA replay ratios are within run-to-run variation
and should be treated as no meaningful measurable end-to-end speedup.

The preprocessing cases use an empty adjacency file to isolate loading and the
deleted setup work. AP-MI cases use an unreachable MI threshold so network
writing and DPI do not obscure the one-time cost. The BRCA100 input is an
already DPI-reduced adjacency replay, so it is a loading/setup check rather
than a first-pass network-inference benchmark.

Peak RSS was effectively unchanged, as expected: the deleted calculation used
only `Theta(N)` transient doubles for its sort and quartile buffers. The
candidate/baseline median RSS ratios across the six cases ranged from 0.989 to
0.999.

Machine-readable measurements are in:

- `raw_timings_2026-08-08.csv`
- `results_2026-08-08.csv`

The runner emits `raw_timings.csv` and `summary.csv`; those outputs were renamed
to the dated snapshot names above when recorded in the repository.

## Correctness validation

All compared network-data bytes were identical. A separate retained-edge
validator checks five modes rather than relying on the header-only timing
outputs:

- bootstrap resampling with ties and repeated observations;
- conditioning followed by bootstrap resampling;
- bootstrap resampling with nonzero noise correction;
- all-gene inference;
- BRCA100 adjacency replay.

The final candidate also passed all 66 repository tests. Forty focused tests
passed under combined AddressSanitizer and UndefinedBehaviorSanitizer, including
the fixed MI references and expression, adjacency, DPI, and subnetwork tests.

## Reproduce

Generate the ignored rank-cache inputs as described in
`benchmarks/rank_cache/README.md`, build matched baseline and candidate binaries,
and stage the inputs on a native Linux filesystem. Then run:

```bash
python3 benchmarks/unused_bandwidth/run_benchmark.py \
  --baseline /path/to/baseline/sjaracne.exe \
  --candidate /path/to/candidate/sjaracne.exe \
  --rank-data /path/to/rank-cache/generated \
  --brca-data /path/to/staged/brca-data \
  --output /tmp/sjaracne-unused-bandwidth \
  --repetitions 3
```

For exact retained-network validation:

```bash
python3 benchmarks/unused_bandwidth/validate_outputs.py \
  --baseline /path/to/baseline/sjaracne.exe \
  --candidate /path/to/candidate/sjaracne.exe \
  --fixtures benchmarks/rank_cache/fixtures \
  --brca-data /path/to/staged/brca-data
```
