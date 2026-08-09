# Compact sorted-adjacency benchmark

Measured results are recorded in
[`RESULTS_2026-08-08.md`](RESULTS_2026-08-08.md), with the individual runs in
[`raw_timings.csv`](raw_timings.csv) and medians in
[`summary.csv`](summary.csv).

This harness compares the `std::map` adjacency rows at commit `9f11af6` with
the compact, target-sorted row implementation. It measures the entire SJARACNe
process, not an isolated container microbenchmark.

The workload has three parts:

1. the existing five deterministic `benchmarks/dpi_intersection` adjacency
   replays (three sparse selected-hub cases, one degree-skewed all-gene case,
   and one dense tied-MI control);
2. a replay of the real BRCA100 TF adjacency; and
3. a full BRCA100 seed-1 MI plus DPI run, which exercises adjacency construction
   as well as traversal and serialization.

The BRCA replay is a realistic imported topology, but the tracked
`tests/inputs/adjmat_dir/TF_run.adj` is already DPI-pruned. It must not be
presented as a benchmark of the original pre-DPI graph. The full run provides
the end-to-end construction measurement.

## Protocol

- Compile the baseline from exactly `9f11af6` and the candidate with the same
  compiler and release flags.
- Run one SJARACNe process at a time.
- Put binaries, generated fixtures, expression data, and result files on WSL's
  native ext4 filesystem (for example, under `/tmp`), rather than `/mnt/c` or
  `/mnt/d`.
- Use one warmup per implementation for replay cases. The full BRCA warmup is
  disabled by default because it is expensive.
- Use at least three measured repetitions. First-run order alternates both by
  case and repetition, so neither executable is always run first.
- Compare `network_data_sha256` across every measured baseline and candidate
  output. Header lines beginning with `>` are excluded because they carry run
  metadata; the complete-file hash is still recorded for inspection.

The Python timer records end-to-end `time.perf_counter()` wall time. GNU
`/usr/bin/time` independently records elapsed, user CPU, system CPU, and
maximum resident set size. On Linux, GNU time reports maximum RSS in KiB.

## Generate the existing five synthetic fixtures

No second fixture generator is needed. Reuse the DPI-intersection generator:

```bash
python3 benchmarks/dpi_intersection/generate_benchmarks.py \
  --expression /tmp/compact-adjacency/expression_g05000_n0100.exp \
  --output-dir /tmp/compact-adjacency/fixtures
```

The 5,000-gene expression file can be prepared with the existing rank-cache
benchmark tooling, then copied to `/tmp` with the other inputs.

## Run the matched benchmark

Example using the frozen binaries for this optimization:

```bash
python3 benchmarks/compact_adjacency/run_benchmark.py \
  --baseline /tmp/sjaracne-sorted-row-baseline-9f11af6/sjaracne.exe \
  --candidate /tmp/sjaracne-sorted-row-candidate/sjaracne.exe \
  --synthetic-expression /tmp/compact-adjacency/expression_g05000_n0100.exp \
  --fixtures-dir /tmp/compact-adjacency/fixtures \
  --brca-expression /tmp/compact-adjacency/BRCA100.exp \
  --brca-hubs /tmp/compact-adjacency/BRCA100_TF.txt \
  --brca-adjacency /tmp/compact-adjacency/TF_run.adj \
  --config-dir /tmp/compact-adjacency/config \
  --results-dir /tmp/compact-adjacency/results \
  --warmups 1 --full-warmups 0 --repeats 3
```

For the full BRCA run, the TF and subnetwork lists are both
`BRCA100_TF.txt`, with `p=1e-7`, DPI tolerance `0`, adaptive partitioning,
bootstrap sample `1`, seed `1`, and `N=40`. The adjacency replay uses `p=1`
so the loader retains every already-thresholded edge in the stored network.

The runner writes:

- `run_metadata.json`: paths and SHA-256 hashes of both binaries;
- `raw_timings.csv`: every measured timing, RSS value, command, and output hash;
- `summary.csv`: per-case median wall/CPU/RSS, speedup, RSS ratio, and RSS
  reduction; and
- `logs/`: standard output, standard error, and raw GNU-time records.

Network outputs are deleted after hashing unless `--keep-outputs` is supplied.
The runner stops immediately if a command fails or DPI does not execute. It
fails before writing a successful summary if any measured network-data hash
differs between the two implementations.
