# Single-source DPI skip benchmark

This benchmark measures the fast path that skips DPI when fewer than two source
rows are available. With one source row, SJARACNe has only hub-to-target edges and
cannot form the three-edge triangle required by DPI.

The fixture contains 19,936 genes and one source row with 19,935 neighbors. Its MI
values are strictly decreasing so the previous implementation performs the failed
target-to-target lookup loop instead of taking the equal-MI early exit.

```bash
python3 benchmarks/dpi_skip/generate_single_source_adjacency.py \
  --expression <expression_g19936_n0500.exp> \
  --output <single_source_g19936.adj> \
  --hub-output <source_hub.txt>
```

Both executables were built with GCC 13.3.0 and `-O3` under WSL on an Intel Core
i7-10700F. Runs used `-p 1 -e 0 -s <source_hub.txt>`. Times include reading the
27 MB expression matrix, parsing the adjacency file, DPI, and writing the network.

## Results (2026-08-08)

| Implementation | Wall times (s) | Median (s) |
|---|---|---:|
| Previous DPI loop | 1.85, 1.88, 1.90 | 1.88 |
| Source-row fast path | 1.42, 1.39, 1.49 | 1.42 |

The end-to-end speedup was **1.32x**, saving about 0.46 seconds. The retained
network row was byte-identical before and after the change (SHA-256
`b607275efe4932aa224405ed39d4da916cccdefbdac9b47736eb82d22678d40f`).

This is intentionally not presented as a general SJARACNe speedup. It removes an
impossible DPI pass from single-source networks; expression loading and MI
calculation are unaffected.
