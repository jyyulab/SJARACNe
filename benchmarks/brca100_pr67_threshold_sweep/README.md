# BRCA100 PR67 per-subsample threshold sweep

This benchmark varies only PR67's seed-level AP-MI `-p` value.  It keeps the
following fixed:

- PR67 commit `7633ebb4a0d966dbda15a4e32d0efa492fb71aeb`;
- estimator-matched model `m=80`, `Npar=40`, SHA-256
  `e3a8522682a8ea239821aaa10b12db72d00e07bfdcad43599d8e76a06be80944`;
- fixed 80% sampling without replacement;
- BRCA100 expression, TF and SIG lists;
- seeds 1 through 100;
- adaptive-partitioning MI, `Npar=40`, and DPI epsilon 0; and
- consensus-network `p=1e-5`.

Biological-reference validation is intentionally deferred.  This sweep asks
whether a less stringent seed-level operating point restores usable driver
coverage and graph connectivity while retaining reasonable consensus support.
The results are descriptive network QC, not proof of biological accuracy.

## Grid

| Seed-level p | AP-MI cutoff | Calibration status | Purpose |
|---:|---:|---|---|
| `1e-7` | 0.322465 | GPD extrapolation | original PR67 default |
| `1e-6` | 0.288001 | GPD extrapolation | coarse grid |
| `1e-5` | 0.247951 | GPD extrapolation | coarse grid |
| `2e-5` | 0.234669 | directly tested; validation boundary | first directly tested point |
| `5e-5` | 0.216162 | directly tested on held-out stream | coarse grid |
| `1e-4` | 0.201409 | directly tested on held-out stream | coarse grid |
| `2e-4` | 0.185973 | directly tested on held-out stream | coarse grid |
| `3e-4` | 0.176615 | interpolated inside accepted validation range | near PR66 density anchor |
| `3.528045626e-4` | 0.172803 | interpolated inside accepted validation range | exact PR66-cutoff match |

The final probability is calculated from the PR67 GPD model so that its MI
cutoff equals PR66's legacy affine cutoff at `m=80`.  PR66 remains context, not
ground truth.

Every point is inferred directly with `-p ... -M ...`.  The benchmark does not
filter or replay a looser adjacency matrix: stored adjacency MIs have limited
decimal precision, and DPI is applied after seed-level thresholding.

## Persistent execution

Run from WSL using the already pinned NetBID2 environment.  The work root is on
the persistent native-ext4 home filesystem, not `/tmp`.

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
RUN="$ROOT/benchmarks/brca100_netbid_qc/netbid2-r"
WORK="$HOME/sjaracne-benchmarks/brca100-pr67-threshold-sweep-20260818"
SWEEP="$ROOT/benchmarks/brca100_pr67_threshold_sweep/run_sweep.py"

mkdir -p "$WORK"

"$RUN" python "$SWEEP" --phase prepare --work-root "$WORK"

# Small direct smoke test before the full sweep.
"$RUN" python "$SWEEP" --phase infer \
  --points p1e-07,p2e-05,p_pr66_cutoff_match \
  --drivers tf --seed-start 1 --seed-end 1 --workers 3 \
  --work-root "$WORK"

# Full resumable run. Existing smoke outputs are hash-validated and skipped.
"$RUN" python "$SWEEP" --phase infer \
  --points all --drivers all --seed-start 1 --seed-end 100 --workers 12 \
  --work-root "$WORK"

# Equivalent fixed launcher used for the recorded long inference run:
bash "$ROOT/benchmarks/brca100_pr67_threshold_sweep/run_full_inference.sh"

# Consensus is serialized and keeps the seed-level p fixed per arm.
"$RUN" python "$SWEEP" --phase consensus \
  --points all --drivers all --seed-start 1 --seed-end 100 --workers 1 \
  --work-root "$WORK"
```

The runner writes an immutable `sweep_design.json`, one `point_manifest.json`
per probability, checksummed seed metadata, a cumulative `run_manifest.tsv`,
consensus manifests, and an invocation history.  Completed jobs are accepted on
resume only after their command fingerprint, headers, structure, and adjacency
SHA-256 are revalidated.

## Interpretation

Selection should prioritize the smallest, more stringent directly validated
`p` at or just beyond the joint TF/SIG coverage-connectivity knee.  At minimum,
report consensus edges, active-driver fraction, zero-filled target sizes,
incident-node fraction, weak components, largest-component fraction, MI,
support, and adjacent-threshold overlap.  Scale-free adjusted R-squared alone
must not select the operating point: the original extremely sparse PR67 SIG
network had a high value despite severe driver and connectivity loss.

The sweep does not establish biological optimality, empirical FDR, or
independent-cohort NetBID reproducibility.  Those require a later biological
validation phase.
