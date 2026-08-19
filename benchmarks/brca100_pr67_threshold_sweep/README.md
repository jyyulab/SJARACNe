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

## Results

The completed, independently audited sweep is summarized in
[results_2026-08-19/RESULTS.md](results_2026-08-19/RESULTS.md). It nominates
`p=3e-4` only as a provisional BRCA100 topology operating point and recommends
a focused refinement between `2e-4` and `3e-4` before any default change.

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

Here, Gate 2 means validation against a second, independent stream of
canonical rank permutations generated under the AP-MI independence null.  It
is not a held-out BRCA cohort, a biological-reference test, or downstream
NetBID reproducibility.

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

# Recompute the data-section hashes and prove, for all 400 anchor networks,
# that p=1e-7 reproduces the prior PR67 arm and that the cutoff-match point
# reproduces PR66.
"$RUN" python \
  "$ROOT/benchmarks/brca100_pr67_threshold_sweep/validate_anchor_equivalence.py" \
  --sweep-work-root "$WORK" \
  --prior-work-root \
    "$HOME/sjaracne-benchmarks/brca100-netbid-qc-20260817-rerun"

# Consensus is serialized and keeps the seed-level p fixed per arm.
"$RUN" python "$SWEEP" --phase consensus \
  --points all --drivers all --seed-start 1 --seed-end 100 --workers 1 \
  --work-root "$WORK"

# Reconstruct retained-edge support from the same 100 seed networks.
"$RUN" python \
  "$ROOT/benchmarks/brca100_pr67_threshold_sweep/run_support_summaries.py" \
  --benchmark-repo "$ROOT" --work-root "$WORK"

# Import every consensus network into NetBID2 and write compact,
# machine-readable summaries and zero-filled driver target-size tables.
"$RUN" python \
  "$ROOT/benchmarks/brca100_pr67_threshold_sweep/run_netbid_qc.py" \
  --points all --drivers all --html-points none --work-root "$WORK"

# Produce matched topology/support tables and deterministic plots, with PR66
# shown only as context.
"$RUN" python \
  "$ROOT/benchmarks/brca100_pr67_threshold_sweep/analyze_sweep.py" \
  --work-root "$WORK" \
  --pr66-work-root \
    "$HOME/sjaracne-benchmarks/brca100-netbid-qc-20260817-rerun"

# Alternatively, run this fixed launcher instead of the five downstream
# commands above. Do not run both sequences against the same completed root.
bash \
  "$ROOT/benchmarks/brca100_pr67_threshold_sweep/run_full_downstream.sh"

# Recorded HTML reports cover the selected point, its stricter neighbor, and
# the PR66-cutoff match. All 18 summaries are revalidated/reused first; optional
# HTML provenance is written separately from the stable summary record.
SELECTED_HTML_POINTS="p2e-04,p3e-04,p_pr66_cutoff_match"
"$RUN" python \
  "$ROOT/benchmarks/brca100_pr67_threshold_sweep/run_netbid_qc.py" \
  --points all --drivers all --html-points "$SELECTED_HTML_POINTS" \
  --work-root "$WORK"
```

The runner writes an immutable `sweep_design.json`, one `point_manifest.json`
per probability, checksummed seed metadata, a cumulative `run_manifest.tsv`,
consensus manifests, and an invocation history.  Completed jobs are accepted on
resume only after their command fingerprint, headers, structure, and adjacency
SHA-256 are revalidated.

The NetBID2 runner discovers points only from immutable
`results/<point>/point_manifest.json` files. Compact artifacts live under each
arm's `netbid2_qc/`; optional reports live separately under
`netbid2_qc_html/`. Both modes use input-and-environment fingerprints, pending
manifests, atomic directory promotion, and complete file inventories for safe
resume and narrow-window crash recovery. The all-arm summary run writes stable
provenance to `results/netbid2_qc_manifest.json`. A later HTML run requires that
exact complete summary aggregate, never rewrites it, and writes optional-report
provenance to `results/netbid2_qc_html_manifest.json` instead.

## Compact results package

After the full analysis succeeds, create a reviewable package outside the live
work root.  The destination and its sibling `.partial` path must not exist.

```bash
PACKAGE="$ROOT/benchmarks/brca100_pr67_threshold_sweep/results_2026-08-19"
"$RUN" python \
  "$ROOT/benchmarks/brca100_pr67_threshold_sweep/package_results.py" \
  --work-root "$WORK" --output-root "$PACKAGE"
```

The packager requires the exact 9-point by 2-driver design, 1,800 completed
seed runs, both completed full-run invocations, 400 anchor comparisons, all
consensus/support artifacts, all 18 NetBID2 summaries, and a completed analysis
with exact PR66 cutoff-match evidence.  It copies the compact analysis tables
and plots, anchor evidence, and immutable design/build/point/arm/aggregate
manifests.  Optional HTML manifests are included only when their root aggregate
is complete.  Each optional per-arm record must reproduce the exact
`run_netbid_qc.py` input fingerprint and command contract, and its three shared
TSVs must be byte-identical to the already validated stable summary outputs.

Raw adjacencies, full consensus networks, support tables, NetBID2 data products,
and HTML reports are never copied.  `omitted_artifacts.json` preserves their
validated hashes and sizes.  The packager rereads every omitted regular file;
it rejects a current SHA-256 mismatch and also rejects a byte-count mismatch
where the producing manifest records one.  `.gitattributes` disables text
conversion, and `SHA256SUMS` covers every package file except itself.  The
script intentionally
does not invent a `RESULTS.md`; add the reviewed numerical report later and
regenerate `SHA256SUMS` with the `write_sha256s()` helper in
`package_results.py`.

## Interpretation

The pragmatic topology screen was declared before inspecting the intermediate
sweep points, but after observing the PR66 and PR67 endpoint networks; it is
therefore endpoint-informed, not a neutral preregistration.  It nominates the
smallest, most stringent `p` in the held-out validation range for which **both**
TF and SIG consensus networks have
at least 90% active candidate drivers, at least 70% of all expression nodes
incident to an edge, and at least 95% of incident nodes in the largest weak
component.  These are engineering non-collapse floors, not statistically or
biologically optimal cutoffs.  If no point passes, the sweep reports that no
operating point was found rather than relaxing the floors after inspecting the
results.  Any nominated value is a provisional topology operating point.  The
coverage-connectivity knee, adjacent-threshold overlap, and support are
secondary diagnostics rather than additional pass/fail gates.

At minimum, report consensus edges, active-driver fraction, zero-filled target
sizes, incident-node fraction, weak components, largest-component fraction,
MI, support, and adjacent-threshold overlap.  Scale-free adjusted R-squared
alone must not select the operating point: the original extremely sparse PR67
SIG network had a high value despite severe driver and connectivity loss.
Support summaries are conditional on an edge surviving the consensus filter;
they are not edge-inclusion probabilities or uncertainty intervals.  Plot IQR
bands describe distributions across the retained edges, not confidence
intervals or variation across independent benchmark datasets.

The sweep does not establish biological optimality, empirical FDR, or
independent-cohort NetBID reproducibility.  Those require a later biological
validation phase.
