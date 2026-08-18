# PR67 threshold-sweep analysis contract

`analyze_sweep.py` analyzes only the PR67 per-subsample p-value sweep. It does
not choose a threshold and does not perform biological-reference or downstream
NetBID reproducibility validation.

## Required work-root layout

```text
<work-root>/
  inputs/
    BRCA100.exp
    BRCA100_TF.txt
    BRCA100_SIG.txt
  results/
    run_manifest.tsv                         # preferred aggregated seed evidence
    <p_key>/
      point_manifest.json                    # authoritative p metadata
      tf/
        seed_metadata/*.json                 # fallback seed evidence
        consensus/consensus_network_ncol_.txt
        consensus/consensus_support.tsv
      sig/
        seed_metadata/*.json
        consensus/consensus_network_ncol_.txt
        consensus/consensus_support.tsv
```

Each `point_manifest.json` must contain `p_value`. The analysis also preserves
`p_token`, `mi_cutoff`, `validation_class`, `tail_extrapolated`, commit, and null
model hash fields when present. The folder name is never parsed as the p-value.

Seed evidence may instead be supplied by a per-point or per-arm
`run_manifest.tsv`, or by a nonempty `seeds` list in `point_manifest.json`. All
TF/SIG and threshold arms must resolve to exactly the same seed set. When seed
JSON commands contain `-p`, every command is checked against the corresponding
point manifest.

All directed-edge endpoints and candidate drivers are checked against the
fixed expression matrix. This supplies the denominator for the reported
incident-node fraction.

Run:

```bash
python benchmarks/brca100_pr67_threshold_sweep/analyze_sweep.py \
  --work-root /path/to/brca100-pr67-threshold-sweep
```

The default output is `<work-root>/results/analysis`. Use
`--pr66-work-root /path/to/prior-brca100-work-root` to add PR66 horizontal
references and PR66-to-sweep overlap tables. PR66 is context only; it is not a
sweep point. The default within-sweep anchor is p=1e-7 and can be changed with
`--anchor-p`.

Machine-readable outputs include:

- `network_summary.tsv`: edge count, fixed-list driver coverage, weak-component
  connectivity, target sizes, support, and consensus-MI summaries.
- `adjacent_overlap.tsv`: each p against the next larger (less stringent) p.
- `anchor_overlap.tsv`: every p against the selected within-sweep anchor.
- `point_manifest_summary.tsv` and `seed_manifest_summary.tsv`: provenance and
  matched-design checks.
- Optional `pr66_context_summary.tsv` and `pr66_context_overlap.tsv`.

Plots are emitted as deterministic PNG and SVG files under `analysis/plots`.
The x-coordinate is `log10(p)`; manifest labels are used only as tick labels.
