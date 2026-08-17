# BRCA100 NetBID2 network QC benchmark

This benchmark compares cumulative SJARACNe methods on the same BRCA100 input:

1. `12113fb`: pre-PR66 baseline, using the legacy N-out-of-N bootstrap with replacement and legacy affine MI cutoff.
2. `5809183`: PR66, using fixed 80% subsampling without replacement and the legacy affine cutoff evaluated at m=80.
3. `7633ebb`: PR67 stacked on PR66, using fixed 80% subsampling without replacement and the estimator-matched m=80, Npar=40 AP-MI null.

The comparison must keep expression data, TF/SIG lists, 100 resampling seeds,
`p=1e-7`, `Npar=40`, DPI, consensus, compiler, and workflow settings fixed.
NetBID2 QC should be run separately for TF and SIG consensus networks.

NetBID2's plots are descriptive QC, not proof that one inferred network is
biologically more correct. The benchmark should therefore also report edge,
node, and driver counts; edge Jaccard; per-driver target-size correlations;
retained-edge MI/rank correlations; target-size distributions; and consensus
support distributions.

## Isolated environment

The worktree-local environment does not modify the existing ArchR environment.
It pins R 4.4.3 through Micromamba and installs NetBID2 2.2.0 from exact commit
`5defa454d600b94f5dd6d1f9f4428f99759a6821`.

Create the environment from WSL:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
mkdir -p "$ROOT/.micromamba/bin" "$ROOT/.micromamba/home" "$ROOT/.micromamba/cache"
cp /mnt/d/GitHub/Xenium_benchmark/.micromamba/bin/micromamba \
  "$ROOT/.micromamba/bin/micromamba"
"$ROOT/.micromamba/bin/micromamba" create -y \
  -r "$ROOT/.micromamba" \
  -p "$ROOT/.micromamba/envs/netbid2-r44" \
  -f "$ROOT/benchmarks/brca100_netbid_qc/environment-netbid2-r44.yml"
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" Rscript \
  "$ROOT/benchmarks/brca100_netbid_qc/install_netbid2.R"
```

For exact reproduction of the solved Linux packages, use the checked-in
explicit lock instead of re-solving the YAML:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
"$ROOT/.micromamba/bin/micromamba" create -y \
  -r "$ROOT/.micromamba" \
  -p "$ROOT/.micromamba/envs/netbid2-r44" \
  -f "$ROOT/benchmarks/brca100_netbid_qc/environment-netbid2-r44-linux-64.lock"
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" Rscript \
  "$ROOT/benchmarks/brca100_netbid_qc/install_netbid2.R"
```

Verify the pinned package and required network-QC exports:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" Rscript \
  "$ROOT/benchmarks/brca100_netbid_qc/verify_netbid2.R"
```

Exercise the actual NetBID2 importer and HTML QC renderer on its bundled TF
consensus network:

```bash
ROOT=/mnt/d/GitHub/SJARACNe-brca100-netbid-qc
"$ROOT/benchmarks/brca100_netbid_qc/netbid2-r" Rscript \
  "$ROOT/benchmarks/brca100_netbid_qc/smoke_test_netbid2_qc.R"
```

The generated Linux package lock and R session information are recorded after
environment creation. Large environment, work, and QC output directories are
ignored by Git.
