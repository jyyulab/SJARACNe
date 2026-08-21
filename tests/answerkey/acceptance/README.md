# Acceptance fixtures

`cnn_5.txt` is the active end-to-end acceptance result. It is generated from
five fixed-size 80% subsamples without replacement.

`sjaracne_workflow.yml` is the former five-run legacy-bootstrap job. The
remaining parameter, bootstrap, and consensus files are historical snapshots
from a 100-run legacy bootstrap-with-replacement analysis. None of those legacy
files is used by the current acceptance test or describes the new workflow
default.
