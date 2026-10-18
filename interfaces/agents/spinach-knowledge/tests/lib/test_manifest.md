# tests/lib/test_manifest.m

`manifest=test_manifest()` returns the central structure array of regression identifiers, display names and test function handles used by `run_tests` and `run_test`. Tests must be registered here to participate in those entry points.

The `kernel/flux_relaxation` entry invokes `test_flux_relaxation`, including activation, independent transfer, nonzero Redfield additivity and zero-flux Hilbert checks.
