# test_clexpn_experiments

Syntax: `result=test_clexpn_experiments()`

Tests independent dense matrix propagation for free evolution, Hahn, CPMG and ESEEM, zero-time convention, bath pump boundaries, complex receiver phase and column-major reduced density.

Finite-drive regressions include the two-spin zero-drift z-drive acceptance probe, transverse preparation, Hilbert/Liouville equivalence, overlapping and partial intervals, detuning, fractional placement, noncommuting coupled drift and simultaneous GKSL relaxation. Unsupported ideal events and projected drives are rejected.

Preparation selectors verified for the example migration: `transition.ms=[-0.5 0.5]` for E4, `[-1 1]` for E3, `[-0.5 0.5;-0.5 0.5]` for two electrons, and `[-0.5 0.5;3.5 3.5]` for E/51V. `transition.spins` follows the observed-spin order. These select the equal two-manifold superpositions without manually assembled vectors.
