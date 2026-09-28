# tests/kernel/test_multiprop_adaptive.m

- Signature: `result=test_multiprop_adaptive()`

## Purpose

Tests adaptive repeated propagator application by comparing `multiprop()` with explicit matrix-power references. Returns a regression test result with explanatory messages.

## Numerical / algorithmic content

- Checks binary adaptive squaring for non-normal sparse and diagonal sparse propagators acting on state vectors, and for a square Liouville-space vector stack. These states propagate by left multiplication, `P^N*rho`.
- Checks that zero propagator applications leave the state unchanged and that a one-dimensional wavefunction follows the vector branch.
- Checks Hilbert-space density-matrix propagation for unitary, sparse non-unitary, and diagonal sparse propagators against `P^N*rho*(P^N)'`.
- Checks that `prop_chop` removes small elements generated during propagator squaring, and that row vectors are rejected as invalid Spinach state vectors.

## Outputs

- `result` — regression test result with explanatory messages.