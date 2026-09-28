# tests/kernel/test_dynamic_iserstep_highorder.m

- Signature: `result=test_dynamic_iserstep_highorder()`

## Purpose

Regression test of high-order `iserstep` branches on a one-proton Zeeman-Hilbert fixture with a coherent 2×2 density matrix and a nonlinear, non-commuting Hamiltonian.

## Numerical / algorithmic content

- Uses `dt=5e-3`; the zero-time `LG4A` result is checked against the input state at `1e-15` tolerance.
- Compares `LG4A`, `RKMK4`, `RKMK-DP5`, `RKMK-DP8`, and `RKMK-RKF45` with a reference formed from two `RKMK-DP8` half-steps; branch comparison tolerances range from `5e-8` to `5e-10`.
- Checks trace and Hermiticity to `1e-12`.

## Outputs

- `result` — regression-test result with explanatory messages.
