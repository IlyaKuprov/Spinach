# tests/kernel/test_step_zero_time.m

- Signature: `result=test_step_zero_time()`

## Purpose

Tests zero-duration propagation: a propagator over zero time must leave the density matrix unchanged.

## Numerical / algorithmic content

A one-proton Hilbert-space spin system is constructed. The test sets `rho=S.x+2*S.z` and `H=3*S.x+5*S.z`, then calls `step(spin_system,H,rho,0)`. It compares the propagated state with `rho` using `test_close` with absolute and relative tolerances of `1e-15`.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announce and register the zero-duration propagation test.
- Build a one-proton spin system with `zeeman-hilb` formalism and no approximation.
- Propagate the density matrix for zero time and check the identity limit.