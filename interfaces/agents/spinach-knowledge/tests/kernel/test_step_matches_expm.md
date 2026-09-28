# tests/kernel/test_step_matches_expm.m

- Signature: `result=test_step_matches_expm()`

## Purpose

Tests Hilbert-space propagation against matrix exponentiation.

## Physical / mathematical content

- Checks the Spinach sign convention for density-matrix evolution: `rho(t)=exp(-iHt) rho(0) exp(+iHt)`.

## Numerical / algorithmic content

- Builds a one-proton spin system in the `zeeman-hilb` formalism with no approximation.
- Sets `H=2*pi*123*S.z`, `rho=S.x+0.25*S.y`, and `dt=2.5e-3` using `S=pauli(2)`.
- Computes `P=expm(-1i*H*dt)` and compares `step(spin_system,H,rho,dt)` with `P*rho*P'` using absolute and relative tolerances of `1e-13`.

## Outputs

- `result` - regression test result with explanatory messages.