# tests/kernel/test_dynamic_equilibrium_frontends.m

- Signature: `result=test_dynamic_equilibrium_frontends()`

## Purpose

Tests the `thermalize()`, `steady()`, and `residual()` dynamic front ends on compact Liouville-space systems. Returns a regression test result with explanatory messages.

## Physical / mathematical content

- In a one-spin spherical-tensor Liouville-space system, the IME thermalisation check verifies that the requested equilibrium state is stationary: `R_ime*rho_eq = 0`.
- The steady-state checks use a constructed contractive affine propagator `P` with a known fixed point `rho_ss`, verified by `P*rho_ss = rho_ss`.
- A heteronuclear `1H`–`13C` system with a weak order matrix tests the effect of residual ordering on a coupling tensor.

## Numerical / algorithmic content

- `thermalize()` in DiBari mode is compared with the defining product `R*propagator(spin_system,H_left,1i*beta)`, where `beta = hbar/(kbol*temperature)` and `temperature = 300.0`.
- `steady()` is called in both `'newton'` and `'squaring'` modes; each returned state is compared with the constructed fixed point.
- After `residual()` is applied, the coupling tensor is checked for preservation of its isotropic trace, equality of its `xx` and `yy` entries, and removal of off-diagonal entries.

## Outputs

- `result` — regression test result with explanatory messages. The fixed-point and tensor comparisons are recorded with tolerances of `1e-14` for the IME, DiBari, and constructed-fixed-point checks; `1e-12` for the Newton result; `1e-10` for the squaring result; and the tolerances specified by the individual residual-tensor checks.

## Implementation structure

- Build a one-spin spherical-tensor Liouville-space system for the thermalisation and steady-state checks.
- Check the IME stationary-state condition and the DiBari relaxation–propagator product.
- Construct and verify a contractive affine fixed point, then compare the Newton and squaring `steady()` results with it.
- Build a heteronuclear system with coordinates and weak order, apply `residual()`, and check the coupling tensor's trace, axial `xy` degeneracy, and zero off-diagonal components.