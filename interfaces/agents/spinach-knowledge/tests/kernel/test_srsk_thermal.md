# tests/kernel/test_srsk_thermal.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_srsk_thermal.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_srsk_thermal.m)

## Purpose

Regression test for once-only thermalisation of additive SRSK (scalar-relaxation-of-the-second-kind) relaxation. The test verifies that SRSK adds zero-destination rates before final thermalisation, so that the total generator is thermalised exactly once.

## Behaviour

The test builds a two-spin system (`1H`, `14N`) with a fast quadrupolar source (`r1`/`r2` rates of `1e5` rad/s) coupled to a proton via scalar coupling. The scalar coupling is swept over `100`, `-100`, and `0` Hz, covering positive, negative, and zero coupling. An oriented quadrupole interaction (Euler angles `[0.3 0.7 0.2]`) supplies a complex, noncommuting thermalisation Hamiltonian as a contrary control; the test asserts that this Hamiltonian has a non-negligible imaginary Frobenius norm and does not commute with the relaxation superoperator.

For each coupling value, the test:

- Checks the isotropic `relaxation` call (no orientation argument) against a `thermalize` reference for both `IME` and `dibari` equilibrium methods, requiring exact single thermalisation.
- Computes the analytic SRSK rate augmentation from the scalar-coupling formula: `rate_long = (4/3)*(2*pi*J)^2*1e-5/(1+delta_omega^2*1e-10)` and `rate_trans = (2/3)*(2*pi*J)^2*1e-5 + rate_long/2`, where `delta_omega` is the difference of base Larmor frequencies.
- Builds a reference system with `t1_t2` only and explicitly augmented `r1`/`r2` rates, and compares against the SRSK generator for each spherical-tensor retention policy in `{'labframe','diagonal','kite','secular'}`.
- Verifies that the observed generator matches the once-only thermalised reference, that the laboratory-frame equilibrium state is stationary (`observed*rho_eq` close to zero), and that the trace functional is preserved (`unit'*observed` close to zero).

A second stage adds a spectator cavity mode (`C3`, frequency `1e9`, lifetime `0.5`) with independently specified mode dissipation. Mode damping/dephasing combinations `[0 3; 2 0; 2 3; 0 0]` are tested for couplings `100` and `0` Hz. For each combination the test checks:

- Mode dissipators preserve trace in all cases.
- Amplitude damping (`damp > 0`) is non-unital (does not annihilate the identity) and retains finite-temperature thermal occupation (differs from the zero-temperature dissipator).
- Pure dephasing and zero dissipation are unital (annihilate the identity).
- The combined spin-plus-mode generator equals spin thermalisation followed by a single mode dissipator, i.e. SRSK recursion excludes mode dissipation before outer thermalisation.
- The no-SRSK production path with analytically augmented spin rates reproduces the same spin-boson generator.

All comparisons use `test_close` with tolerances `1e-8` (absolute) and `1e-12` (relative), or `1e-10` for mode trace checks.

## Inputs and outputs

```matlab
result = test_srsk_thermal()
```

**Outputs:**

- `result` — regression check accumulator (from `new_test_result`) recording pass/fail status for rates, retention, equilibrium stationarity, trace preservation, and once-only thermalisation checks.

**Inputs:** none.

## References

- [Spinach library](https://spindynamics.org/) — the library this test belongs to.
- SRSK (scalar relaxation of the second kind) relaxation theory, as implemented in Spinach's `relaxation` module.
