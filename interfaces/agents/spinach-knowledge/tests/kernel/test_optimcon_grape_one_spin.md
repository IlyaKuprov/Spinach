# tests/kernel/test_optimcon_grape_one_spin.m

- Signature: `result=test_optimcon_grape_one_spin()`

## Purpose

Checks a minimal one-spin Hilbert-space optimal-control setup and the fidelity and gradient returned by `grape_hilb()`.

## Test

- Creates a one-spin `zeeman-hilb` system with a zero drift, one `S.y` control operator, `S.x` initial operator, `S.z` target operator, and two pulse intervals of `0.02` and `0.03`. Configures `optimcon()` with `lbfgs`, no penalties, and zero optimisation iterations.
- Verifies that `optimcon()` registers one control, preserves the timing grid, and assigns one waveform value per interval. It also checks that omitting the method defaults to `lbfgs` when a distortion is supplied, while a distortion with the exact-Hessian `newton` method is rejected.
- For waveform `[7 11]`, compares the `grape_hilb()` fidelity with independent Hilbert-space propagation using `expm()` to `1e-13`, and its gradient with centred finite differences using a `1e-6` step to `1e-8`.
- Checks that the forward trajectory is not stored when plotting is disabled.

## Output

- `result` — regression test result with explanatory messages.