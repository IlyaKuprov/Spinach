# tests/kernel/test_dynamic_cubic_mex_suite.m

- Signature: `result=test_dynamic_cubic_mex_suite()`

## Purpose

Tests the cubic-polynomial MEX helper used by `eigenfields()`.

## Physical / mathematical content

Checks roots in the unit interval, including endpoint and repeated roots. The identically zero polynomial is ignored by the eigenfields root filter.

## Numerical / algorithmic content

The test uses the production root tolerance `sqrt(eps)` and compares `cubic_roots` results with explicit references for three roots, a triple root, a double root with an endpoint root, quadratic and linear degeneracies, and constant and zero polynomials. It checks extreme coefficient scales (`1e300` and `1e-300`) and a derivative-root use case against its two analytical turning points. It also performs 200 seeded randomized cubic round trips: well-separated roots are converted to coefficients with `poly`, recovered with `cubic_roots`, and checked against the original roots. The previous RNG state is restored.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

The function announces the cubic-polynomial root MEX test target, creates a `new_test_result` for `kernel/dynamic_cubic_mex_suite`, runs the reference and randomized checks, and records their outcomes with `test_close` and `test_true`.

Source contact: ilya.kuprov@weizmann.ac.il