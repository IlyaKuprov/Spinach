# examples/fundamentals/convention_tests/blinv_test.m

- Signature: `blinv_test()`

## Purpose

Checks rank-2 Blicharski invariants against the corresponding spherical-tensor coefficient products, for individual tensors and a cross-product of two tensors.

## Checks

The test draws two random 3×3 tensors, converts each to spherical coefficients with `mat2sphten`, and obtains the self-invariants with `blinv` and the cross-invariant with `blprod`. It compares the rank-2 coefficient expressions with (2/3) times each self-invariant and with (2/3) times the cross-invariant. Every residual must be no larger than 10 eps; otherwise the corresponding test fails.
