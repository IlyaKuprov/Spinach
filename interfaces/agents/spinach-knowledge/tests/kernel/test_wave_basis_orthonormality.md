# tests/kernel/test_wave_basis_orthonormality.m

- Signature: `result=test_wave_basis_orthonormality()`

## Purpose

Tests that the sine, cosine, and Legendre waveform bases returned by Spinach have orthonormal columns, as required by pulse optimisation.

## Physical / mathematical content

Orthonormal basis columns give independent waveform coefficients.

## Numerical / algorithmic content

For each basis, the test constructs `B=wave_basis(basis_type,5,32)` and checks that its Gram matrix `B'*B` equals `eye(5)` within absolute and relative tolerances of `1e-12`.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

The test announces its target, creates a regression test result, then checks the `sine_waves`, `cosine_waves`, and `legendre` bases in a loop.