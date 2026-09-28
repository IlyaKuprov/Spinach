# examples/fundamentals/tensor_structures/amensolve_test_1.m

- Signature: `amensolve_test_1()`

## Purpose

Tests `ttclass/amensolve` against dense references. The main cases construct positive-definite tensor-train systems from `B'*B` plus a diagonal shift, form the right-hand side from a known tensor-train solution, and compare AMEn's answer with the exact dense solution and dense residual. The suite also checks repeatability, a zero-enrichment regression, and finite output for a nonsymmetric smoke case.

## Physical / mathematical content

- This is a numerical tensor-train linear-solver test, not a spin-dynamics example.
- The positive-definite cases use a known solution `x_exact`, with `A=B'*B+diag_shift*I` and `y=A*x_exact`; the dense solution and residual provide independent checks of the computed result.
- The nonsymmetric case is only a finite-output smoke test; the source does not impose the positive-definite accuracy contract on it.

## Numerical / algorithmic content

- Cases span small and medium systems and a dense-reference system of dimension 2000. Relative solution error and relative residual are checked against case-dependent tolerances derived from each case's `tol`.
- A medium case is run twice with the random-number generator reset to test reproducibility. The zero-enrichment case checks that output and dense-reference error/residual remain finite; it does not apply the enriched-run accuracy threshold.
- The test also verifies that the result is one tensor train, preserves the expected mode sizes, and contains only finite values.

## Implementation structure

- Build test-case specifications with dimensions, term counts, diagonal shifts, tolerances, and AMEn options.
- Construct tensor-train operators and right-hand sides, then materialise dense references for comparison.
- Run `amensolve`, check dense solution error, residual, output structure, dimensions, and finiteness, and report each case.
- Run the reproducibility, zero-enrichment, and nonsymmetric smoke checks, then print the final success message.
