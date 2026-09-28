# examples/fundamentals/tensor_structures/amensum_test_1.m

- Signature: `amensum_test_1()`

## Purpose

Tests `ttclass/amensum` by summing buffered rank-one tensor trains and comparing the result with their exact dense sum. It checks relative Frobenius error, output structure, physical dimensions, and finiteness. The source notes that AMEn summation is approximate and its motivating paper focuses on enrichment-assisted updates, so the strict accuracy checks apply to enriched runs; zero enrichment is retained as a finite-output regression smoke test.

## Physical / mathematical content

- The inputs are sums of Kronecker products with generated coefficients, including signed coefficients. The dense reference is assembled term by term from the same factors.
- No spin system or physical time evolution is constructed.

## Numerical / algorithmic content

- The accuracy suite covers `small_exact`, `medium_balanced`, `large_2000` (dimension 2000), and `signed_coeffs` cases, with case-specific term counts, tolerances, and AMEn options.
- Relative Frobenius error is compared with `max(50*tol,1e-12)`. A repeated medium case resets the random-number generator before each solve and compares the resulting errors for reproducibility.
- With `enrichment_rank=0`, the test checks that the relative error is finite but deliberately does not require the enriched-run accuracy threshold.

## Implementation structure

- Generate normalised complex Kronecker factors and signed coefficients, store the rank-one terms in a `ttclass`, and assemble their dense sum.
- Call `amensum` for each specified case and check its dense-reference error, single-train output, mode sizes, and finiteness.
- Run the reproducibility and zero-enrichment smoke checks, then report success.
