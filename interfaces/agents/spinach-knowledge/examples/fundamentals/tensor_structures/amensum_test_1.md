# examples/fundamentals/tensor_structures/amensum_test_1.m

- MATLAB implementation: [examples/fundamentals/tensor_structures/amensum_test_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/tensor_structures/amensum_test_1.m)

- Signature: `amensum_test_1()`
- Source: [`examples/fundamentals/tensor_structures/amensum_test_1.m`](../../../../../../examples/fundamentals/tensor_structures/amensum_test_1.m)

## Purpose

Tests `ttclass/amensum` by summing buffered rank-one tensor trains and comparing the result with the exact dense sum. The source explicitly treats AMEn summation as approximate: strict accuracy checks are for enrichment-assisted runs, while zero enrichment is a finite-output smoke test. It does not identify a particular paper or provide a DOI.

## Tensor construction and mathematical scope

This is a numerical tensor algebra test, not a magnetic-resonance spin model: no spin system, Hamiltonian, or physical basis is specified. Each buffered term is represented by a Kronecker product of one random matrix factor per core, weighted by a signed coefficient `(-1)^n * (0.5 + rand())`. These factors and coefficients are stored as a `ttclass`; the reference is assembled by summing the corresponding dense Kronecker products. Relative error is measured in the Frobenius norm.

## Cases and settings

The main cases are seeded from `rng(1)`; the builder resets the random generator per case. Dimensions are tensor mode sizes and carry no physical units.

| Case | Mode sizes | Terms | Tolerance | AMEn options (sweeps, initial rank, enrichment rank) |
|---|---:|---:|---:|---|
| `small_exact` | [5 4 3 2] | 6 | 1e-12 | 80, 2, 4 |
| `medium_balanced` | [10 10 10] | 24 | 1e-10 | 80, 2, 4 |
| `large_2000` | [20 10 10] | 32 | 1e-8 | 120, 2, 6 |
| `signed_coeffs` | [5 5 5 4] | 18 | 1e-9 | 120, 3, 5 |

All cases set verbosity to 0. The separate zero-enrichment smoke uses mode sizes [8 10 10], 20 terms, seed 77, tolerance 1e-8, 160 sweeps, initial rank 2, and enrichment rank 0. No physical units apply to these parameters.

## Use and checks

Run `amensum_test_1` in the Spinach MATLAB environment. The test calls `amensum(x_ref,tol,opts)`, materialises the resulting tensor train with `full`, and checks relative Frobenius error against the dense sum, a single-train output, matching mode sizes, finiteness, and reproducibility. The source sets the enriched-run error limit to `max(50*tol,1e-12)`; the no-enrichment case only checks that its error is finite. These source assertions are not evidence of a test run or pass.
