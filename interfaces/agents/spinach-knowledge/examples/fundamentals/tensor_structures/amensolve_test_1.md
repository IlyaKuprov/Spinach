# examples/fundamentals/tensor_structures/amensolve_test_1.m

- MATLAB implementation: [examples/fundamentals/tensor_structures/amensolve_test_1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/tensor_structures/amensolve_test_1.m)

- Signature: `amensolve_test_1()`
- Source: [`examples/fundamentals/tensor_structures/amensolve_test_1.m`](../../../../../../examples/fundamentals/tensor_structures/amensolve_test_1.m)

## Purpose

Exercises `ttclass/amensolve` on buffered tensor-train linear systems, comparing the approximate solve with a constructed exact solution and dense references. It covers three Hermitian positive-definite cases, a repeatability check, a zero-enrichment finite-output regression, and a nonsymmetric smoke case.

## Tensor construction and mathematical scope

This is an algebraic tensor-train test, not a specified magnetic-resonance model: it defines no spins, Hamiltonian, or physical basis. For each positive-definite case, random rectangular core factors form buffered `B` and random vector factors form a known buffered tensor-train solution `x_exact`. The operator is constructed as `shrink(B' * B + diag_shift * unit_like(B' * B))`; the right-hand side is `shrink(A * x_exact)`. Dense `A`, `x_exact`, and `y` are materialised for references.

## Cases and solver settings

The script seeds the random generator with `rng(1)`. Mode sizes are algebraic dimensions, not physical units:

| Case | Mode sizes | Dense vector dimension | Operator terms | RHS terms | Diagonal shift | Solver tolerance |
|---|---:|---:|---:|---:|---:|---:|
| `small_exact` | [5 4 3 2] | 120 | 3 | 2 | 5e-2 | 1e-10 |
| `medium_balanced` | [10 10 10] | 1000 | 3 | 2 | 8e-2 | 1e-8 |
| `large_2000` | [20 10 10] | 2000 | 2 | 2 | 1e-1 | 2e-8 |

The per-case AMEn options set sweep limits (80, 120, or 140), initial guess ranks (2), enrichment ranks (4 or 6), rank caps (24 or 32), dense-local size limits (450–600), local iteration limits (150–260), and verbosity 0. The source also applies case-specific solution-error and residual limits derived from the requested tolerance; these are test acceptance criteria, not observed results.

## Use and checks

Run `amensolve_test_1` in the Spinach MATLAB environment. The test calls `amensolve(A_tt,y_tt,tol,opts)`, converts the returned train with `full`, and checks relative solution error and relative residual against dense references. It also asserts a single-train output, preserved mode sizes, finite values, and repeatability. The zero-enrichment regression requires finite error and residual values; the nonsymmetric smoke uses mode sizes [6 6 5], four operator and four RHS terms, seed 123, and tolerance 1e-8, and checks finite output and residual.

No physical units are assigned to dimensions, shift, or tolerance. Source assertions describe intended checks only; they do not establish that this test was run or passed.
