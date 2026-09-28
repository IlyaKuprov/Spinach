# examples/fundamentals/clebsch_gordan/cg_test_large.m

- Signature: `cg_test_large()`

## Purpose

Checks Spinach Clebsch-Gordan coefficients against the arbitrary-precision reference values in the large test table, which are reported as Mathematica results.

## Physical and mathematical content

The test exercises Clebsch-Gordan coefficient evaluation for the cases stored in `cg_test_table_large.mat`; it is a numerical reference comparison rather than a spectrum simulation.

## Numerical and algorithmic content

For each row, the script computes a coefficient and compares it with the stored reference. An absolute difference below `2*eps` is reported as PASS; otherwise the script raises an error with the difference and fails.

## Implementation structure

The function loads `cg_test_table_large.mat`, loops over its rows, calls `clebsch_gordan` using the row's inputs, and checks the result against the reference in that row.
