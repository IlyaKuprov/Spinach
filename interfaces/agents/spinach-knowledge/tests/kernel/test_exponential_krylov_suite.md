# tests/kernel/test_exponential_krylov_suite.m

- Signature: `result=test_exponential_krylov_suite()`

## Purpose

Regression tests for Krylov and exponential-integral routines.

## Tests

- Runs `arnoldi` on a 3-by-3 diagonal matrix and checks its basis and recurrence identities.
- Checks `expdrop` and `expmint` cases, including a diagonal input and zero time.
- Checks `expmint2` on a nested-integral case.

The suite does not test dynamical propagation or performance on large matrices.

## Outputs

- result -regression test result with explanatory messages
- The test checks Arnoldi basis identities, Chebyshev coefficients,
- exponential drop boundary values, and Van Loan exponential-integral
- helpers against small closed-form references.
