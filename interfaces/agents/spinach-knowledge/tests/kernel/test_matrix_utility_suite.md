# tests/kernel/test_matrix_utility_suite.m

- Signature: `result=test_matrix_utility_suite()`

## Purpose

Tests small matrix utility functions used throughout the Spinach kernel for algebra, sparsity, block assembly, and indexing operations.

## Numerical / algorithmic content

- Checks that `acomm(A,B)` equals `A*B+B*A` and that `cheap_norm(A)` returns the matrix one-norm for a CPU matrix.
- Tests `cheap_norm` on non-negative, tall, wide, real, and complex polyadic matrices, including block estimation, rectangular sign-history dimensions, conjugate-transpose and unit-phase handling, and the Higham–Tisseur Algorithm 2.4 zero-sign convention.
- Checks that the block polyadic estimate is positive and does not exceed the exact one-norm for a test matrix. The test saves and restores the random-number-generator state.
- Checks `isreal` dispatch for polyadic matrices; `iseye`, `istraceless`, and `krondelta` predicates; sparse block-diagonal assembly; column and row replication; and removal of sub-tolerance sparse elements by `clean_up`.

## Outputs

- `result` — regression test result with explanatory messages.
