# tests/kernel/test_matrix_utility_suite.m

- Signature: `result=test_matrix_utility_suite()`

## Purpose

Tests small matrix utility functions. Syntax: result=test_matrix_utility_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks low-level matrix helpers that Spinach uses throughout the
- kernel for algebra, sparsity, block assembly, and indexing operations.

## Implementation structure

- Tests small matrix utility functions. Syntax:
- result=test_matrix_utility_suite()
- result -regression test result with explanatory messages
- The test checks low-level matrix helpers that Spinach uses throughout the
- kernel for algebra, sparsity, block assembly, and indexing operations.
- Announce the test target
- State the utility target of the test
- Define small test matrices
- Check anticommutator and cheap norm
- Check polyadic norm estimation paths
- Check rectangular polyadic sign-history dimensions
- Check the Higham-Tisseur zero-sign convention
