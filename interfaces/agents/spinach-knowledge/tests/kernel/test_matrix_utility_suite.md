# tests/kernel/test_matrix_utility_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_matrix_utility_suite.m`
- Signature: `result=test_matrix_utility_suite()`
- Total lines: 118

## Purpose

Tests small matrix utility functions. Syntax: result=test_matrix_utility_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `acomm()`, `cheap_norm()`, `polyadic()`, `onCleanup()`, `rng()`, `test_true()`, `clear()`, `iseye()`, `speye()`, `istraceless()`, `krondelta()`, `sp_block_diag()`, `repcols()`, `reprows()`.
