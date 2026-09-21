# tests/kernel/test_dynamic_remaining_regularisation_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_remaining_regularisation_suite.m`
- Signature: `result=test_dynamic_remaining_regularisation_suite()`
- Total lines: 62

## Purpose

Tests remaining regularisation and inverse-problem utilities. Syntax: result=test_dynamic_remaining_regularisation_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Outputs

- result -regression test result with explanatory messages
- The test checks L-curve corner detection, positivity-constrained
- Tikhonov inversion, and L1 sparsity targeting on compact analytical
- inverse problems.

## Implementation structure

- Tests remaining regularisation and inverse-problem utilities. Syntax:
- result=test_dynamic_remaining_regularisation_suite()
- result -regression test result with explanatory messages
- The test checks L-curve corner detection, positivity-constrained
- Tikhonov inversion, and L1 sparsity targeting on compact analytical
- inverse problems.
- Announce the test target
- State the utility target of the test
- Check L-curve analysis on a synthetic corner near lambda equals one
- Check positivity-constrained Tikhonov inversion against a scalar analytic solution
- Check L1 sparsity targeting on an identity sensing matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `lcurve()`, `test_true()`, `tikhonov()`, `test_close()`, `minimising()`, `tikhol1n()`, `nnz()`, `x_l1()`.
