# tests/kernel/test_dynamic_remaining_regularisation_suite.m

- Signature: `result=test_dynamic_remaining_regularisation_suite()`

## Purpose

Tests remaining regularisation and inverse-problem utilities. Syntax: result=test_dynamic_remaining_regularisation_suite()

## Physical / mathematical content

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
