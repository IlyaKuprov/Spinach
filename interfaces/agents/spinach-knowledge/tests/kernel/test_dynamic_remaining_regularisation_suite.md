# tests/kernel/test_dynamic_remaining_regularisation_suite.m

- Signature: `result=test_dynamic_remaining_regularisation_suite()`

## Purpose

Tests remaining regularisation and inverse-problem utilities. Syntax: result=test_dynamic_remaining_regularisation_suite()


## Numerical / algorithmic content

- Uses a synthetic L-curve with a known corner, scalar analytic references for positivity-constrained Tikhonov inversion, and an identity sensing matrix to test L1 sparse support and finite output metrics.
## Outputs

- `result` — regression test result with explanatory messages.
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
