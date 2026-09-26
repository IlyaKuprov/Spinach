# tests/kernel/test_dynamic_trajectory_frontends.m

- Signature: `result=test_dynamic_trajectory_frontends()`

## Purpose

Tests trajectory-analysis dynamic front-end kernels. Syntax: result=test_dynamic_trajectory_frontends()

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test exercises trajan() plotting branches and trajsimil() scoring
- branches on a compact two-spin spherical-tensor trajectory.

## Implementation structure

- Tests trajectory-analysis dynamic front-end kernels. Syntax:
- result=test_dynamic_trajectory_frontends()
- result -regression test result with explanatory messages
- The test exercises trajan() plotting branches and trajsimil() scoring
- branches on a compact two-spin spherical-tensor trajectory.
- Announce the test target
- State the dynamic trajectory target of the test
- Force invisible figures during plotting checks
- Build the trajectory used by plotting and similarity checks
- Check all trajan() property branches
- Check all trajsimil() scoring families
- Check correlation-order analysis with an explicit time axis
