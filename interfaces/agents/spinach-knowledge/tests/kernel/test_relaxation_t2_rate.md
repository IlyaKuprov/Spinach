# tests/kernel/test_relaxation_t2_rate.m

- Signature: `result=test_relaxation_t2_rate()`

## Purpose

Tests phenomenological T2 relaxation rate. Syntax: result=test_relaxation_t2_rate()

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks that the t1_t2 relaxation model gives transverse L+ order
- the negative generator eigenvalue corresponding to the specified R2 rate.

## Implementation structure

- Tests phenomenological T2 relaxation rate. Syntax:
- result=test_relaxation_t2_rate()
- result -regression test result with explanatory messages
- The test checks that the t1_t2 relaxation model gives transverse L+ order
- the negative generator eigenvalue corresponding to the specified R2 rate.
- Announce the test target
- State the relaxation target of the test
- Build a one-spin relaxation system
- Build relaxation superoperator and transverse state
- Check that L+ is an eigenstate with the negative R2 generator eigenvalue
