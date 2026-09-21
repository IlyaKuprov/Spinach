# tests/kernel/test_relaxation_t2_rate.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_relaxation_t2_rate.m`
- Signature: `result=test_relaxation_t2_rate()`
- Total lines: 48

## Purpose

Tests phenomenological T2 relaxation rate. Syntax: result=test_relaxation_t2_rate()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `relaxation()`, `state()`, `test_close()`.
