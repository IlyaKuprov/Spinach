# tests/kernel/test_ctx_singlerot_acquire.m

- Signature: `result=test_ctx_singlerot_acquire()`

## Purpose

Tests the single-rotor context with acquire(). Syntax: result=test_ctx_singlerot_acquire()

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin MAS calculation through
- singlerot() and checks the returned time-domain trace for basic
- physical and dimensional invariants.

## Implementation structure

- Tests the single-rotor context with acquire(). Syntax:
- result=test_ctx_singlerot_acquire()
- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin MAS calculation through
- singlerot() and checks the returned time-domain trace for basic
- physical and dimensional invariants.
- Announce the test target
- State the single-rotor target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a tiny MAS acquisition
- Run the production single-rotor context
- Check the number of acquired points
