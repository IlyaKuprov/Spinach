# tests/kernel/test_ctx_doublerot_acquire.m

- Signature: `result=test_ctx_doublerot_acquire()`

## Purpose

Tests the double-rotor context with acquire(). Syntax: result=test_ctx_doublerot_acquire()

## Physical / mathematical content

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin double-rotation calculation
- through doublerot() and checks the returned time-domain trace for basic
- physical and dimensional invariants.

## Implementation structure

- Tests the double-rotor context with acquire(). Syntax:
- result=test_ctx_doublerot_acquire()
- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin double-rotation calculation
- through doublerot() and checks the returned time-domain trace for basic
- physical and dimensional invariants.
- Announce the test target
- State the double-rotor target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a tiny double-rotation acquisition
- Run the production double-rotor context
- Check the number of acquired points
