# tests/kernel/test_ctx_floquet_acquire.m

- Signature: `result=test_ctx_floquet_acquire()`

## Purpose

Tests the Floquet context with acquire(). Syntax: result=test_ctx_floquet_acquire()

## Physical / mathematical content

- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin MAS calculation through
- floquet() and checks the returned time-domain trace for basic physical
- and dimensional invariants.

## Implementation structure

- Tests the Floquet context with acquire(). Syntax:
- result=test_ctx_floquet_acquire()
- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin MAS calculation through
- floquet() and checks the returned time-domain trace for basic physical
- and dimensional invariants.
- Announce the test target
- State the Floquet-context target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a tiny Floquet acquisition
- Run the production Floquet context
- Check the number of acquired points
