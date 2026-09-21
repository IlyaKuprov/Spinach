# tests/kernel/test_ctx_floquet_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_floquet_acquire.m`
- Signature: `result=test_ctx_floquet_acquire()`
- Total lines: 68

## Purpose

Tests the Floquet context with acquire(). Syntax: result=test_ctx_floquet_acquire()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `floquet()`, `acquire()`, `test_spin_system()`, `state()`, `test_close()`, `fid()`, `test_true()`, `all()`.
