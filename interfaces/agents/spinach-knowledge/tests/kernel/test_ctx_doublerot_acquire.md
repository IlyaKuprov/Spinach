# tests/kernel/test_ctx_doublerot_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_doublerot_acquire.m`
- Signature: `result=test_ctx_doublerot_acquire()`
- Total lines: 71

## Purpose

Tests the double-rotor context with acquire(). Syntax: result=test_ctx_doublerot_acquire()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `doublerot()`, `acquire()`, `test_spin_system()`, `state()`, `test_close()`, `fid()`, `test_true()`, `all()`.
