# tests/kernel/test_ctx_singlerot_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_singlerot_acquire.m`
- Signature: `result=test_ctx_singlerot_acquire()`
- Total lines: 68

## Purpose

Tests the single-rotor context with acquire(). Syntax: result=test_ctx_singlerot_acquire()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `singlerot()`, `acquire()`, `test_spin_system()`, `state()`, `test_close()`, `fid()`, `test_true()`, `all()`.
