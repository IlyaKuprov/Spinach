# tests/kernel/test_ctx_gridfree_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_gridfree_acquire.m`
- Signature: `result=test_ctx_gridfree_acquire()`
- Total lines: 66

## Purpose

Tests the grid-free Fokker-Planck context with acquire(). Syntax: result=test_ctx_gridfree_acquire()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin MAS calculation through
- gridfree() and checks the returned time-domain trace for basic physical
- and dimensional invariants.

## Implementation structure

- Tests the grid-free Fokker-Planck context with acquire(). Syntax:
- result=test_ctx_gridfree_acquire()
- result -regression test result with explanatory messages
- The test runs a tiny anisotropic one-spin MAS calculation through
- gridfree() and checks the returned time-domain trace for basic physical
- and dimensional invariants.
- Announce the test target
- State the grid-free target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a tiny grid-free acquisition
- Run the production grid-free context
- Check the number of acquired points

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `gridfree()`, `acquire()`, `test_spin_system()`, `state()`, `test_close()`, `fid()`, `test_true()`, `all()`.
