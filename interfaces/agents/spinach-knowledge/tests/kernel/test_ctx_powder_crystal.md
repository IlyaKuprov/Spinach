# tests/kernel/test_ctx_powder_crystal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_powder_crystal.m`
- Signature: `result=test_ctx_powder_crystal()`
- Total lines: 58

## Purpose

Tests static powder and crystal contexts at one orientation. Syntax: result=test_ctx_powder_crystal()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test uses the single_crystal grid so that powder() and crystal()
- represent the same Euler orientation of the same anisotropic one-spin
- Hamiltonian.

## Implementation structure

- Tests static powder and crystal contexts at one orientation. Syntax:
- result=test_ctx_powder_crystal()
- result -regression test result with explanatory messages
- The test uses the single_crystal grid so that powder() and crystal()
- represent the same Euler orientation of the same anisotropic one-spin
- Hamiltonian.
- Announce the test target
- State the static-context target of the test
- Build a one-spin anisotropic Liouville-space system
- Set up a short static acquisition
- Run both production static contexts
- Check that one-point powder averaging is a crystal calculation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `powder()`, `crystal()`, `test_spin_system()`, `state()`, `test_close()`.
