# tests/kernel/test_hilbert_state.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_hilbert_state.m`
- Signature: `result=test_hilbert_state()`
- Total lines: 44

## Purpose

Tests Hilbert-space state generation. Syntax: result=test_hilbert_state()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks that state() returns the expected density matrices in
- Hilbert space for a one-spin system.

## Implementation structure

- Tests Hilbert-space state generation. Syntax:
- result=test_hilbert_state()
- result -regression test result with explanatory messages
- The test checks that state() returns the expected density matrices in
- Hilbert space for a one-spin system.
- Announce the test target
- State the physical target of the test
- Build a one-proton Hilbert-space spin system
- Textbook spin-half reference matrices
- Check density matrices generated from state labels

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `state()`, `test_spin_system()`, `pauli()`, `test_close()`.
