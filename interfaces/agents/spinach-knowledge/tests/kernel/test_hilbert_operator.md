# tests/kernel/test_hilbert_operator.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_hilbert_operator.m`
- Signature: `result=test_hilbert_operator()`
- Total lines: 48

## Purpose

Tests Hilbert-space operator generation. Syntax: result=test_hilbert_operator()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks that operator() builds the correct one-spin Hilbert-space
- angular momentum matrices from human-readable labels.

## Implementation structure

- Tests Hilbert-space operator generation. Syntax:
- result=test_hilbert_operator()
- result -regression test result with explanatory messages
- The test checks that operator() builds the correct one-spin Hilbert-space
- angular momentum matrices from human-readable labels.
- Announce the test target
- State the physical target of the test
- Build a one-proton Hilbert-space spin system
- Textbook spin-half reference matrices
- Check label-to-matrix mapping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `operator()`, `test_spin_system()`, `pauli()`, `test_close()`.
