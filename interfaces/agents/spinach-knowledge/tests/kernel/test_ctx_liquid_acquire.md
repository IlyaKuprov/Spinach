# tests/kernel/test_ctx_liquid_acquire.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_ctx_liquid_acquire.m`
- Signature: `result=test_ctx_liquid_acquire()`
- Total lines: 62

## Purpose

Tests the liquid context against the direct acquire path. Syntax: result=test_ctx_liquid_acquire()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test runs a one-spin offset FID through liquid() and compares it
- with the same Hamiltonian, relaxation, and kinetics objects passed
- directly to acquire().

## Implementation structure

- Tests the liquid context against the direct acquire path. Syntax:
- result=test_ctx_liquid_acquire()
- result -regression test result with explanatory messages
- The test runs a one-spin offset FID through liquid() and compares it
- with the same Hamiltonian, relaxation, and kinetics objects passed
- directly to acquire().
- Announce the test target
- State the liquid-context target of the test
- Build a one-spin Liouville-space system
- Set up a short offset acquisition
- Run the production liquid context
- Build the same direct acquire input path

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `liquid()`, `acquire()`, `test_spin_system()`, `state()`, `assume()`, `hamiltonian()`, `relaxation()`, `kinetics()`, `frqoffset()`, `test_close()`.
