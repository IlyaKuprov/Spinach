# tests/kernel/test_ctx_liquid_acquire.m

- Signature: `result=test_ctx_liquid_acquire()`

## Purpose

Tests the liquid context against the direct acquire path. Syntax: result=test_ctx_liquid_acquire()

## Physical / mathematical content

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
