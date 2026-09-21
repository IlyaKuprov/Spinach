# tests/kernel/test_zeeman_hamiltonian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_zeeman_hamiltonian.m`
- Signature: `result=test_zeeman_hamiltonian()`
- Total lines: 42

## Purpose

Tests the one-spin Zeeman Hamiltonian. Syntax: result=test_zeeman_hamiltonian()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks Spinach's NMR convention for a positive chemical shift:
- the rotating-frame Hamiltonian contribution is -2*pi*nu*Lz.

## Implementation structure

- Tests the one-spin Zeeman Hamiltonian. Syntax:
- result=test_zeeman_hamiltonian()
- result -regression test result with explanatory messages
- The test checks Spinach's NMR convention for a positive chemical shift:
- the rotating-frame Hamiltonian contribution is -2*pi*nu*Lz.
- Announce the test target
- State the Hamiltonian target of the test
- Build a one-proton Hilbert-space spin system with a 1 ppm shift
- Build Spinach and reference Hamiltonians
- Check the physical frequency and sign convention

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `hamiltonian()`, `assume()`, `ppm2hz()`, `operator()`, `test_close()`.
