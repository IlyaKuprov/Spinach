# tests/kernel/test_scalar_coupling_hamiltonian.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_scalar_coupling_hamiltonian.m`
- Signature: `result=test_scalar_coupling_hamiltonian()`
- Total lines: 46

## Purpose

Tests the two-spin scalar-coupling Hamiltonian. Syntax: result=test_scalar_coupling_hamiltonian()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Outputs

- result -regression test result with explanatory messages
- The test checks that an isotropic scalar coupling J produces the textbook
- Hamiltonian 2*pi*J*(Ix*Sx+Iy*Sy+Iz*Sz).

## Implementation structure

- Tests the two-spin scalar-coupling Hamiltonian. Syntax:
- result=test_scalar_coupling_hamiltonian()
- result -regression test result with explanatory messages
- The test checks that an isotropic scalar coupling J produces the textbook
- Hamiltonian 2*pi*J*(Ix*Sx+Iy*Sy+Iz*Sz).
- Announce the test target
- State the Hamiltonian target of the test
- Build a two-proton Hilbert-space spin system with a 10 Hz J coupling
- Build Spinach and textbook Hamiltonians
- Check the scalar-coupling Hamiltonian

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_spin_system()`, `hamiltonian()`, `assume()`, `operator()`, `test_close()`.
