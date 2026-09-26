# tests/kernel/test_scalar_coupling_hamiltonian.m

- Signature: `result=test_scalar_coupling_hamiltonian()`

## Purpose

Tests the two-spin scalar-coupling Hamiltonian. Syntax: result=test_scalar_coupling_hamiltonian()

## Physical / mathematical content

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
