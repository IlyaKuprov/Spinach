# examples/fundamentals/convention_tests/rotations_2.m

- Signature: `rotations_2()`

## Purpose

A rotations test comparing the Hamiltonians for a manually rotated (at the interaction specification level) spin system with the Hamiltonian that has been rotated using Spinach operator rotation functionality.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- A rotations test comparing the Hamiltonians for a manually rotated (at
- the interaction specification level) spin system with the Hamiltonian
- that has been rotated using Spinach operator rotation functionality.
- Generate random matrices
- % Kernel level rotation
- Magnet field
- Basis set
- A pair of spins at a distance, A
- Spinach housekeeping, A
- Hamiltonian, A
- % Input level rotation
- A pair of spins at a distance, B
