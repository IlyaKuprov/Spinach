# examples/fundamentals/symmetry_1.m

- Signature: `symmetry_1()`

## Purpose

Liouvillian symmetrization for a radical pair with four equivalent nuclei under the S4 permutation group.

## Physical / mathematical content

The zero-field model contains two electron spins and four equivalent protons; the four protons are assigned to an S4 symmetry group. The basis uses the full set of symmetry sectors rather than restricting to A1g alone.

## Numerical / algorithmic content

Constructs the Hamiltonian superoperator, concatenates the symmetry-irrep projectors into a transformation matrix, and displays sparsity plots for the original Liouvillian and its symmetry-transformed form.

## Implementation structure

Defines the spin system and interactions, builds the spherical-tensor Liouville basis with S4 symmetry, makes the lab-frame assumption, then compares the two sparsity patterns.
