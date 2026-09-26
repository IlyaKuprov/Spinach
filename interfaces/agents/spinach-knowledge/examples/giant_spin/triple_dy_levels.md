# examples/giant_spin/triple_dy_levels.m

- Signature: `triple_dy_levels()`

## Purpose

Eight lowest energy levels as a function of the applied magnetic fi- eld in a triple Dy triangular complex -see Figure 12 in Ligand field parameters and g-tensor for the J=15/2 ground term were computed using the SINGLE_ANISO routine in MOLCAS. Calculation time: hours

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Eight lowest energy levels as a function of the applied magnetic fi-
- eld in a triple Dy triangular complex -see Figure 12 in
- Ligand field parameters and g-tensor for the J=15/2 ground term were
- computed using the SINGLE_ANISO routine in MOLCAS.
- Calculation time: hours
- Three J=15/2 dysprosium atoms
- g-tensor eigenvalues
- Spin-orbit corrections
- to the DD couplings
- g-tensor eigenvectors
- g-tensor matrix
- Triangle arrangement
