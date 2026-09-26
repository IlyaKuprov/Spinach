# examples/giant_spin/triple_tb_eqmag_temp.m

- Signature: `triple_tb_eqmag_temp()`

## Purpose

Simulation of the temperature dependence of the magnetisation of a triple-Tb triangular complex -see Figure S27 and S28 in the Supple- mentary Information of the following paper Ligand field parameters and g-tensor for the J=6 ground term were computed using the SINGLE_ANISO routine in MOLCAS. Calculation time: hours

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simulation of the temperature dependence of the magnetisation of a
- triple-Tb triangular complex -see Figure S27 and S28 in the Supple-
- mentary Information of the following paper
- Ligand field parameters and g-tensor for the J=6 ground term were
- computed using the SINGLE_ANISO routine in MOLCAS.
- Calculation time: hours
- Three J=6 terbium atoms
- g-tensor eigenvalues
- Spin-orbit corrections
- to the DD couplings
- g-tensor eigenvalues (rows)
- g-tensor matrices
