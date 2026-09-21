# examples/giant_spin/triple_tb_eqmag_field.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/giant_spin/triple_tb_eqmag_field.m`
- Signature: `triple_tb_eqmag_field()`
- Total lines: 222

## Purpose

Simulation of the field dependence of the magnetisation of a triple- Tb triangular complex -see Figure S27 and S28 in the Supplementary Information of the following paper Ligand field parameters and g-tensor for the J=6 ground term were computed using the SINGLE_ANISO routine in MOLCAS. Calculation time: hours

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Simulation of the field dependence of the magnetisation of a triple-
- Tb triangular complex -see Figure S27 and S28 in the Supplementary
- Information of the following paper
- Ligand field parameters and g-tensor for the J=6 ground term were
- computed using the SINGLE_ANISO routine in MOLCAS.
- Calculation time: hours
- Three J=6 terbium atoms
- g-tensor eigenvalues
- Spin-orbit corrections
- to the DD couplings
- g-tensor eigenvalues (rows)
- g-tensor matrices

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `icm2hz()`, `dcm2euler()`, `wigner()`, `stev2sph()`, `nan()`, `kfigure()`, `create()`, `basis()`, `eqmag()`, `mag()`, `kxlabel()`, `kylabel()`, `load()`.
