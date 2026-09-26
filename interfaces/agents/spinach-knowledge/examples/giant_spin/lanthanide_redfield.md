# examples/giant_spin/lanthanide_redfield.m

- Signature: `lanthanide_redfield()`

## Purpose

Relaxation rate of Gd(III) as a function of zero-field splitting, computed using Redfield theory. The correla- tion time (1 fs) refers to vibrational dynamics of the ligand cage. Calculation time: minutes

## Physical / mathematical content

- Giant-spin examples. The effective model treats lanthanides or high-spin centres using crystal-field / Stevens-operator Hamiltonians, Zeeman splitting, and magnetisation dynamics.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Relaxation rate of Gd(III) as a function of zero-field
- splitting, computed using Redfield theory. The correla-
- tion time (1 fs) refers to vibrational dynamics of the
- ligand cage.
- Calculation time: minutes
- Spin system properties
- Magnet field
- Relaxation parameters
- Basis set
- Linearly spaced B20 in cm^-1
- Loop over B20 values
- Giant spin Hamiltonian parameters
