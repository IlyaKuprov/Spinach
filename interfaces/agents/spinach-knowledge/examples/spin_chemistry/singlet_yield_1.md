# examples/spin_chemistry/singlet_yield_1.m

- Signature: `singlet_yield_1()`

## Purpose

Liquid state magnetic field effect simulation on a radical pair with four nuclei using exponential recombination ki- netics model. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Liquid state magnetic field effect simulation on a radical
- pair with four nuclei using exponential recombination ki-
- netics model.
- Calculation time: seconds
- Unit magnet (field sweep)
- System specification
- Basis set
- Fields and kinetics parameters
- Disable ZTE
- Spinach housekeeping
- Simulation
- Plot the answer
