# examples/spin_chemistry/singlet_yield_anisotropy_1.m

- Signature: `singlet_yield_anisotropy_1()`

## Purpose

Singlet yield anisotropy calculation for a radical pair using exponential recombination kinetics model. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Singlet yield anisotropy calculation for a radical pair
- using exponential recombination kinetics model.
- Calculation time: seconds
- Unit magnet (field sweep)
- System specification
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Preprocessing
- Plotting
