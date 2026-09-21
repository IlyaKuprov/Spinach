# examples/spin_chemistry/singlet_yield_anisotropy_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_anisotropy_2.m`
- Signature: `singlet_yield_anisotropy_2()`
- Total lines: 73

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
- Isotopes
- Basis set
- Rotation matrices
- Interaction eigenvalues
- Coupling tensors
- Zeeman interactions
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `powder()`, `cell2mat()`, `get_hull()`, `kfigure()`, `trisurf()`.
