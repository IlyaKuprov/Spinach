# examples/spin_chemistry/singlet_yield_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_2.m`
- Signature: `singlet_yield_2()`
- Total lines: 55

## Purpose

Liquid state magnetic field effect simulation on a radical pair with six equivalent nuclei using exponential recombi- nation kinetics model. Full S6 symmatry is used. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Liquid state magnetic field effect simulation on a radical
- pair with six equivalent nuclei using exponential recombi-
- nation kinetics model. Full S6 symmatry is used.
- Calculation time: seconds
- Unit magnet (field sweep)
- System specification
- Basis set
- Fields and kinetics parameters
- Spinach housekeeping
- Simulation
- Plot the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `mt2hz()`, `create()`, `basis()`, `liquid()`, `kfigure()`, `kxlabel()`, `kylabel()`.
