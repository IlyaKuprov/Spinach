# examples/spin_chemistry/singlet_yield_3.m

- Signature: `singlet_yield_3()`

## Purpose

Figure 1 from the paper by Timmel, Till, Brocklehurst, McLauchlan and Hore: Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Figure 1 from the paper by Timmel, Till, Brocklehurst, McLauchlan
- and Hore:
- Calculation time: seconds
- Unit magnet (field sweep)
- Spin system
- Basis set
- Couplings
- Sequence parameters
- Spinach housekeeping
- Simulation
- Plotting
