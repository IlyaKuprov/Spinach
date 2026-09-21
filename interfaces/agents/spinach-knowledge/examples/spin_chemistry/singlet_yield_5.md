# examples/spin_chemistry/singlet_yield_5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_5.m`
- Signature: `singlet_yield_5()`
- Total lines: 47

## Purpose

Figure 5 from the paper by Timmel, Till, Brocklehurst, McLauchlan and Hore: Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Figure 5 from the paper by Timmel, Till, Brocklehurst, McLauchlan
- and Hore:
- Calculation time: seconds
- Unit magnet (field sweep)
- System specification
- Basis set
- Magnetic field and kinetics
- Spinach run
- Simulation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `liquid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
