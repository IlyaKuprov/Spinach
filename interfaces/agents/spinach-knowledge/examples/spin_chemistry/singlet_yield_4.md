# examples/spin_chemistry/singlet_yield_4.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_4.m`
- Signature: `singlet_yield_4()`
- Total lines: 54

## Purpose

Figure 3 from the paper by Till, Timmel, Brocklehurst and Hore: Note: the original paper only uses electron Zeeman operators for the field sweep, and therefore misses the effects associa- ted with the rise in the nuclear Zeeman interaction on the high field side of the resulting plot. Calculation time: seconds

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Figure 3 from the paper by Till, Timmel, Brocklehurst and Hore:
- Note: the original paper only uses electron Zeeman operators for
- the field sweep, and therefore misses the effects associa-
- ted with the rise in the nuclear Zeeman interaction on the
- high field side of the resulting plot.
- Calculation time: seconds
- Unit magnet (field sweep)
- Spin system
- Basis set
- Couplings
- Sequence parameters
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gauss2mhz()`, `create()`, `basis()`, `liquid()`, `kfigure()`, `kylabel()`, `kxlabel()`.
