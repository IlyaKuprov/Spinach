# examples/spin_chemistry/singlet_yield_anisotropy_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/singlet_yield_anisotropy_3.m`
- Signature: `singlet_yield_anisotropy_3()`
- Total lines: 72

## Purpose

Singlet yield anisotropy calculation for a model radical pair reaction, Haberkorn recombination model. Run time: minutes on NVidia Titan V card, hours on CPU.

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Singlet yield anisotropy calculation for a model radical pair
- reaction, Haberkorn recombination model.
- Run time: minutes on NVidia Titan V card, hours on CPU.
- Earth field
- Isotopes
- Basis set
- Hyperfine coupling tensors
- Zeeman interactions
- Kinetics parameters
- Sequence parameters
- Enable GPU arithmetic
- sys.enable={'gpu'};

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mt2hz()`, `create()`, `basis()`, `powder()`, `kfigure()`, `cell2mat()`, `kxlabel()`, `kylabel()`.
