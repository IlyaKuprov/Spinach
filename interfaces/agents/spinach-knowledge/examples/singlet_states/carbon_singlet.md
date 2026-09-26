# examples/singlet_states/carbon_singlet.m

- Signature: `carbon_singlet()`

## Purpose

Singlet relaxation rate for the two triple bond carbons in cis-dimethylbut-2-ynedioate. Magnetic parameters com- puted with DFT. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Singlet relaxation rate for the two triple bond carbons
- in cis-dimethylbut-2-ynedioate. Magnetic parameters com-
- puted with DFT.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis set
- Relaxation superoperator accuracy
- Spinach housekeeping
- Relaxation superoperator
- Action on longitudinal magnetization
- Action on a singlet state
