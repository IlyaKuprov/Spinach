# examples/relaxation_theory/hfc_relaxation_1.m

- Signature: `hfc_relaxation_1()`

## Purpose

Computes and prints the full Redfield superoperator for an electron- nucleus system with an anisotropic hyperfine coupling in liquid state. Hyperfine coupling is computed from the Cartesian coordinates using the point dipole approximation. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Computes and prints the full Redfield superoperator for an electron-
- nucleus system with an anisotropic hyperfine coupling in liquid state.
- Hyperfine coupling is computed from the Cartesian coordinates using
- the point dipole approximation.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Relaxation superoperator
- Textbook rates
- Textbook and Spinach R1 for first spin
