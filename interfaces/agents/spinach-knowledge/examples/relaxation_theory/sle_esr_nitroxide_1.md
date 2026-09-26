# examples/relaxation_theory/sle_esr_nitroxide_1.m

- Signature: `sle_esr_nitroxide_1()`

## Purpose

Comparison between nitroxide simulation using SLE formalism and Redfield relaxation theory. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Comparison between nitroxide simulation using SLE formalism
- and Redfield relaxation theory.
- Calculation time: seconds
- Spin system properties
- Magnet induction
- Proximity cut-off
- Basis set
- SLE housekeeping
- SLE parameters
- SLE simulation
- SLE plotting
- BRW parameters
