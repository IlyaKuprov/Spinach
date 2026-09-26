# examples/relaxation_theory/hfc_antisymm_1.m

- Signature: `hfc_antisymm_1()`

## Purpose

Longitudinal and transverse relaxation rates in a system with a significant antisymmetry in the hyperfine tensor. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Longitudinal and transverse relaxation rates in a system
- with a significant antisymmetry in the hyperfine tensor.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Relaxation superoperator
- Textbook rates
- Textbook and Spinach R1 for first spin
- Textbook and Spinach R1 for second spin
- Textbook and Spinach R2 for first spin
