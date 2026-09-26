# examples/relaxation_theory/aniso_diff_test_2.m

- Signature: `aniso_diff_test_2()`

## Purpose

Relaxation superoperator calculation for an anisotropically shielded two-spin system with an anisotropic rotational diffusion tensor. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Relaxation superoperator calculation for an anisotropically shielded
- two-spin system with an anisotropic rotational diffusion tensor.
- Calculation time: seconds
- Magnet field (Tesla)
- Isotopes
- Chemical shift tensors (ppm)
- Scalar couplings (Hz)
- Difusion tensor eigenvalues
- Relaxation theory
- Basis set
- Spinach housekeeping
- Relaxation superoperator
