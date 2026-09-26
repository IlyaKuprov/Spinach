# examples/relaxation_theory/t1t2_strychnine.m

- Signature: `t1t2_strychnine()`

## Purpose

Relaxation analysis for strychnine, dipolar processes only. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Relaxation analysis for strychnine, dipolar processes only.
- Calculation time: seconds.
- Spin system properties
- Magnet field
- Basis set
- Relaxation theory parameters
- Distance cut-off
- Spinach housekeeping
- Relaxation analysis
