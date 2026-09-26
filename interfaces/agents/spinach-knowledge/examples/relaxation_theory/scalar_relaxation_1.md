# examples/relaxation_theory/scalar_relaxation_1.m

- Signature: `scalar_relaxation_1()`

## Purpose

Redfield superoperator for the scalar relaxation of the first kind in a two-proton system with a noisy J-coupling. This si- tuation occurs in aziridines, where the slow nitrogen inversi- on jitters scalar couplings on a millisecond time scale. Set to demonstrate the effect described in: Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Redfield superoperator for the scalar relaxation of the first
- kind in a two-proton system with a noisy J-coupling. This si-
- tuation occurs in aziridines, where the slow nitrogen inversi-
- on jitters scalar couplings on a millisecond time scale. Set
- to demonstrate the effect described in:
- Calculation time: seconds
- System specification
- Basis set
- Relaxation superoperator
- Spinach housekeeping
- Show a spy plot of R
