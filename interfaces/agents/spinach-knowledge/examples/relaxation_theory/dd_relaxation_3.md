# examples/relaxation_theory/dd_relaxation_3.m

- Signature: `dd_relaxation_3()`

## Purpose

Extreme narrowing limit case comparison between the dipolar relaxation rates in proton-proton and proton-deuterium system. The rate must sca- le with the square of the magnetogyric ratio and with S(S+1), where S is the quantum number of the partner spin. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Extreme narrowing limit case comparison between the dipolar relaxation
- rates in proton-proton and proton-deuterium system. The rate must sca-
- le with the square of the magnetogyric ratio and with S(S+1), where S
- is the quantum number of the partner spin.
- Calculation time: seconds
- % System with two protons
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Relaxation superoperator
- Longitudinal relaxation rate
