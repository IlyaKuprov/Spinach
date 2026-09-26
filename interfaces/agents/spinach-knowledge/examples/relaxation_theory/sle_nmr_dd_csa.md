# examples/relaxation_theory/sle_nmr_dd_csa.m

- Signature: `sle_nmr_dd_csa()`

## Purpose

15N-1H DD-CSA cross-correlation in a protein amide bond spin system using SLE formalism and Redfield relaxation theory. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- 15N-1H DD-CSA cross-correlation in a protein amide bond
- spin system using SLE formalism and Redfield relaxation
- theory.
- Calculation time: seconds
- System specification
- Proximity cut-off
- Basis set
- SLE housekeeping
- SLE parameters
- SLE simulation
- SLE plotting
- BRW parameters
