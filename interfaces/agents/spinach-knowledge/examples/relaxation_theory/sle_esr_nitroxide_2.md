# examples/relaxation_theory/sle_esr_nitroxide_2.m

- Signature: `sle_esr_nitroxide_2()`

## Purpose

Slow motion regime simulation of an ESR spectrum of a nitroxide radical. Set to reproduce Figure 2 from the paper by Concilio et al.: Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Slow motion regime simulation of an ESR spectrum of a nitroxide radical.
- Set to reproduce Figure 2 from the paper by Concilio et al.:
- Calculation time: seconds
- Magnet field
- Isotopes
- Coupling Matrices
- Zeeman Interactions
- Basis set
- Spinach housekeeping
- SLE parameters
- SLE simulation
- SLE plotting
