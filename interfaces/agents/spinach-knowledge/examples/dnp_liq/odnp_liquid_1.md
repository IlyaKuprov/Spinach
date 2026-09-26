# examples/dnp_liq/odnp_liquid_1.m

- Signature: `odnp_liquid_1()`

## Purpose

Overhauser type DNP in liquid phase at room temperature, using a continu- ous on-resonance CW irradiation of the electron ESR signal. The simulati- on uses Redfield theory to account for the dipolar cross-relaxation. Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Overhauser type DNP in liquid phase at room temperature, using a continu-
- ous on-resonance CW irradiation of the electron ESR signal. The simulati-
- on uses Redfield theory to account for the dipolar cross-relaxation.
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Coordinates (Angstrom)
- Complete basis set
- Relaxation theory
- Spinach housekeeping
- Experiment paramaters
- Simulation
