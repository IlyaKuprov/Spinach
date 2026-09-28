# examples/dnp_liq/ccdnp/rates_si_sys_b.m

- Signature: `rates_si_sys_b()`

## Purpose

Self-and cross-relaxation rates in cross-correlated DNP, considering a system with two electrons connected by exchange coupling, both coupled to a nucleus by dipolar couplings. Further particulars in: https://doi.org/10.1016/j.jmr.2021.106940

Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, electron-electron scalar exchange, electron-nuclear dipolar couplings, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then calculate the relaxation superoperator and print selected rate matrix elements.

## Implementation structure

- Set up the spin system and basis with Zeeman interactions, electron-electron scalar exchange, and coordinate-derived electron-nuclear dipolar couplings.
- Compute the relaxation superoperator.
- Print selected self- and cross-relaxation rates from its matrix elements.
