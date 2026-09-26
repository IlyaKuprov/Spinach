# examples/dnp_liq/ccdnp/rates_si_sys_b.m

- Signature: `rates_si_sys_b()`

## Purpose

Self-and cross-relaxation rates in cross-correlated DNP, considering a system with two electrons connected by exchange coupling, both cou- pled to a nucleus by dipolar couplings. Further particulars in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Self-and cross-relaxation rates in cross-correlated DNP, considering
- a system with two electrons connected by exchange coupling, both cou-
- pled to a nucleus by dipolar couplings. Further particulars in:
- Calculation time: seconds
- Magnet field
- Spin system
- Zeeman interactions
- Exchange coupling
- coordinates
- Basis set
- Relaxation theory
- Spinach housekeeping
