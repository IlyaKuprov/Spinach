# examples/dnp_liq/ccdnp/freq_scan_si_sys_a.m

- Signature: `freq_scan_si_sys_a()`

## Purpose

Steady state nuclear magnetisation as a function of microwave frequency offset and the magnet field in a DNP experiment with two electrons connected by exchange coupling, both coupled to a nucleus by dipolar couplings. Further particulars in: https://doi.org/10.1016/j.jmr.2021.106940

Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, electron-electron scalar exchange, electron-nuclear dipolar couplings, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation uses MATLAB's `parfor` loop to parallelise the calculation across the field grid.

## Implementation structure

- Set up the spin system, interactions, and basis.
- Build the microwave-frequency and magnetic-field grids.
- Use `parfor` over the field grid to compute normalized steady-state DNP responses.
- Plot the real response against magnetic field and microwave-frequency offset.
