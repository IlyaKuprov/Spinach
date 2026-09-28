# examples/dnp_liq/ccdnp/states_vs_tauc_si_sys_b.m

- Signature: `states_vs_tauc_si_sys_b()`

## Purpose

Steady-state amplitudes of selected spin observables as a function of rotational correlation time in a DNP experiment with two electrons connected by exchange coupling, both coupled to a nucleus by dipolar couplings. Further particulars in: https://doi.org/10.1016/j.jmr.2021.106940

Calculation time: seconds

## Physical / mathematical content

- This liquid-state DNP example examines steady-state signal amplitudes as rotational correlation time changes, with electron-electron scalar exchange, electron-nuclear dipolar couplings, and Redfield relaxation.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The implementation uses a `parfor` loop over the correlation-time grid.

## Implementation structure

- Hold the magnetic field and microwave frequency fixed while varying the correlation-time grid, with selected detection operators.
- For each grid point, use `parfor` to update the correlation time, create and basis the spin system, and evaluate `liquid(spin_system,@dnp_freq_scan,locpar,'esr')`.
- Plot absolute steady-state signal amplitudes against rotational correlation time.
