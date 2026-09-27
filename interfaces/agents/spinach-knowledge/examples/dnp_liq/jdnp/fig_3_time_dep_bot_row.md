# examples/dnp_liq/jdnp/fig_3_time_dep_bot_row.m

- Signature: `fig_3_time_dep_bot_row()`
- Reference: [Physical Chemistry Chemical Physics, DOI: 10.1039/D1CP04186J](https://doi.org/10.1039/d1cp04186j)
- Calculation time: seconds (per source comment)

## Purpose

Plots the time-dependent proton DNP signal for the two-electron system from `system_specification()`, at three static fields. The source identifies this as the bottom-row simulation; the paired top-row example removes one electron to demonstrate the contrast.

## Calculation

For fields 0.034, 0.34, and 3.4 T, the script sets the microwave offset to the trityl–free-electron frequency difference and sets the scalar electron–electron coupling to the sum of the isotropic electron and proton Zeeman frequencies. It builds the Spinach system and basis, adds the microwave drive and offset to the ESR Hamiltonian, includes the relaxation superoperator, and propagates from thermal equilibrium while observing proton (L_z).

Each 300 ms trajectory is sampled at 1 ms intervals and normalized by the proton equilibrium signal, then plotted in its own field panel.
