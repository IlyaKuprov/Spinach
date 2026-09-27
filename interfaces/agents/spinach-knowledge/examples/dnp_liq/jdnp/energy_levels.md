# examples/dnp_liq/jdnp/energy_levels.m

- Signature: `energy_levels()`

## Purpose

Plots the two-electron energy levels as the exchange term is varied, illustrating the progression from Zeeman-dominated to exchange-dominated behaviour.

## Model and calculation

The script uses a 14.1 T field and two electron spins with deliberately exaggerated scalar g values of 1.9 and 2.1. It constructs a full Hilbert-space basis, forms the lab-frame Zeeman Hamiltonian and the pairwise (mathbf L_1cdotmathbf L_2) operator, then diagonalises (H_Z-omega_Jmathbf L_1cdotmathbf L_2) at 100 values of (omega_J) spanning (-3omega_E) to (+3omega_E), where (omega_E) is the electron Larmor frequency. Sorted energies are plotted in units of (omega_E) against (omega_J/omega_E).
