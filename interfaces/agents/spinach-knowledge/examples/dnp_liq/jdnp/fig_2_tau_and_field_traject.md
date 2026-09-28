# examples/dnp_liq/jdnp/fig_2_tau_and_field_traject.m

- Signature: `fig_2_tau_and_field_traject()`
- Reference: [Physical Chemistry Chemical Physics, DOI: 10.1039/D1CP04186J](https://doi.org/10.1039/d1cp04186j)
- Calculation time: seconds (per source comment)

## Purpose

Plots time-dependent proton DNP for six static fields and four rotational correlation times using the spin system from `system_specification()`.

## Calculation

At each field (0.5, 3.4, 7.0, 11.7, 14.1, and 23.5 T), the script sets the microwave offset from the trityl and free-electron frequencies and sets the electron–electron scalar coupling to the sum of the isotropic electron and proton Zeeman frequencies. For each correlation time (300, 400, 500, and 600 ps), it constructs the Spinach system and basis, forms the ESR Hamiltonian plus relaxation superoperator and microwave terms, and propagates the thermal-equilibrium state while observing proton (L_z).

The 200 ms trajectories use 1 ms time steps. The proton signal is normalized by its equilibrium expectation value and plotted in a separate field panel, with one curve per correlation time.
