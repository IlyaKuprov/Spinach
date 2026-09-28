# examples/dnp_sol/steady_state/top_q_con_time_ensemble_b1.m

- Signature: `top_q_con_time_ensemble_b1()`

## Purpose

Calculates the proton longitudinal expectation value after TOP DNP steady-state simulations as a function of total contact time, averaging separately over two electron-Rabi-frequency distributions.

## Physical / mathematical content

The model contains an electron and a proton at 3.5 Å in a Q-band field (1.2142 T), with an 80 K spin temperature. The electron Zeeman tensor is specified by principal values [2.00319, 2.00319, 2.00258]; the proton uses a [0, 0, 5] ppm shift. It uses the full sphten-liouv basis, diagonal relaxation superoperators, and dibari equilibrium. The powder grid is `rep_2ang_800pts_sph`.

For each ensemble, the script samples five B1 values by Gaussian quadrature: 10–20 MHz and 25–35 MHz. It evaluates loop counts 1 through 256, using 10 ns pulses separated by 14 ns delays and ensemble-specific electron offsets (95 MHz and 92 MHz). At each point it calls `powder` with `@topdnp_steady` and the `esr` mode, then takes the quadrature-weighted B1 average. The output plot compares the proton Iz expectation versus total contact time for the two ensembles and is saved as `top_q_con_time_ensemble_b1.fig`.
