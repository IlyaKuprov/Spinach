# examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1.m

- Signature: `top_q_rep_time_ensemble_b1()`

## Purpose

Calculates steady-state proton longitudinal expectation versus repetition time, averaging over a microwave B1 distribution. The source comments estimate the calculation takes minutes.

## Physical / mathematical content

The model is an electron-proton pair at 3.5 Å in a Q-band field (1.2142 T), with spin temperature 80 K, trityl electron g principal values [2.00319, 2.00319, 2.00258], and proton shift [0, 0, 5] ppm. It uses the full sphten-liouv basis, distance- and orientation-dependent proton relaxation via `r1n_dnp`, diagonal relaxation, dibari equilibrium, and the `rep_2ang_800pts_sph` powder grid.

The script samples five B1 values by Gauss-Legendre quadrature over 10–20 MHz and evaluates 30 logarithmically spaced repetition times from 10 μs to 1 ms. Each calculation uses 300 TOP DNP blocks, 10 ns pulses, 14 ns delays, a −13 MHz added shift, and 95 MHz electron offset; shot spacing is repetition time minus the pulse-train duration. It calls `powder` with `@topdnp_steady` in `esr` mode and averages over B1 weights. The output plots proton Iz expectation versus repetition time and is saved as `top_q_rep_time_ensemble_b1.fig`.
