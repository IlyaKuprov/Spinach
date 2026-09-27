# examples/dnp_sol/steady_state/top_q_rep_time_ensemble_b1_r.m

- Signature: `top_q_rep_time_ensemble_b1_r()`

## Purpose

Calculates steady-state proton longitudinal expectation versus repetition time, averaging over electron-proton distance and microwave B1 distributions. The source comments estimate the calculation takes hours.

## Physical / mathematical content

The model is an electron-proton pair at Q-band (1.2142 T), with spin temperature 80 K, trityl electron g principal values [2.00319, 2.00319, 2.00258], and proton shift [0, 0, 5] ppm. Three Gauss-Legendre distance nodes span 3.5–20 Å, and five B1 nodes span 10–20 MHz. The distance-dependent proton relaxation rate uses `r1n_dnp`; the full sphten-liouv basis, diagonal relaxation, dibari equilibrium, and `rep_2ang_800pts_sph` powder grid are used.

The script evaluates 30 logarithmically spaced repetition times from 10 μs to 1 ms. For each distance and B1 node it uses 300 TOP DNP blocks, 10 ns pulses, 14 ns delays, a −13 MHz added shift, and 95 MHz electron offset; shot spacing is repetition time minus the pulse-train duration. It calls `powder` with `@topdnp_steady` in `esr` mode, averages over B1 weights, then averages over distance with the radial Jacobian r². The output plots proton Iz expectation versus repetition time and is saved as `top_q_rep_time_ensemble_b1_r.fig`.
