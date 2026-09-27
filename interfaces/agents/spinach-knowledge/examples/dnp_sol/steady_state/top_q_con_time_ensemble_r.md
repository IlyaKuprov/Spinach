# examples/dnp_sol/steady_state/top_q_con_time_ensemble_r.m

- Signature: `top_q_con_time_ensemble_r()`

## Purpose

Calculates proton longitudinal expectation versus total TOP DNP contact time for two parameter sets, averaging each result over an electron-proton distance distribution. The source comments estimate the calculation takes hours.

## Physical / mathematical content

The model contains an electron and a proton in a Q-band field (1.2142 T), with trityl electron g principal values [2.00319, 2.00319, 2.00258], a [0, 0, 5] ppm proton shift, and spin temperature 80 K. Three Gauss-Legendre distance nodes span 3.5–20 Å; each pair is placed on z, and the proton relaxation rate depends on distance and orientation through `r1n_dnp`. The calculation uses the full sphten-liouv basis, diagonal relaxation, dibari equilibrium, and the `rep_2ang_800pts_sph` powder grid.

For each distance, the script runs loop counts 1–256 with 10 ns pulses and 14 ns delays. Parameter set A uses 18 MHz electron irradiation, 95 MHz offset, and shot spacing 102 μs minus the pulse-train duration; set B uses 33 MHz, 92 MHz, and 153 μs minus that duration. Each point calls `powder` with `@topdnp_steady` in `esr` mode. The distance average includes the radial Jacobian r² and quadrature weights. The plot compares proton Iz versus contact time for the two parameter sets and is saved as `top_q_con_time_ensemble_r.fig`.
