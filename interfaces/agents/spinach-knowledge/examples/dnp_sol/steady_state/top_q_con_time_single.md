# examples/dnp_sol/steady_state/top_q_con_time_single.m

- Signature: `top_q_con_time_single()`

## Purpose

Calculates proton longitudinal expectation versus total TOP DNP contact time for two fixed irradiation parameter sets at a fixed electron-proton separation. The source comments estimate the calculation takes hours.

## Physical / mathematical content

The model is an electron-proton pair at 3.5 Å in a Q-band field (1.2142 T), with spin temperature 80 K, trityl electron g principal values [2.00319, 2.00319, 2.00258], and proton shift [0, 0, 5] ppm. It uses distance- and orientation-dependent proton relaxation through `r1n_dnp`, the full sphten-liouv basis, diagonal relaxation, dibari equilibrium, and the `rep_2ang_800pts_sph` powder grid.

For loop counts 1–256, the script uses 10 ns pulses and 14 ns delays and runs two steady-state calculations per count with `powder`, `@topdnp_steady`, and `esr` mode. Set A uses 18 MHz irradiation, 95 MHz electron offset, and shot spacing 102 μs minus the pulse-train duration; set B uses 33 MHz, 92 MHz, and 153 μs minus that duration. The plot compares proton Iz expectation against total contact time and is saved as `top_q_con_time_single.fig`.
