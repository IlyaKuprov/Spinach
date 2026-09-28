# examples/dnp_sol/steady_state/top_q_rep_time_ensemble_r.m

- Signature: `top_q_rep_time_ensemble_r()`

## Purpose

Calculates steady-state proton longitudinal expectation versus repetition time while averaging over electron-proton distance. The source comments estimate the calculation takes minutes.

## Physical / mathematical content

The model is an electron-proton pair in a Q-band field (1.2142 T), with spin temperature 80 K, trityl electron g principal values [2.00319, 2.00319, 2.00258], and proton shift [0, 0, 5] ppm. Three Gauss-Legendre distance nodes span 3.5–20 Å. The proton relaxation rate depends on distance and orientation via `r1n_dnp`; the full sphten-liouv basis, diagonal relaxation, dibari equilibrium, and `rep_2ang_800pts_sph` powder grid are used.

The script evaluates 30 logarithmically spaced repetition times from 10 μs to 1 ms. Each uses a fixed 18 MHz electron irradiation, 95 MHz offset, 300 TOP DNP blocks, 10 ns pulses, and 14 ns delays; shot spacing is repetition time minus the pulse-train duration. Each steady-state point calls `powder` with `@topdnp_steady` in `esr` mode. The distance average includes quadrature weights and the radial Jacobian r². The output plots proton Iz expectation versus repetition time and is saved as `top_q_rep_time_ensemble_r.fig`.
