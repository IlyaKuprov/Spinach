# examples/dnp_sol/steady_state/top_q_nutation_ensemble_b1_r.m

- Signature: `top_q_nutation_ensemble_b1_r()`

## Purpose

Builds steady-state TOP DNP field profiles across six selected nutation frequencies while averaging over electron-proton distance and B1 distributions. The source comments estimate the calculation takes minutes.

## Physical / mathematical content

The model is an electron-proton pair at Q-band (1.2142 T), with spin temperature 80 K, trityl electron g principal values [2.00319, 2.00319, 2.00258], and proton shift [0, 0, 5] ppm. Three Gauss-Legendre distance nodes span 3.5–20 Å, and the five B1 nodes span 0.2ν–1.2ν for each selected nutation frequency ν. The script evaluates ten microwave offsets from 88 to 97 MHz, with 300 TOP DNP blocks, 10 ns pulses, 14 ns delays, and shot spacing set to 153 μs minus the pulse-train duration. Distance-dependent proton relaxation uses `r1n_dnp`; the calculation uses the full sphten-liouv basis, diagonal relaxation, dibari equilibrium, and `rep_2ang_800pts_sph` powder grid.

For each ν in [6.8, 9.6, 13.5, 17.5, 25, 36] MHz, the local function `top_field_profile_b1_r` evaluates the steady state with `powder` and `@topdnp_steady` in `esr` mode. It averages over B1 quadrature weights and over distance using the radial Jacobian r², then adds a curve to the 3-D plot. The figure is saved as `top_q_nutation_ensemble_b1_r.fig`.
