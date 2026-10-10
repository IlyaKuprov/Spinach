# examples/dnp_sol/steady_state/top_q_con_time_ensemble_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/top_q_con_time_ensemble_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/top_q_con_time_ensemble_r.m)

- Signature: `top_q_con_time_ensemble_r()`

## Question and model

How does the steady-state proton longitudinal polarisation vary with TOP contact time for two fixed irradiation conditions after averaging over electron–proton distance? Unlike the `ensemble_b1_r` variant, this script has no B1 ensemble. It uses a Q-band magnet setting of 1.2142, trityl electron g principal values `[2.00319 2.00319 2.00258]`, proton shift values `[0 0 5]`, Euler angles `(pi/180)*{[0 10 0],[0 0 10]}`, and spin temperature 80 K.

## Scan and averaging

The scan uses `nloops=1:256` TOP blocks, with 10 ns pulse duration and 14 ns delay (24 ns per block). The electron–proton distance is integrated with 3 Gauss–Legendre points from 3.5 to 20 Å. At each distance and loop count, `powder(spin_system,@topdnp_steady,localpar,'esr')` calculates each fixed setting: A uses 18 MHz irradiation, 95 MHz electron offset, and 102 μs minus pulse-train duration for shot spacing; B uses 33 MHz, 92 MHz, and 153 μs minus pulse-train duration. The detected operator is proton `Lz`, and `r1n_dnp` provides the distance- and orientation-dependent proton R1. The source assigns R1 entries `1e3` and R2 values `200e3` and `50e3` (units are not annotated), uses spins `E` and `1H` with grid `rep_2ang_800pts_sph`, and sets `addshift=-13e6`.

The basis is `sphten-liouv` without approximation; propagator chopping tolerance is `1e-12`, and `hygiene` is disabled. The relaxation model is `t1_t2`, with diagonal relaxation terms retained and the `dibari` equilibrium setting.

## Output and limits

Each setting is averaged over distance using the quadrature weights and radial Jacobian `r^2`. The real proton `I_z` expectation is plotted against total contact time, with curves labelled TOP, 18 MHz and TOP, 33 MHz, and saved as `top_q_con_time_ensemble_r.fig`. The source estimates hours of calculation; only a figure is saved, not a numeric results table.
