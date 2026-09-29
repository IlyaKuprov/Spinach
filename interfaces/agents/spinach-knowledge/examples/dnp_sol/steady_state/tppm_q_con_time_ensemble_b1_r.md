# examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1_r.m)

## Purpose

Choose this variant when the steady-state proton signal versus TPPM contact time must include both an electron-nutation-frequency distribution and an electron–proton distance distribution. Relative to tppm_q_con_time_ensemble_b1.m, the added scan dimension is separation; relative to tppm_q_con_time_ensemble_r.m, this also averages B1.

## Model and sequence

Call the no-argument function **tppm_q_con_time_ensemble_b1_r()** from MATLAB with Spinach available. It uses a Q-band electron–proton (E, 1H) model: sys.magnet=1.2142; trityl g principal values [2.00319 2.00319 2.00258]; proton values [0 0 5] (source-described ppm guess); Euler angles (pi/180)*{[0 10 0],[0 0 10]}; spin temperature 80; sphten-liouv / none basis; prop_chop=1e-12; and hygiene disabled. At each distance node it sets coordinates (0,0,0) and (0,0,r) and evaluates r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r,bet) as the proton R1 function. Relaxation is t1_t2, with r1_rates={1000 r1n_rate}, r2_rates={200000 50e3}, diagonal relaxation retention, and dibari equilibrium.

`gaussleg(3.5,20,3)` produces four distance nodes (ångström, as commented in the source); `gaussleg(25e6,35e6,5)` produces six B1 nodes in Hz. The initial irr_powers=33e6 Hz setting is replaced by each B1 node. At each pair of nodes it runs powder(spin_system,@xixdnp_steady,localpar,'esr') for loop counts 1:256. Sequence settings are pulse_dur=16e-9 s, grid='rep_2ang_800pts_sph', phase 120*pi/180 for the second pulse, addshift=-13e6, and el_offs=2e6; shot spacing is 816e-6 less the two-train pulse duration.

## Result and scope

The source first averages over B1 weights, then averages over distance weights multiplied by r^2 (the radial Jacobian), normalising each weighted sum. It plots the real proton Lz expectation value against total contact time 2*pulse_dur*loop_counts in microseconds and saves tppm_q_con_time_ensemble_b1_r.fig. The source describes the calculation as taking hours. This is a four-node distance and six-node B1 quadrature (24 distance/B1 pairs), not a continuous distribution; no separate numerical data file is saved. Units for magnet, temperature, relaxation-rate entries, addshift, and el_offs are not stated in this source.
