# examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1.m

- MATLAB implementation: [examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_b1.m)

## Purpose

Use this example to examine the steady-state proton signal versus TPPM contact time when the electron nutation frequency is distributed, with one fixed electron–proton separation. The scan is over TPPM loop count; the plotted horizontal axis is total contact time.

## Model and sequence

Call the no-argument function **tppm_q_con_time_ensemble_b1()** from MATLAB with Spinach available. It builds an electron–proton (E, 1H) Q-band model (sys.magnet=1.2142), with trityl electron Zeeman principal values [2.00319 2.00319 2.00258], proton values [0 0 5] (described in the source as a ppm guess), Euler angles (pi/180)*{[0 10 0],[0 0 10]}, spin temperature 80, and coordinates (0,0,0) and (0,0,3.500) (the source does not label the coordinate unit in this file). The relaxation model is t1_t2, rlx_keep='diagonal', and equilibrium='dibari'; r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet) supplies the orientation- and separation-dependent proton R1, while r1_rates={1000 r1n_rate} and r2_rates={200000 50e3}. Basis settings are sphten-liouv / none, with sys.tols.prop_chop=1e-12 and hygiene disabled.

The TPPM calculation starts with irr_powers=33e6 Hz, which is replaced at each B1 node; it uses spins={E,1H}, pulse_dur=16e-9 s, and electron nutation frequency sampled at six gaussleg(25e6,35e6,5) Hz nodes, grid='rep_2ang_800pts_sph', second-pulse phase 120*pi/180, addshift=-13e6, and el_offs=2e6. For each B1 node, powder(spin_system,@xixdnp_steady,localpar,'esr') evaluates loop counts 1:256; the shot spacing is 816e-6 minus the duration of the two pulse trains.

## Result and scope

B1-node results are combined with the Gauss–Legendre weights and normalised. The figure plots the real part of the proton Lz expectation value against 2*pulse_dur*loop_counts in microseconds and saves tppm_q_con_time_ensemble_b1.fig. The source describes the calculation as taking hours. It is a finite six-node B1 quadrature at one fixed separation; the function saves a figure, not a separate numerical data file. The source does not state units for the magnet, spin temperature, relaxation-rate entries, addshift, or el_offs.
