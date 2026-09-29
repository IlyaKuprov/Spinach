# examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/tppm_q_con_time_ensemble_r.m)

## Purpose

Use this example for the steady-state proton signal versus TPPM contact time averaged over electron–proton separation, with a single fixed electron nutation frequency. It varies loop count for the contact-time curve and integrates over the distance ensemble; unlike the ensemble_b1 variants, it does not sample B1.

## Model and sequence

Call the no-argument function **tppm_q_con_time_ensemble_r()** from MATLAB with Spinach available. The model contains an electron and proton (E, 1H) at Q band (sys.magnet=1.2142) with trityl g principal values [2.00319 2.00319 2.00258], proton values [0 0 5] (source-described ppm guess), Euler angles (pi/180)*{[0 10 0],[0 0 10]}, and spin temperature 80. It uses sphten-liouv / none basis, prop_chop=1e-12, and disables hygiene. For each distance node, coordinates are (0,0,0) and (0,0,r); r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r,bet) supplies orientation- and distance-dependent proton R1. Relaxation is t1_t2, with r1_rates={1e3 r1n_rate}, r2_rates={200e3 50e3}, rlx_keep='diagonal', and equilibrium='dibari'.

The distance quadrature is gaussleg(3.5,20,3) Angstrom. For each node, the experiment uses irr_powers=33e6 Hz, pulse_dur=16e-9 s, grid='rep_2ang_800pts_sph', second-pulse phase 120*pi/180, addshift=-13e6, and el_offs=2e6. It runs powder(spin_system,@xixdnp_steady,localpar,'esr') over loop counts 1:256; shot spacing is 816e-6 minus the duration of the two pulse trains.

## Result and scope

The node results are averaged with the Gauss–Legendre distance weights multiplied by r^2 for the radial Jacobian and normalised by the weighted sum. The plot is the real proton Lz expectation value versus total contact time 2*pulse_dur*loop_counts in microseconds; the figure is saved as tppm_q_con_time_ensemble_r.fig. The source describes the calculation as taking hours. The distance distribution is represented by four quadrature nodes; no separate numerical data file is saved. The source does not state units for magnet, temperature, relaxation-rate entries, addshift, or el_offs.
