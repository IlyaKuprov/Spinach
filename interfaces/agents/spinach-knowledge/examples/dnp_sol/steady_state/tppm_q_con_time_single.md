# examples/dnp_sol/steady_state/tppm_q_con_time_single.m

- MATLAB implementation: [examples/dnp_sol/steady_state/tppm_q_con_time_single.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/tppm_q_con_time_single.m)

## Purpose

Use this shorter baseline to inspect the steady-state proton signal versus TPPM contact time at one electron–proton separation and one electron nutation frequency. It sweeps the number of TPPM loops, rather than averaging over distance or B1; the source estimates minutes of calculation time.

## Model and sequence

Call the no-argument function **tppm_q_con_time_single()** from MATLAB with Spinach available. It builds an electron–proton (E, 1H) Q-band model (sys.magnet=1.2142) with trityl g principal values [2.00319 2.00319 2.00258], proton values [0 0 5] (described in the source as a ppm guess), Euler angles (pi/180)*{[0 10 0],[0 0 10]}, spin temperature 80, and coordinates (0,0,0) and (0,0,3.500). r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet) supplies the proton R1 function. Relaxation is t1_t2, r1_rates={1e3 r1n_rate}, r2_rates={200e3 50e3}, rlx_keep='diagonal', and equilibrium='dibari'. The basis is sphten-liouv / none; prop_chop=1e-12.

The TPPM parameters are spins={E,1H}, electron nutation frequency 33e6 Hz, pulse_dur=16e-9 s, grid='rep_2ang_800pts_sph', second-pulse phase 120*pi/180, addshift=-13e6, and el_offs=2e6. For each loop count 1:256, the example calls powder(spin_system,@xixdnp_steady,parameters,'esr'); shot spacing is 816e-6 less the two-train pulse duration.

## Result and scope

The saved figure, tppm_q_con_time_single.fig, plots the real proton Lz expectation value against total contact time 2*pulse_dur*loop_counts in microseconds. The function does not save a separate numerical data file. This is a single geometry and single B1 setting, not an ensemble calculation. The source does not label the coordinate unit in this file, and it does not state units for magnet, temperature, relaxation-rate entries, addshift, or el_offs.
