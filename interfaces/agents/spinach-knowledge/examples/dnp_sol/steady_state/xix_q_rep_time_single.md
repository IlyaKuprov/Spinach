# examples/dnp_sol/steady_state/xix_q_rep_time_single.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_q_rep_time_single.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_rep_time_single.m)

## Purpose

Run the no-argument MATLAB function `xix_q_rep_time_single()` with Spinach and its example helpers available; it evaluates a steady-state XiX DNP response as repetition time is varied for an electron–proton pair. The source estimates a calculation time of seconds.

## Model and scan

The setup labels its magnet Q-band and sets `sys.magnet=1.2142`; the spin-temperature value is 80. It uses an electron trityl g-tensor [2.00319 2.00319 2.00258], a proton Zeeman-shift vector [0 0 5] (the source labels the proton value a ppm guess), and Euler rotations of [0 10 0] and [0 0 10] degrees. The two spin coordinates are at z = 0 and 3.500; the source does not annotate a unit for these coordinates. Proton longitudinal relaxation is supplied by `r1n_dnp` as a function of orientation and the fixed electron–nuclear separation; its call uses `2.00230`, `1e-3`, and `52` as additional source parameters. The source sets `r1_rates={1000 r1n_rate}`, `r2_rates={200000 50e3}`, `t1_t2` relaxation, and `dibari` equilibrium.

The driver uses the `rep_2ang_800pts_sph` powder grid, `sphten-liouv` with no basis approximation, and a propagator chop tolerance of `1e-12`. Its XiX settings are a 48 ns pulse, 36 blocks, an inverted second-pulse phase (`pi`), electron nutation frequency 18e6 Hz, additional shift −13e6, and electron offset −39e6. It evaluates 30 logarithmically spaced repetition times from 1e−5 to 1e−3 s; for each point, shot spacing is the repetition time minus `2*nloops*pulse_dur`. A `parfor` loop calls `powder(spin_system,@xixdnp_steady,localpar,'esr')` at each setting.

## Output and dependencies

The plotted quantity is the real part of the returned DNP value, labelled as the proton `Lz` expectation value, against repetition time in ms. The figure is saved as `xix_q_rep_time_single.fig`. The driver uses Spinach system/basis/detection and powder functions plus the example helpers `r1n_dnp` and `xixdnp_steady`; it writes a MATLAB figure rather than a separate numeric data file.
