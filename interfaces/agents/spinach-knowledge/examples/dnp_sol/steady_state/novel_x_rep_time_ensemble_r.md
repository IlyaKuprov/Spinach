# examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_r.m

- Signature: `novel_x_rep_time_ensemble_r()`

## Purpose

Calculates steady-state NOVEL DNP signal versus repetition time, averaged over an electron–proton distance distribution. It compares cases without and with a flipback pulse; the source estimates minutes of calculation time.

## Model and method

The X-band model is an electron–proton pair at 0.34 T with trityl g values [2.00319, 2.00319, 2.00258] and spin temperature 80. It uses the `sphten-liouv` formalism without basis approximation, propagator chop tolerance 10^-12, and `t1_t2` relaxation with diagonal retention and DiBari equilibrium. The nuclear R1 rate is evaluated by a distance- and orientation-dependent `r1n_dnp` function handle.

Three Gauss–Legendre nodes span distances of 3.5–20 Å. At each distance, the source creates the electron–proton coordinates and runs powder-averaged steady-state calculations on the `rep_2ang_800pts_sph` grid. The microwave nutation frequency is 15 MHz; the contact pulse is 500 ns, the NOVEL flip-pulse setting is enabled, the added shift is -3.3 MHz, and electron offset is zero. Thirty logarithmically spaced repetition times span 10^-4–10^-2 s. Both flipback conditions are calculated, and the distance average applies quadrature weights multiplied by r^2.

## Output

The figure plots the real proton longitudinal expectation value versus repetition time in milliseconds for the two flipback conditions and saves it as `novel_x_rep_time_ensemble_r.fig`.
