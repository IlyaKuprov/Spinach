# examples/dnp_sol/steady_state/novel_x_rep_time_single.m

- Signature: `novel_x_rep_time_single()`

## Purpose

Calculates steady-state NOVEL DNP proton signal versus repetition time for a single electron–proton pair. It compares cases without and with a flipback pulse; the source estimates minutes of calculation time.

## Model and method

The X-band model uses a 0.34 T field, trityl g values [2.00319, 2.00319, 2.00258], spin temperature 80, and an electron–proton separation of 3.5 Å. It uses the `sphten-liouv` formalism without basis approximation, propagator chop tolerance 10^-12, and `t1_t2` relaxation with diagonal retention and DiBari equilibrium. The nuclear R1 rate is computed through a distance- and orientation-dependent `r1n_dnp` function handle.

The microwave nutation frequency is 15 MHz; the contact pulse is 500 ns, a NOVEL flip-pulse setting is enabled, the added shift is -3.3 MHz, and electron offset is zero. Thirty logarithmically spaced repetition times span 10^-4–10^-2 s. At each point the source runs powder-averaged steady-state calculations on the `rep_2ang_800pts_sph` grid for both flipback conditions, using `noveldnp_steady` with the ESR context.

## Output

The figure plots the real proton longitudinal expectation value versus repetition time in milliseconds for the two flipback conditions and saves it as `novel_x_rep_time_single.fig`.
