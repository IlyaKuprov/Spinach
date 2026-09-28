# examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1_r.m

- Signature: `novel_x_rep_time_ensemble_b1_r()`

## Purpose

Calculates steady-state NOVEL DNP signal versus repetition time while averaging over electron–proton distance and microwave B1 distributions. It compares cases without and with a flipback pulse; the source estimates hours of calculation time.

## Model and method

The X-band model is an electron–proton pair at 0.34 T with trityl g values [2.00319, 2.00319, 2.00258] and spin temperature 80. It uses the `sphten-liouv` formalism without basis approximation, propagator chop tolerance 10^-12, and `t1_t2` relaxation with diagonal retention and DiBari equilibrium. The nuclear R1 rate is evaluated using a distance- and orientation-dependent `r1n_dnp` function handle.

Three Gauss–Legendre distance nodes span 3.5–20 Å and five B1 nodes span 14–16 MHz. Thirty logarithmically spaced repetition times span 10^-4–10^-2 s. At each distance, the source builds the pair coordinates and calculates powder-averaged steady states on the `rep_2ang_800pts_sph` grid, for both flipback conditions at each B1 value. The proton contact duration is 500 ns, the NOVEL flip-pulse setting is enabled, the added shift is -3.3 MHz, and electron offset is zero. It integrates over B1 with quadrature weights and over distance with quadrature weights multiplied by r^2.

## Output

The figure plots the real proton longitudinal expectation value versus repetition time in milliseconds for the two flipback conditions and saves it as `novel_x_rep_time_ensemble_b1_r.fig`.
