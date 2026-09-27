# examples/dnp_sol/steady_state/novel_x_rep_time_ensemble_b1.m

- Signature: `novel_x_rep_time_ensemble_b1()`

## Purpose

Calculates steady-state NOVEL DNP signal versus repetition time, averaged over a five-point microwave B1 distribution. It compares acquisition without and with a flipback pulse; the source estimates hours of calculation time.

## Model and method

The model is an electron–proton pair at 0.34 T, with the proton 3.5 Å from the electron, trityl electron g values [2.00319, 2.00319, 2.00258], and spin temperature 80. It uses the `sphten-liouv` formalism without basis approximation, a propagator chop tolerance of 10^-12, and `t1_t2` relaxation with diagonal retention and DiBari equilibrium. The nuclear R1 rate is supplied by a distance- and orientation-dependent `r1n_dnp` function handle.

Five Gauss–Legendre nodes span B1 values from 14 to 16 MHz. Thirty logarithmically spaced repetition times span 10^-4 to 10^-2 s. For each B1 value, the electron pulse duration is recalculated; a 500 ns contact pulse, NOVEL flip-pulse setting, -3.3 MHz added shift, and zero electron offset are used. At each repetition time, powder-averaged steady-state calculations with and without flipback are run on the `rep_2ang_800pts_sph` grid. The B1 results are combined using the quadrature weights.

## Output

The figure plots the real proton longitudinal expectation value versus repetition time in milliseconds for the two flipback conditions and saves it as `novel_x_rep_time_ensemble_b1.fig`.
