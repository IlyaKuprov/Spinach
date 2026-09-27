# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_b1_r.m

- Signature: `xix_q_rep_time_ensemble_b1_r()`

## Purpose

Calculates steady-state XiX DNP proton signal versus repetition time, averaging over both electron–proton distance and microwave B1. The source estimates a calculation time of hours.

## Physical and numerical setup

The electron–proton pair is simulated at 80 K and 1.2142 T. Three Gauss–Legendre nodes sample distances from 3.5 to 20 Å, and five nodes sample B1 from 10 to 20 MHz. Thirty repetition times are logarithmically spaced from 10 μs to 1 ms. The pulse train has 36 XiX blocks of 48 ns pulses and uses a fixed −39 MHz microwave offset.

## Calculation and output

For each distance and B1 value, the code sets the coordinates and distance-dependent relaxation, then evaluates the steady state with `powder(...,@xixdnp_steady,...,'esr')`. A `parfor` loop distributes repetition-time calculations. The results are quadrature-averaged over B1 and over distance with the radial (r^2) Jacobian. The script plots the real proton (I_z) expectation versus repetition time and saves `xix_q_rep_time_ensemble_b1_r.fig`.
