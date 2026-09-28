# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r.m

- Signature: `xix_q_rep_time_ensemble_r()`

## Purpose

Calculates steady-state XiX DNP proton signal versus repetition time, averaged over electron–proton distance. The source estimates a calculation time of minutes.

## Physical and numerical setup

The spin system is an electron–proton pair at 80 K and 1.2142 T. Three Gauss–Legendre nodes cover distances from 3.5 to 20 Å; the electron nutation frequency is fixed at 18 MHz. Thirty repetition times are logarithmically spaced from 10 μs to 1 ms. The XiX train has 36 blocks with 48 ns pulses and a fixed microwave offset of −39 MHz.

## Calculation and output

At each distance the source updates the coordinates and the distance-/orientation-dependent nuclear relaxation, then evaluates the steady state with `powder(...,@xixdnp_steady,...,'esr')`. It averages over distance quadrature weights including the radial (r^2) Jacobian, plots the real proton (I_z) expectation versus repetition time, and saves `xix_q_rep_time_ensemble_r.fig`.
