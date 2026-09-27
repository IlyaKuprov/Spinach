# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T2e.m

- Signature: `xix_q_rep_time_ensemble_r_T2e()`

## Purpose

Compares XiX DNP proton repetition-time profiles for five electron transverse relaxation times (T_{2e}), with an electron–proton distance average. The source estimates a calculation time of hours.

## Physical and numerical setup

The script compares (T_{2e}=50, 15, 5, 1.5, 0.5) μs and sets the electron R2 rate to (1/T_{2e}); the proton R2 rate is 50e3 s⁻¹. Three Gauss–Legendre nodes sample distances from 3.5 to 20 Å. The spin system is at 80 K and 1.2142 T, with 18 MHz electron nutation frequency. Thirty repetition times span 10 μs to 1 ms logarithmically; the XiX train has 36 blocks of 48 ns pulses and a −39 MHz offset.

## Calculation and output

For each T2e and distance, the source evaluates the steady state with `powder(...,@xixdnp_steady,...,'esr')`. Distance results are averaged with quadrature weights and the radial (r^2) Jacobian. The script plots the negative real proton (I_z) expectation against repetition time, adds a T2e legend, and saves `xix_q_rep_time_ensemble_r_T2e.fig`.
