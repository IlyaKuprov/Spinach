# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1n.m

- Signature: `xix_q_rep_time_ensemble_r_T1n()`

## Purpose

Compares XiX DNP proton repetition-time profiles for five nuclear longitudinal relaxation times (T_{1n}), with an electron–proton distance average. The source estimates a calculation time of hours.

## Physical and numerical setup

The script compares (T_{1n}=50, 5, 0.5, 0.05, 0.005) s and sets the nuclear R1 rate to (1/T_{1n}), while the electron R1 rate is fixed at 1e3 s⁻¹. Three Gauss–Legendre nodes sample distances from 3.5 to 20 Å. The spin system is at 80 K and 1.2142 T, with 18 MHz electron nutation frequency. Thirty repetition times span 10 μs to 1 ms logarithmically; the XiX train has 36 blocks of 48 ns pulses and a −39 MHz offset.

## Calculation and output

For each T1n and distance, the source evaluates the steady state with `powder(...,@xixdnp_steady,...,'esr')`. It averages over the distance quadrature with the radial (r^2) Jacobian, plots the negative real proton (I_z) expectation against repetition time, and saves `xix_q_rep_time_ensemble_r_T1n.fig` with a legend for the five T1n values.
