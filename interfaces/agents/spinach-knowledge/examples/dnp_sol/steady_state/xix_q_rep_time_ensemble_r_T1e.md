# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T1e.m

- Signature: `xix_q_rep_time_ensemble_r_T1e()`

## Purpose

Compares XiX DNP proton repetition-time profiles for five electron longitudinal relaxation times (T_{1e}), with an electron–proton distance average. The source estimates a calculation time of hours.

## Physical and numerical setup

The script compares (T_{1e}=10, 3, 1, 0.3, 0.1) ms. For each value it uses (R_{1e}=1/T_{1e}), evaluates the source's distance-/orientation-dependent nuclear relaxation function, and samples three distance nodes from 3.5 to 20 Å. The spin system is at 80 K and 1.2142 T; the electron nutation frequency is 18 MHz. Thirty repetition times span 10 μs to 1 ms logarithmically. The XiX train has 36 blocks of 48 ns pulses and a −39 MHz offset.

## Calculation and output

For each T1e and distance, the steady state is calculated over repetition time with `powder(...,@xixdnp_steady,...,'esr')`; distance results are averaged with quadrature weights and the radial (r^2) Jacobian. The script overlays the negative real proton (I_z) expectation against repetition time, adds a T1e legend, and saves `xix_q_rep_time_ensemble_r_T1e.fig`.
