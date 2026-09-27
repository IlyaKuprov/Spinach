# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_b1.m

- Signature: `xix_q_rep_time_ensemble_b1()`

## Purpose

Calculates the steady-state XiX DNP proton signal versus shot repetition time, averaged over a microwave B1 distribution. The source estimates a calculation time of minutes.

## Physical and numerical setup

The electron–proton pair is simulated at 80 K and 1.2142 T with a fixed 3.5 Å separation. Five Gauss–Legendre nodes sample electron nutation frequencies from 10 to 20 MHz. Thirty repetition times are logarithmically spaced from 10 μs to 1 ms. The pulse train uses 36 XiX blocks, 48 ns pulses, and a fixed microwave offset of −39 MHz.

## Calculation and output

For each B1 node, the script runs the steady-state powder calculation over repetition times with MATLAB `parfor`. Shot spacing is the repetition time minus the total pulse-train duration; each result comes from `powder(...,@xixdnp_steady,...,'esr')`. It averages the signal using the B1 quadrature weights, plots the real proton (I_z) expectation against repetition time, and saves `xix_q_rep_time_ensemble_b1.fig`.
