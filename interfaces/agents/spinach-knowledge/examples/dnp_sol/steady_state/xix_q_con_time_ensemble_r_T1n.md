# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1n.m

Signature: `xix_q_con_time_ensemble_r_T1n()`

This example compares steady-state XiX DNP proton-contact curves while varying nuclear T1. Each curve averages over an electron–proton distance ensemble; the MATLAB source estimates the calculation time as hours.

## Setup and scan

The model is an E–¹H pair with the source's Q-band magnet setting `sys.magnet=1.2142`, spin temperature `80`, trityl g values `[2.00319 2.00319 2.00258]`, and proton Zeeman entry `[0 0 5]` (described in the source as a ppm guess). Euler-angle entries are `[0 10 0]` and `[0 0 10]` degrees. The basis is the full `sphten-liouv` basis (`approximation='none'`); propagator chopping tolerance is `1e-12`; the source also sets `sys.disable={'hygiene'}`.

The scan uses nuclear T1 values `[50 5 0.5 0.05 0.005]` seconds. The source sets `inter.r1_rates={1e3,1/T1n}`, `inter.r2_rates={200e3,50e3}`, diagonal relaxation retention, and Di Bari equilibrium. Unlike the T2-ensemble siblings, this variant sets a fixed nuclear T1 rate rather than calling `r1n_dnp`.

For each T1, `gaussleg(3.5,20,3)` supplies four distance nodes over 3.5–20 Å; each is placed on the z axis and evaluated separately. The source integrates those values using the returned weights and the radial `r^2` Jacobian. The contact scan is 1–64 XiX loops, each using two 48 ns pulses; plotted contact time is `2*pulse_dur*nloops`, displayed in μs. The pulse train uses 153 μs shot spacing minus its own duration. Other fixed experiment settings are 18e6 electron nutation frequency (source specifies Hz), `phase=pi` (inverted second pulse), `addshift=-13e6`, `el_offs=61e6`, and grid `rep_2ang_800pts_sph`.

## Run dependencies and output

Run from a Spinach MATLAB environment with `gaussleg`, `powder`, and the usual Spinach system/basis/state and plotting functions available. The steady-state experiment is supplied by [`xixdnp_steady`](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/hyperpol/xixdnp_steady.m); unlike its T2e/T2n siblings, this file does not call `r1n_dnp`. The source file is [`xix_q_con_time_ensemble_r_T1n.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1n.m).

The no-argument function returns no MATLAB output. It plots the real part of the distance-averaged proton `Lz` expectation value against contact time and saves `xix_q_con_time_ensemble_r_T1n.fig` in the current directory. The plotted family differs only by nuclear T1; the four-node distance quadrature and remaining model/protocol settings stay fixed.
