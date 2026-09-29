# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1_r.m)

- Signature: `xix_q_con_time_ensemble_b1_r()`

## Purpose

Calculates steady-state proton signal versus XiX contact time with both electron B1/Rabi-frequency and electron–proton distance ensembles. Unlike `xix_q_con_time_ensemble_b1`, it varies distance as well as B1; unlike `xix_q_con_time_ensemble_r`, it does not use a single fixed B1 value.

## Model and settings

The function builds an electron–`1H` Q-band pair at `sys.magnet=1.2142` (1.2142 T), with electron Zeeman principal values `[2.00319 2.00319 2.00258]` and proton shift `[0 0 5]` (the source calls this a ppm guess). The Euler-angle triplets are `[0 10 0]` degrees, converted to radians. Spin temperature is 80 K. The basis is `sphten-liouv` with no approximation, propagator chop tolerance `1e-12`, and hygiene disabled.

The distance quadrature is `gaussleg(3.5,20,3)` (Å); the B1 quadrature is `gaussleg(10e6,20e6,5)` (Hz). For each distance node, the electron and proton coordinates are reset to `[0 0 0]` and `[0 0 r]`. Relaxation is `t1_t2`, diagonal-only, with `dibari` equilibrium; electron R1/R2 are 1000/200000 and proton R2 is 50000. Proton R1 is a distance- and orientation-dependent function handle calling `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r(n),bet)`. Rate-value units are not annotated in the source.

## Experiment and scan

At each distance and B1 node the function scans `nloops=1:64`. XiX uses a 48 ns pulse duration and π phase for the second pulse; total contact time is `2*nloops*48e-9` s (96 ns to 6.144 μs). It uses grid `rep_2ang_800pts_sph`, `addshift=-13e6`, `el_offs=61e6`, and shot spacing of 153 μs minus the total pulse duration.

## Calculation and output

The proton detector is `state(spin_system,'Lz','1H')`; each scan point is computed with `powder(spin_system,@xixdnp_steady,localpar,'esr')`. The code first averages over B1 using its quadrature weights, then averages over distance using weights multiplied by the radial (r^2) Jacobian and normalises by the weighted (r^2) sum. It plots the real proton (L_z) expectation value against total contact time in μs and saves `xix_q_con_time_ensemble_b1_r.fig`; there is no explicit MATLAB return value. The source comment estimates hours of calculation time, not a measured runtime here.
