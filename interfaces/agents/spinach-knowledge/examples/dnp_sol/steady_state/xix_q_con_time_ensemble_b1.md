# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_con_time_ensemble_b1.m)

- Signature: `xix_q_con_time_ensemble_b1()`

## Purpose

Calculates the steady-state proton signal versus XiX contact time while averaging over an electron microwave-field (B1/Rabi-frequency) ensemble. This variant holds the electron–proton separation fixed at 3.5 Å; it does not scan distance.

## Model and settings

The function builds an electron–`1H` pair for the Q-band setting `sys.magnet=1.2142` (1.2142 T). The electron Zeeman principal values are `[2.00319 2.00319 2.00258]`; the proton shift is `[0 0 5]`, described in the source as a ppm guess. Both Euler-angle triplets are `[0 10 0]` degrees, converted to radians. The pair coordinates are `[0 0 0]` and `[0 0 3.5]` Å. Spin temperature is 80 K.

The basis is `sphten-liouv` with `approximation='none'`; the propagator chop tolerance is `1e-12`, and hygiene is disabled. Relaxation uses `t1_t2`, diagonal relaxation terms, and the `dibari` equilibrium. The electron rates are set to 1000 (R1) and 200000 (R2); proton R2 is 50000. Proton R1 is supplied as a function handle calling `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r_en,bet)`, so it depends on the helper's distance and orientation inputs. The source does not annotate units for these rate values.

## Experiment and scan

`gaussleg(10e6,20e6,5)` supplies six quadrature nodes and weights for the B1 ensemble (the source labels the range in Hz). The electron nutation frequency is set to each node. For every field node, the function scans `nloops=1:64`; each loop contributes two 48 ns pulses, giving contact times `2*nloops*48e-9` s (96 ns to 6.144 μs). The second pulse phase is π. The powder grid is `rep_2ang_800pts_sph`; the source also sets `addshift=-13e6`, `el_offs=61e6`, and shot spacing to 153 μs minus the total pulse duration.

## Calculation and output

The proton detector is `state(spin_system,'Lz','1H')`. Each field/contact-time point is evaluated by `powder(spin_system,@xixdnp_steady,localpar,'esr')`; the resulting signal is reduced over B1 using the quadrature weights. The function plots the real proton (L_z) expectation value against total contact time in μs and saves `xix_q_con_time_ensemble_b1.fig`. It returns no explicit MATLAB output. The source comment estimates calculation time as hours; this is not a measured runtime here.
