# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T2n.m

- Signature: `xix_q_con_time_ensemble_r_T2n()`

## Purpose

Compares steady-state XiX proton-polarisation contact-time curves across five nuclear T2 values, averaging over electron–proton distance for each curve.

## Model and scan

The source varies proton T2 across 2000, 200, 20, 2 and 0.2 μs; proton T1 is set by the distance- and orientation-dependent `r1n_dnp` rate, and electron T1 is fixed at 1 ms. Each run models a trityl–proton pair at 1.2142 T and 80 K with three Gauss–Legendre distance nodes from 3.5 to 20 Å. Diagonal relaxation terms are retained and the equilibrium is `dibari`; the full spherical-tensor Liouville basis is used without basis approximation. Distance averaging includes the radial `r^2` Jacobian.

Each curve scans 1–64 XiX loops using 48 ns pulses, inverted second-pulse phase, 18 MHz electron nutation frequency and an 800-point two-angle spherical powder grid. The source sets −13 MHz added shift, +61 MHz electron offset, and 153 μs shot spacing less total pulse duration; steady states are evaluated with `powder(...,@xixdnp_steady,...,'esr')`. Distance-averaged real proton `Lz` expectation-value curves are overlaid, labelled by T2n, and saved as `xix_q_con_time_ensemble_r_T2n.fig`.
