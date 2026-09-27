# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1n.m

- Signature: `xix_q_con_time_ensemble_r_T1n()`

## Purpose

Compares steady-state XiX proton-polarisation contact-time curves across five nuclear T1 values, with an electron–proton distance average in each curve.

## Model and scan

The source varies proton T1 across 50, 5, 0.5, 0.05 and 0.005 s; electron T1 is fixed at 1 ms. Each run models a trityl–proton pair at 1.2142 T and 80 K, with three Gauss–Legendre distance nodes from 3.5 to 20 Å. The proton T1 rate depends on distance and orientation through `r1n_dnp`; T2 rates, diagonal relaxation retention and `dibari` equilibrium are set explicitly. It uses the full spherical-tensor Liouville basis without basis approximation, and applies distance weights with the radial `r^2` Jacobian.

Each curve scans 1–64 XiX loops with 48 ns pulses, inverted second-pulse phase, 18 MHz electron nutation frequency and an 800-point two-angle spherical powder grid. The source sets −13 MHz added shift, +61 MHz electron offset, and 153 μs shot spacing less total pulse duration; steady states are calculated with `powder(...,@xixdnp_steady,...,'esr')`. The distance-averaged real proton `Lz` expectation-value curves are overlaid, labelled by proton T1, and saved as `xix_q_con_time_ensemble_r_T1n.fig`.
