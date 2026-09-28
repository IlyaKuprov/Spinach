# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r_T1e.m

- Signature: `xix_q_con_time_ensemble_r_T1e()`

## Purpose

Compares steady-state XiX proton-polarisation contact-time curves across five electron T1 values while averaging over electron–proton distance.

## Model and scan

The source repeats the calculation for T1e = 10, 3, 1, 0.3 and 0.1 ms. Each run models a trityl electron and proton at 1.2142 T and 80 K, with three Gauss–Legendre distance nodes from 3.5 to 20 Å. The distance-dependent proton T1 rate is evaluated with `r1n_dnp`; electron T1 is set to the selected value, and T2 rates, diagonal relaxation retention and `dibari` equilibrium are specified. The full spherical-tensor Liouville basis is used without basis approximation. Distance averaging includes the radial `r^2` Jacobian.

Each curve scans 1–64 XiX loops with 48 ns pulses, inverted second-pulse phase, 18 MHz electron nutation frequency, and an 800-point two-angle spherical powder grid. The source uses −13 MHz added shift, +61 MHz electron offset, 153 μs shot spacing minus total pulse duration, and computes each steady state with `powder(...,@xixdnp_steady,...,'esr')`. It overlays the distance-averaged real proton `Lz` expectation-value curves against total contact time, labels them by T1e, and saves `xix_q_con_time_ensemble_r_T1e.fig`.
