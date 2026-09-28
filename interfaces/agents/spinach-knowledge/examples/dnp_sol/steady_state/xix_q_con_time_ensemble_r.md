# examples/dnp_sol/steady_state/xix_q_con_time_ensemble_r.m

- Signature: `xix_q_con_time_ensemble_r()`

## Purpose

Calculates steady-state proton polarisation versus XiX contact time, averaged over an electron–proton distance distribution.

## Model and distance average

The source uses a trityl electron–proton pair at 1.2142 T and 80 K. Three Gauss–Legendre distance nodes span 3.5–20 Å. For each distance it updates the pair coordinates and the orientation-dependent proton T1 rate via `r1n_dnp`; the model also specifies T2 rates, diagonal relaxation retention and `dibari` equilibrium. The spin system uses the full spherical-tensor Liouville basis with no basis approximation. The distance average applies the quadrature weights multiplied by the radial `r^2` Jacobian.

## XiX scan and output

The experiment detects proton `Lz` on an 800-point two-angle spherical powder grid. It computes steady state for 1–64 XiX loops, using 48 ns pulses and an inverted second-pulse phase; the contact time is twice the loop count times the pulse duration. The source sets an 18 MHz electron nutation frequency, −13 MHz added shift, +61 MHz electron offset, and 153 μs shot spacing less total pulse duration. Each case uses `powder(...,@xixdnp_steady,...,'esr')`. The distance-averaged real proton expectation value is plotted against total contact time and saved as `xix_q_con_time_ensemble_r.fig`.
