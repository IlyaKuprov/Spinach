# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1n.m

- Signature: `xix_q_field_profile_ensemble_r_T1n()`

## Purpose

Simulate steady-state XiX DNP field profiles at Q-band for several proton longitudinal relaxation times, averaging over an electron–proton distance ensemble. The source estimates a calculation time of minutes.

## Physical / mathematical content

The spin system contains an electron and a proton at a 1.2142 T field and 80 K. Electron and proton Zeeman interactions use the specified g-tensor and proton chemical-shift values. Electron–proton separation varies across a three-point Gauss–Legendre ensemble from 3.5 to 20; each result is weighted by the quadrature weight and the radial Jacobian, \(r^2\). The five proton \(T_{1n}\) values are 50, 5, 0.5, 0.05, and 0.005 seconds.

## Numerical / algorithmic content

For each relaxation time and distance, the helper constructs a Spinach system in an unrestricted spherical-tensor Liouville basis with `t1_t2` relaxation, then detects proton `Lz`. It evaluates `xixdnp_steady` through `powder(...,'esr')` at 201 microwave resonance offsets from −100 to 100 MHz, using the specified powder grid and XiX pulse parameters. After radial averaging, it plots the real steady-state proton polarisation against offset.

## Implementation structure

The entry-point function prepares the figure, calls `xix_field_profile_ensemble_r(T1n)` for each relaxation time, adds a legend, and saves `xix_q_field_profile_ensemble_r_T1n.fig`. The local helper sets the spin system and experiment parameters, runs the distance-resolved steady-state simulations, performs the weighted average, and plots each curve.
