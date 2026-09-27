# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r_T2n.m

- Signature: `xix_q_rep_time_ensemble_r_T2n()`

## Purpose

Simulate how proton transverse relaxation time affects steady-state XiX DNP as a function of repetition time, averaged over an electron–proton distance ensemble. The source notes a calculation time of hours.

## Physical / mathematical content

- A Q-band electron–proton system at 1.2142 T and 80 K uses a trityl electron g-tensor and a proton chemical-shift estimate. The relaxation model includes distance- and orientation-dependent proton longitudinal relaxation, with proton transverse relaxation set by each selected $T_{2n}$.
- The distance ensemble is integrated using Gauss–Legendre weights and an $r^2$ Jacobian. The plotted observable is the negative real part of the distance-averaged, steady-state proton polarization.

## Numerical / algorithmic content

- The calculation uses a `sphten-liouv` basis without approximation and evaluates `xixdnp_steady` through `powder(...,'esr')` on the `rep_2ang_800pts_sph` grid.
- For each of five $T_{2n}$ values (2000, 200, 20, 2, and 0.2 μs), it sweeps 30 logarithmically spaced repetition times from $10^{-5}$ to $10^{-3}$ s and sets shot spacing to the repetition time minus the duration of 36 two-pulse XiX blocks.

## Implementation structure

- The main function initializes the figure, calls `xix_rep_time_ensemble_r(T2n)` for each relaxation time, adds a legend, and saves `xix_q_rep_time_ensemble_r_T2n.fig`.
- The helper sets up the spin system and pulse parameters, obtains a three-point distance quadrature from 3.5 to 20 Å, runs the steady-state powder simulation at every distance and repetition time, averages over distance, and plots the resulting curve.
