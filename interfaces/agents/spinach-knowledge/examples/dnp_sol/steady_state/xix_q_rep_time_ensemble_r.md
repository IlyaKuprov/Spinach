# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r.m

- Signature: `xix_q_rep_time_ensemble_r()`
- Source: [xix_q_rep_time_ensemble_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_r.m)

## Purpose and protocol variant

Scans XiX steady-state proton polarisation against shot repetition time, averaging over electron–proton distance while keeping microwave nutation power fixed. The source estimates minutes.

## Inputs and scan axes

- Four Gauss–Legendre distance nodes over 3.5–20 Å; the distance average includes the radial `r^2` Jacobian.
- Thirty logarithmically spaced repetition times, `logspace(-5,-3,30)`; each shot spacing subtracts the two-pulse-per-block XiX train duration from the repetition time. Pulse duration is 48 ns.
- Fixed electron nutation frequency `parameters.irr_powers=18e6` Hz, fixed `parameters.el_offs=-39e6`, and `parameters.addshift=-13e6`.

## Shared spin-system setup

The source sets `sys.magnet=1.2142` for its Q-band setup, uses isotopes `E` and `1H`, electron Zeeman eigenvalues `[2.00319 2.00319 2.00258]`, and proton eigenvalues `[0 0 5]` (the source calls these a ppm guess). Euler angles are `(pi/180)*{[0 10 0],[0 0 10]}`. Spin temperature is set to `80`. The basis is `sphten-liouv` with no approximation, `sys.tols.prop_chop=1e-12`, and `sys.disable={'hygiene'}`. Relaxation uses `t1_t2`, a distance-/orientation-dependent `r1n_dnp(...)` handle with `inter.r1_rates={1e3,r1n_rate}`, `inter.r2_rates={200e3,50e3}`, diagonal relaxation retention, and `dibari` equilibrium. The source does not annotate units beside the temperature or these rate values.

All use the `rep_2ang_800pts_sph` powder grid, detect proton `Lz`, and call `powder(spin_system,@xixdnp_steady,parameters,'esr')`. The XiX block settings are 48 ns pulse duration, 36 blocks, and phase `pi` (the second pulse is inverted); `addshift=-13e6`. Values are reproduced as configured; where the source does not label a unit, none is added.

## Helpers and source links

The run relies on Spinach `create`, `basis`, `state`, and `powder`, plus `gaussleg` for quadrature and `r1n_dnp` for the relaxation handle. The steady-state protocol is provided by `xixdnp_steady`. Plotting uses the repository `kfigure`/`kxlabel`/`kylabel`/`kgrid` helpers. See the linked source and the repository's XiX DNP reference, [XiX DNP publication (repository DOI link)](https://doi.org/10.1021/jacs.1c09900) (listed in `examples/dnp_sol/xix_dnp/xix_paper.url`).

## Output and limits

For each distance node, the code constructs the spin system and runs the 30-point repetition-time scan at fixed B1, then averages the proton `I_z` expectation with the radial distance weights. It plots the real expectation versus repetition time in ms and saves `xix_q_rep_time_ensemble_r.fig`. No separate numeric table is written.

## Clarification

Unlike the B1-ensemble repetition-time variants, this file has no B1 quadrature: it uses one fixed 18 MHz electron nutation frequency and only quadratures the distance.
