# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_b1.m

- Signature: `xix_q_rep_time_ensemble_b1()`
- Source: [xix_q_rep_time_ensemble_b1.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_b1.m)

## Purpose and protocol variant

Scans XiX steady-state proton polarisation against shot repetition time, averaging over a microwave B1 ensemble only. Electron–proton distance and microwave offset are fixed. The source estimates minutes.

## Inputs and scan axes

- Six Gauss–Legendre B1 nodes spanning `10e6`–`20e6` Hz.
- Thirty logarithmically spaced repetition times, `logspace(-5,-3,30)`. They are used in the shot-spacing subtraction with the source's 48 ns pulse duration, i.e. on the same time scale (seconds). For every point, `shot_spacing=rep_time-2*nloops*pulse_dur`.
- Fixed Cartesian coordinate values at z = 0 and z = 3.500 (no unit is annotated beside these coordinates in this source); fixed `parameters.el_offs=-39e6` and `parameters.addshift=-13e6`.

## Shared spin-system setup

The source sets `sys.magnet=1.2142` for its Q-band setup, uses isotopes `E` and `1H`, electron Zeeman eigenvalues `[2.00319 2.00319 2.00258]`, and proton eigenvalues `[0 0 5]` (the source calls these a ppm guess). Euler angles are `(pi/180)*{[0 10 0],[0 0 10]}`. Spin temperature is set to `80`. The basis is `sphten-liouv` with no approximation, `sys.tols.prop_chop=1e-12`, and `sys.disable={'hygiene'}`. Relaxation uses `t1_t2`, a distance-/orientation-dependent `r1n_dnp(...)` handle with `inter.r1_rates={1e3,r1n_rate}`, `inter.r2_rates={200e3,50e3}`, diagonal relaxation retention, and `dibari` equilibrium. The source does not annotate units beside the temperature or these rate values.

All use the `rep_2ang_800pts_sph` powder grid, detect proton `Lz`, and call `powder(spin_system,@xixdnp_steady,parameters,'esr')`. The XiX block settings are 48 ns pulse duration, 36 blocks, and phase `pi` (the second pulse is inverted); `addshift=-13e6`. Values are reproduced as configured; where the source does not label a unit, none is added.

## Helpers and source links

The run relies on Spinach `create`, `basis`, `state`, and `powder`, plus `gaussleg` for quadrature and `r1n_dnp` for the relaxation handle. The steady-state protocol is provided by `xixdnp_steady`. Plotting uses the repository `kfigure`/`kxlabel`/`kylabel`/`kgrid` helpers. See the linked source and the repository's XiX DNP reference, [XiX DNP publication (repository DOI link)](https://doi.org/10.1021/jacs.1c09900) (listed in `examples/dnp_sol/xix_dnp/xix_paper.url`).

## Output and limits

For each B1 node, the script sets `parameters.irr_powers=b1(k)`, evaluates all repetition times with `powder`, and averages over B1 weights. It plots the real proton `I_z` expectation versus repetition time in ms, then saves `xix_q_rep_time_ensemble_b1.fig`. It does not scan distance or offset and does not export a separate numeric table.

## Clarification

The fixed 3.500 Å separation comes directly from the coordinates; the “ensemble” here is B1 alone. The script holds the offset field at `-39e6` rather than sweeping offsets.
