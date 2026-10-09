# examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_b1_r.m

- Signature: `xix_q_rep_time_ensemble_b1_r()`
- Source: [xix_q_rep_time_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_rep_time_ensemble_b1_r.m)

## Purpose and protocol variant

Scans XiX steady-state proton polarisation against shot repetition time while averaging over both electron–proton distance and microwave B1 quadratures. This is the combined-ensemble, repetition-time variant; the source estimates hours.

## Inputs and scan axes

- Four Gauss–Legendre distance nodes over 3.5–20 Å and six B1 nodes over `10e6`–`20e6` Hz (orders 3 and 5 yield one extra node each).
- Thirty logarithmically spaced repetition times, `logspace(-5,-3,30)`; shot spacing is `rep_time-2*nloops*pulse_dur` with pulse duration explicitly set in seconds.
- Fixed `parameters.el_offs=-39e6` and `parameters.addshift=-13e6`; no offset scan.

## Shared spin-system setup

The source sets `sys.magnet=1.2142` for its Q-band setup, uses isotopes `E` and `1H`, electron Zeeman eigenvalues `[2.00319 2.00319 2.00258]`, and proton eigenvalues `[0 0 5]` (the source calls these a ppm guess). Euler angles are `(pi/180)*{[0 10 0],[0 0 10]}`. Spin temperature is set to `80`. The basis is `sphten-liouv` with no approximation, `sys.tols.prop_chop=1e-12`, and `sys.disable={'hygiene'}`. Relaxation uses `t1_t2`, a distance-/orientation-dependent `r1n_dnp(...)` handle with `inter.r1_rates={1e3,r1n_rate}`, `inter.r2_rates={200e3,50e3}`, diagonal relaxation retention, and `dibari` equilibrium. The source does not annotate units beside the temperature or these rate values.

All use the `rep_2ang_800pts_sph` powder grid, detect proton `Lz`, and call `powder(spin_system,@xixdnp_steady,parameters,'esr')`. The XiX block settings are 48 ns pulse duration, 36 blocks, and phase `pi` (the second pulse is inverted); `addshift=-13e6`. Values are reproduced as configured; where the source does not label a unit, none is added.

## Helpers and source links

The run relies on Spinach `create`, `basis`, `state`, and `powder`, plus `gaussleg` for quadrature and `r1n_dnp` for the relaxation handle. The steady-state protocol is provided by `xixdnp_steady`. Plotting uses the repository `kfigure`/`kxlabel`/`kylabel`/`kgrid` helpers. See the linked source and the repository's XiX DNP reference, [XiX DNP publication (repository DOI link)](https://doi.org/10.1021/jacs.1c09900) (listed in `examples/dnp_sol/xix_dnp/xix_paper.url`).

## Output and limits

For each distance and B1 node the function builds the system, evaluates the 30 repetition-time points (the source uses `parfor` over those points), averages over B1 weights, then averages over distance with the radial `r^2` Jacobian. The plotted value is the real proton `I_z` expectation versus repetition time in ms; it saves `xix_q_rep_time_ensemble_b1_r.fig`. The script produces no separate numeric table.

## Clarification

This is the only one of these four repetition-time variants that combines the four-node 3.5–20 Å distance quadrature with the six-node 10–20 MHz B1 quadrature; it evaluates 30 repetition times for each of 24 distance/B1 pairs.
