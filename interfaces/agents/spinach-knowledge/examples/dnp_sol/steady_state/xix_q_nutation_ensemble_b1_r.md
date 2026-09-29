# examples/dnp_sol/steady_state/xix_q_nutation_ensemble_b1_r.m

- Signature: `xix_q_nutation_ensemble_b1_r()`
- Source: [xix_q_nutation_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_nutation_ensemble_b1_r.m)

## Purpose and protocol variant

Builds six steady-state XiX DNP field profiles against microwave offset, one for each electron nutation frequency and its paired shot repetition time. It averages over both the electron–proton distance and microwave B1 quadratures; B1 quadrature limits scale with each selected nutation frequency. The source estimates minutes.

## Inputs and scan axes

- Nutation frequencies `nu=1e6*[6.8 9.6 13.5 17.5 25 36]` Hz, paired in order with `srt=1e-3*[0.051 0.051 0.102 0.153 0.153 0.306]` s.
- Four Gauss–Legendre distance nodes on 3.5–20 Å; six B1 nodes from `0.2*nu` to `1.2*nu` Hz for each frequency. The distance average includes the radial `r^2` Jacobian.
- Thirteen microwave resonance offsets from `linspace(-64e6,-52e6,13)` Hz. The script sets `parameters.addshift=-13e6` and derives shot spacing as `srt-2*nloops*pulse_dur`.

## Shared spin-system setup

The source sets `sys.magnet=1.2142` for its Q-band setup, uses isotopes `E` and `1H`, electron Zeeman eigenvalues `[2.00319 2.00319 2.00258]`, and proton eigenvalues `[0 0 5]` (the source calls these a ppm guess). Euler angles are `(pi/180)*{[0 10 0],[0 0 10]}`. Spin temperature is set to `80`. The basis is `sphten-liouv` with no approximation, `sys.tols.prop_chop=1e-12`, and `sys.disable={'hygiene'}`. Relaxation uses `t1_t2`, a distance-/orientation-dependent `r1n_dnp(...)` handle with `inter.r1_rates={1e3,r1n_rate}`, `inter.r2_rates={200e3,50e3}`, diagonal relaxation retention, and `dibari` equilibrium. The source does not annotate units beside the temperature or these rate values.

All use the `rep_2ang_800pts_sph` powder grid, detect proton `Lz`, and call `powder(spin_system,@xixdnp_steady,parameters,'esr')`. The XiX block settings are 48 ns pulse duration, 36 blocks, and phase `pi` (the second pulse is inverted); `addshift=-13e6`. Values are reproduced as configured; where the source does not label a unit, none is added.

## Helpers and source links

The run relies on Spinach `create`, `basis`, `state`, and `powder`, plus `gaussleg` for quadrature and `r1n_dnp` for the relaxation handle. The steady-state protocol is provided by `xixdnp_steady`. Plotting uses the repository `kfigure`/`kxlabel`/`kylabel`/`kgrid` helpers. See the linked source and the repository's XiX DNP reference, [XiX DNP publication (repository DOI link)](https://doi.org/10.1021/jacs.1c09900) (listed in `examples/dnp_sol/xix_dnp/xix_paper.url`).

## Output and limits

The local function `xix_field_profile_b1_r(nu,srt)` constructs the distance-specific system, runs the powder steady state over B1 nodes and the 13 offsets, then performs B1-weighted and radial distance-weighted averaging. It plots `-real(dnp)` in a 3-D figure (nutation frequency in MHz, resonance offset in MHz, and the steady-state proton `<I_Z>` label) and saves `xix_q_nutation_ensemble_b1_r.fig`. The source saves the figure, not a separate numeric result table.

## Clarification

Despite “ensemble” in the name, the six nutation frequencies are a discrete outer profile scan, each paired with a different `srt`; the B1 quadrature is six nodes scaled separately by each `nu`, not a single fixed B1 interval.
