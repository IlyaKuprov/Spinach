# examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_b1_r.m)

## Purpose

The source header estimates the calculation time as hours.

`xix_w_pulse_dur_ensemble_b1_r()` computes a steady-state XiX DNP proton response map over microwave offset and pulse duration, including both electron–proton distance and electron B1 quadratures. It is the combined-ensemble counterpart of the B1-only and distance-only scans.

## Model and sequence

The two-spin model uses magnet setting `3.4`, isotopes `E` and `1H`, electron Zeeman values `[2.00319 2.00319 2.00258]`, proton shift guess `[0 0 5]` ppm, Euler angles `[0 10 0]` and `[0 0 10]` degrees, and temperature `80`. It uses the full `sphten-liouv` basis, `prop_chop=1e-12`, `t1_t2` relaxation, `r1n_dnp` with the sampled distance, R1 entries `1e3`, R2 entries `200e3` and `50e3`, diagonal retention, and Di Bari equilibrium. As in the source, units are not assigned to the magnet, temperature, distance-quadrature arguments, relaxation entries, or additive shift.

The pulse sequence uses second-pulse phase `pi`, additive shift `-33e6`, and shot spacing `167e-6 - 360e-9`. At each pulse duration it sets `nloops=round(360e-9/(2*pulse_dur))`. The offset scan has 101 points from `-230e6` to `205e6 Hz` (plotted in MHz); the duration scan has 200 points from `2e-9` to `21e-9 s` (plotted in ns).

## Quadratures and output

Distance nodes and weights come from `[r,w]=gaussleg(3.5,20,3)`; the source labels this a distance ensemble but does not state the coordinate unit. At each node the script updates the proton coordinate and evaluates the distance-dependent relaxation rate. The six-node B1 quadrature is `[b1,wb1]=gaussleg(10e6,20e6,5)`, labelled Hz. It loops over distance and B1 nodes, parallelises the duration sweep with `parfor`, and calls `powder(...,@xixdnp_steady,...,'esr')` using `rep_2ang_800pts_sph`.

The accumulated map is averaged over B1 with `wb1`, then over distance with weights `r.^2 .* w`; the `r.^2` factor is the Jacobian noted in the source. The real proton `I_z` expectation is plotted against offset and pulse duration and saved as `xix_w_pulse_dur_ensemble_b1_r.fig`.

## Dependencies and scope

Dependencies are Spinach system, basis, state, powder, and plotting functions, `gaussleg`, `r1n_dnp`, and `xixdnp_steady`. The result includes the two stated quadratures and the configured powder orientation grid; it does not describe a measured experimental profile.
