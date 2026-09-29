# examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_r.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_r.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_w_pulse_dur_ensemble_r.m)

## Purpose

The source header estimates the calculation time as hours.

`xix_w_pulse_dur_ensemble_r()` maps the steady-state XiX DNP proton response over microwave offset and pulse duration, averaging over an electron–proton distance quadrature. The electron nutation frequency is fixed at `20e6 Hz`; there is no B1 ensemble.

## Model and sequence

The two-spin model has magnet setting `3.4`, isotopes `E` and `1H`, electron Zeeman values `[2.00319 2.00319 2.00258]`, proton shift guess `[0 0 5]` ppm, Euler angles `[0 10 0]` and `[0 0 10]` degrees, and temperature `80`. It uses the full `sphten-liouv` basis, `prop_chop=1e-12`, `t1_t2` relaxation with distance-dependent `r1n_dnp`, R1 entries `1e3`, R2 entries `200e3` and `50e3`, diagonal relaxation retention, and Di Bari equilibrium. Units not assigned in the source are left as raw parameter values.

The sequence uses phase `pi` for the second pulse, additive shift `-33e6`, and shot spacing `167e-6 - 360e-9`. For each pulse duration, the XiX block count is `round(360e-9/(2*pulse_dur))`. The offset scan is 101 points from `-230e6` to `205e6 Hz`; the duration scan is 200 points from `2e-9` to `21e-9 s`. The figure displays offset in MHz and duration in ns.

## Distance quadrature and output

The four-node distance quadrature is specified as `[r,w]=gaussleg(3.5,20,3)` (the source does not annotate the coordinate unit). For each node, the script resets the proton coordinate, constructs the distance-dependent relaxation rate, creates the Spinach system, and runs the pulse-duration sweep with `parfor`. Each point calls `powder(...,@xixdnp_steady,...,'esr')` with `irr_powers=20e6` and orientation grid `rep_2ang_800pts_sph`.

The final map is averaged over distance using the quadrature weights multiplied by `r.^2`, the Jacobian identified in the source. It plots the real proton `I_z` expectation against offset and pulse duration and saves `xix_w_pulse_dur_ensemble_r.fig`.

## Dependencies and scope

The entry point uses Spinach system, basis, state, powder, and plotting routines, along with `gaussleg`, `r1n_dnp`, and `xixdnp_steady`. It averages distance only; B1 is held at the single stated value.
