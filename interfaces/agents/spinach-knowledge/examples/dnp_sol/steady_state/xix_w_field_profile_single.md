# examples/dnp_sol/steady_state/xix_w_field_profile_single.m

- MATLAB implementation: [examples/dnp_sol/steady_state/xix_w_field_profile_single.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_w_field_profile_single.m)

## Purpose

The source header estimates the calculation time as seconds.

`xix_w_field_profile_single()` calculates the steady-state XiX DNP proton response versus microwave resonance offset for one electron–proton spin system. It has no distance or microwave-field ensemble. The run still calls `powder` on `rep_2ang_800pts_sph`, so orientations are sampled with that configured grid.

## Model and sequence

The system is an electron (`E`) and `1H` at a W-band magnet setting of `3.4`. The electron Zeeman principal values are `[2.00319 2.00319 2.00258]`; the proton shift is the stated `[0 0 5]` ppm guess. The Euler angles are `[0 10 0]` and `[0 0 10]` degrees, converted to radians in the script. Spin temperature is `80`; the coordinates place the electron at the origin and proton at `[0 0 3.5]` (the coordinate unit is not stated).

The relaxation setup uses `t1_t2`, `r1n_dnp` for the distance/orientation-dependent electron–nuclear contribution, `r1_rates={1e3 r1n_rate}`, `r2_rates={200e3,50e3}`, diagonal relaxation retention, and the Di Bari equilibrium. The function selects the full `sphten-liouv` basis with no approximation, disables hygiene, and sets `prop_chop=1e-12`.

## Run parameters and output

The detected operator is proton `Lz`. The XiX parameters are electron nutation frequency `20e6 Hz`, pulse duration `18e-9 s`, `10` XiX DNP blocks, second-pulse phase `pi` (commented as inverted), shot spacing `167e-6`, and additive shift `-33e6`. The offset scan has 201 points from `-300e6` to `300e6` in the source; the plot expresses the horizontal axis in MHz.

`powder(spin_system,@xixdnp_steady,parameters,'esr')` supplies the steady-state calculation. The script plots the real proton `I_z` expectation against resonance offset and saves `xix_w_field_profile_single.fig`.

## Dependencies and scope

The entry point depends on Spinach system construction, basis, state, powder, and plotting helpers, plus `r1n_dnp` and the XiX sequence function `xixdnp_steady`. The file defines a single parameter set and offset profile; it does not compute distance- or B1-ensemble averages.
