# examples/dnp_sol/steady_state/xix_q_field_profile_single.m

[Source MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_field_profile_single.m) · Signature: `xix_q_field_profile_single()`.

## Use and distinctive protocol

This is the compact steady-state XiX DNP field-profile example for one electron–proton spin system and one fixed separation, without distance-ensemble averaging. The source estimates a run time of seconds. Although the distance is fixed, the calculation still calls `powder` with an orientation grid.

## Setup and scan

The model is `{'E','1H'}`, with `sys.magnet=1.2142` (source comment: Q-band) and `inter.temperature=80` (source comment: spin temperature). The electron Zeeman values are `[2.00319 2.00319 2.00258]`; proton values are `[0 0 5]` (the source calls these a ppm guess). Euler arrays `[0 10 0]` and `[0 0 10]` are multiplied by `pi/180`. The coordinates are `[0 0 0]` and `[0 0 3.500]`; the script derives `r_en` from the proton z coordinate and passes it to `r1n_dnp` together with orientation angle `bet` and the constants `2.00230`, `1e-3`, and `52`.

Relaxation uses `t1_t2`, R1 assignments `{1e3,r1n_rate}`, R2 assignments `{200e3,50e3}`, diagonal retention, and `dibari` equilibrium. The basis is `sphten-liouv` with approximation `none`; it disables `hygiene` and sets propagation chop tolerance `1e-12`. The XiX setup uses grid `rep_2ang_800pts_sph`, 201 offsets from `-100e6` to `100e6` Hz, `18e6` Hz electron nutation frequency, `48e-9` s pulses, 36 blocks, phase `pi`, and `addshift=-13e6`. Shot spacing is computed as `204e-6 - 2*nloops*pulse_dur`.

## Dependencies, output, and limits

The example calls `r1n_dnp`, Spinach `create`, `basis`, `state`, and `powder`, and kernel `xixdnp_steady`. It plots the real proton (I_Z) expectation value against offset in MHz with padded vertical limits, then saves `xix_q_field_profile_single.fig`. The function declares no numeric output. This source explicitly labels offsets and electron nutation frequency in Hz and pulse duration in seconds; units for `sys.magnet`, `inter.temperature`, the coordinate 3.500, relaxation-rate values, and `addshift` are not stated.
