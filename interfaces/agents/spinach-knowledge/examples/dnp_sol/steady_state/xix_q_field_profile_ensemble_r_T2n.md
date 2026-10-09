# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2n.m

[Source MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2n.m) · Signature: `xix_q_field_profile_ensemble_r_T2n()`.

## Use and distinctive variant

This steady-state XiX DNP example varies the proton transverse relaxation time (T_{2n}), averaging each field profile over the electron–proton distance quadrature. The source estimates a run time of minutes.

## Setup and scan

The two-spin model uses `{'E','1H'}`, `sys.magnet=1.2142` (commented as Q-band), and `inter.temperature=80` (commented as spin temperature). Electron Zeeman values are `[2.00319 2.00319 2.00258]`; the proton values are `[0 0 5]` (source: ppm guess). Euler arrays `[0 10 0]` and `[0 0 10]` are scaled by `pi/180`. The basis is `sphten-liouv` with approximation `none`; propagation chop tolerance is `1e-12`, and `hygiene` is disabled.

The four-point distance quadrature is `gaussleg(3.5,20,3)`; for each point the electron is at the origin and the proton coordinate is `[0 0 r]` and radial averaging weighted by quadrature weight times (r^2). The scanned proton (T_{2n}) values are `[2e-3,200e-6,20e-6,2e-6,0.2e-6]` seconds. It uses `t1_t2` relaxation, diagonal relaxation terms, `dibari` equilibrium, R1 assignments `{1e3,r1n_rate}`, and R2 assignments `{200e3,1/T2n}` (electron/proton order). The proton R1 function is `r1n_dnp(sys.magnet,inter.temperature,2.00230,1e-3,52,r,bet)`.

The steady-state XiX calculation uses 201 offsets from `-100e6` to `100e6` Hz, grid `rep_2ang_800pts_sph`, `18e6` Hz electron nutation frequency, `48e-9` s pulses, 36 blocks, phase `pi`, and `addshift=-13e6`. Shot spacing is set by `204e-6 - 2*nloops*pulse_dur`.

## Dependencies, output, and limits

The example depends on `xix_field_profile_ensemble_r`, `r1n_dnp`, `gaussleg`, Spinach `create`, `basis`, `state`, and `powder`, and kernel `xixdnp_steady`. It plots the real proton (I_Z) response versus offset in MHz with vertical limits ([-3×10^{-3},3×10^{-3}]), and saves `xix_q_field_profile_ensemble_r_T2n.fig`. No numerical result is declared as a function output. Units for the magnet value, spin-temperature value, quadrature coordinates, relaxation-rate values, and `addshift` are not specified in this source.
