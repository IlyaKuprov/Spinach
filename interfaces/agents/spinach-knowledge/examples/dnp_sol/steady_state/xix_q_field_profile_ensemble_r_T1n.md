# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1n.m

[Source MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T1n.m) · Signature: `xix_q_field_profile_ensemble_r_T1n()`.

## Use and distinctive variant

This example compares steady-state XiX DNP field profiles while sweeping the proton longitudinal relaxation time (T_{1n}). Each curve is averaged over an electron–proton distance quadrature; the source estimates a run time of minutes.

## Setup and scan

The model contains `E` and `1H`, with `sys.magnet=1.2142` (commented as Q-band) and `inter.temperature=80` (commented as spin temperature). The electron Zeeman values are `[2.00319 2.00319 2.00258]`; the proton values are `[0 0 5]` (the source calls this a ppm guess). Euler-angle arrays `[0 10 0]` and `[0 0 10]` are multiplied by `pi/180`. The spin system uses `sphten-liouv`, approximation `none`, propagation chop tolerance `1e-12`, and disables `hygiene`.

For each point of the four-point distance quadrature, the electron is at the origin and the proton is placed at `[0 0 r]` `gaussleg(3.5,20,3)`; each distance result is combined using its quadrature weight and the radial (r^2) Jacobian. The scan is `T1n=[50 5 0.5 0.05 0.005] seconds`. It sets the relaxation model to `t1_t2`, keeps diagonal relaxation, uses `dibari` equilibrium, R1 assignments `{1e3,1/T1n}`, and fixed R2 assignments `{200e3,50e3}` in electron/proton order.

Each profile uses 201 microwave offsets from `-100e6` to `100e6` Hz, an `18e6` Hz electron nutation frequency, a `48e-9` s pulse, 36 XiX blocks, phase `pi`, and `addshift=-13e6`. Shot spacing is set by `204e-6 - 2*nloops*pulse_dur`. The orientation grid is `rep_2ang_800pts_sph`.

## Dependencies, output, and limits

The example calls local helper `xix_field_profile_ensemble_r`, Spinach's `gaussleg`, `create`, `basis`, `state`, and `powder`, and steady-state kernel `xixdnp_steady`. It plots the real proton (I_Z) response against offset in MHz, with a fixed vertical range of ([-1.3×10^{-3},1.3×10^{-3}]), and saves `xix_q_field_profile_ensemble_r_T1n.fig`. The function has no declared data output. The source labels T1n, pulse duration, and offset units, but does not give units for the magnet value, spin-temperature value, distance coordinates, relaxation-rate values, or `addshift`; those numbers are therefore left unitless here.
