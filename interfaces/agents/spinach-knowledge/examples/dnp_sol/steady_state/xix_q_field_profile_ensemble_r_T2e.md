# examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2e.m

[Source MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/dnp_sol/steady_state/xix_q_field_profile_ensemble_r_T2e.m) · Signature: `xix_q_field_profile_ensemble_r_T2e()`.

## Use and distinctive variant

This steady-state XiX DNP example scans the electron transverse relaxation time (T_{2e}), plotting one field profile per value after averaging over an electron–proton distance quadrature. The source estimates a run time of minutes.

## Setup and scan

The spin system is `{'E','1H'}` at `sys.magnet=1.2142` (source comment: Q-band) and `inter.temperature=80` (source comment: spin temperature). Its electron Zeeman values are `[2.00319 2.00319 2.00258]`; the proton values are `[0 0 5]`, described in the source as a ppm guess. The Euler arrays `[0 10 0]` and `[0 0 10]` are scaled by `pi/180`. It uses `sphten-liouv` with approximation `none`, `sys.tols.prop_chop=1e-12`, and disables `hygiene`.

Distances use `gaussleg(3.5,20,3)`; for each point the electron is at the origin and the proton coordinate is `[0 0 r]`, and the resulting profiles are integrated with the quadrature weights and (r^2) radial Jacobian. The scanned electron (T_{2e}) values are `[50e-6,15e-6,5e-6,1.5e-6,0.5e-6]` seconds. Relaxation is `t1_t2`, diagonal terms are retained, equilibrium is `dibari`, R1 assignments are `{1e3,r1n_rate}`, and R2 assignments are `{1/T2e,50e3}` (electron/proton order). Here `r1n_rate` calls `r1n_dnp` with the magnet, temperature, constants `2.00230`, `1e-3`, `52`, the sampled distance, and orientation angle `bet`.

For each distance, the XiX kernel is run on 201 offsets from `-100e6` to `100e6` Hz with grid `rep_2ang_800pts_sph`, `18e6` Hz electron nutation frequency, `48e-9` s pulses, 36 blocks, phase `pi`, and `addshift=-13e6`. Shot spacing is assigned as `204e-6 - 2*nloops*pulse_dur`.

## Dependencies, output, and limits

Dependencies are local helper `xix_field_profile_ensemble_r`, `r1n_dnp`, `gaussleg`, Spinach `create`, `basis`, `state`, and `powder`, plus `xixdnp_steady`. The plot shows the real proton (I_Z) response versus offset in MHz, uses vertical limits ([-3×10^{-3},3×10^{-3}]), and is saved as `xix_q_field_profile_ensemble_r_T2e.fig`; no numeric array is declared as a function output. The source explicitly gives seconds for T2e and pulse duration and Hz for offsets and nutation frequency. It does not state units for the magnet, spin-temperature, distance, relaxation-rate, or additional-shift values.
