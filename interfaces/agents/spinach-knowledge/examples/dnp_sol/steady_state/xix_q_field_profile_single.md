# examples/dnp_sol/steady_state/xix_q_field_profile_single.m

- Signature: `xix_q_field_profile_single()`

## Purpose

Calculates the steady-state proton signal across a microwave resonance-offset sweep for a single XiX DNP spin system (seconds, according to the source comment).

## Physical and numerical setup

The system contains an electron and a proton at a fixed 3.5 Å separation, at 80 K and a Q-band field of 1.2142 T. The electron Zeeman tensor is set to the trityl values in the source; the proton shift is specified as 5 ppm. The electron nutation frequency is fixed at 18 MHz. Relaxation uses the distance- and orientation-dependent `r1n_dnp` rate function, with R1 rates `{1e3, r1n_rate}` and R2 rates `{200e3, 50e3}`.

## Calculation and output

The script uses the full spherical-tensor Liouville basis and the `rep_2ang_800pts_sph` powder grid. It runs `xixdnp_steady` through `powder(...,'esr')` for 36 XiX blocks, with a 48 ns pulse and 204 μs shot repetition time. The 201 microwave offsets span −100 to +100 MHz. It plots the real proton (I_z) expectation value against offset and saves `xix_q_field_profile_single.fig`.
