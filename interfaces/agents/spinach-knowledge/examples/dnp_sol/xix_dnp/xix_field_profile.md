# examples/dnp_sol/xix_dnp/xix_field_profile.m

- Signature: `xix_field_profile()`

## Purpose

Profiles the proton (I_z) signal at the end of a fixed XiX DNP contact as microwave resonance offset is varied. Further information: https://doi.org/10.1021/jacs.1c09900. Calculation time: minutes for the large powder grid. In this implementation the electron nutation frequency is fixed at 17.8 MHz; only offset is swept.

## Physical / mathematical content

The spin system is a trityl electron and two protons with anisotropic Zeeman interactions, specified coordinates, and spin temperature 80 K. The signal is the real final point of the calculated proton longitudinal-polarization contact curve.

## Numerical / algorithmic content

The script evaluates 120 offsets from -150 to +150 MHz, adding a -13 MHz reference point to the simulated offset. Each point uses 150 XiX blocks of 48 ns pulses and the `rep_2ang_1600pts_sph` powder grid; independent offsets are evaluated with `parfor`.

## Implementation structure

After constructing the full Zeeman–Hilbert basis and proton (L_z) detector, it calls `powder` with `@xixdnp` for each offset and plots the last contact-curve value against the unshifted offset axis.
