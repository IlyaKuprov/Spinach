# examples/dnp_sol/beam_dnp/beam_field_profile.m

- Signature: `beam_field_profile()`

## Purpose

Computes the proton `I_z` signal at the end of a fixed-duration BEAM DNP contact as a function of microwave resonance offset. The associated experiment is described in [Science Advances](https://doi.org/10.1126/sciadv.abq0536); the source notes that the large powder grid makes this a minutes-scale calculation.

## Model and calculation

The X-band model (0.3483 T) contains an electron and two protons, with trityl g-tensor values, proton Zeeman estimates, Cartesian coordinates, and spin temperature 80 K. It uses a full Zeeman-Hilbert basis, proton `Lz` detection, and the BEAM pulse sequence (32 MHz electron nutation frequency, 20.0/28.7 ns pulses, 165 blocks). For each of 120 offsets spanning −60 to +60 MHz, the script adds a −3.3 MHz reference offset, runs the powder calculation on `rep_2ang_800pts_sph`, and records the last contact-curve point. The offset sweep is parallelised with `parfor`.
