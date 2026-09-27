# examples/dnp_sol/tppm_dnp/tppm_field_profile.m

- Signature: `tppm_field_profile()`

## Purpose

Calculates the final proton (I_z) value of a fixed-contact-time TPPM DNP sequence as a function of microwave resonance offset. The source estimates minutes for the powder-averaged calculation.

## Model and method

The Q-band model has one electron and two protons at 80 K. Each calculation uses 17.8 MHz electron nutation frequency, 48 ns pulses, 150 blocks, a 120° second-pulse phase, and a 1600-point spherical powder grid. The script calls `xixdnp` with the TPPM settings, shifts each microwave offset by the −13 MHz reference point, and records the final contact-curve value.

## Scan and output

The 120 offsets span −150 to 150 MHz. Their final proton signals form the plotted field profile.

## Reference

[Redrouthu et al., DOI: 10.1063/5.0153053](https://doi.org/10.1063/5.0153053).
