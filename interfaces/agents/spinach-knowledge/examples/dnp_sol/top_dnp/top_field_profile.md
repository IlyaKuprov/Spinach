# examples/dnp_sol/top_dnp/top_field_profile.m

- Signature: `top_field_profile()`

## Purpose

Calculates the proton (I_z) value at the end of a fixed-contact-time TOP DNP sequence as a function of microwave resonance offset. The source estimates minutes for the powder-averaged calculation.

## Model and method

This Q-band model uses one electron and two protons at 80 K. TOP settings are 17.8 MHz electron nutation frequency, 10 ns pulse duration, 14 ns delay, 300 blocks, and a 1600-point spherical powder grid. For each offset, the script runs `topdnp` with the offset shifted by a −13 MHz reference point and records the final point of the contact curve. The 120 offsets span −150 to 150 MHz; the output is plotted as the field profile.

## Reference

[Redrouthu et al., *Science Advances* (2019), DOI: 10.1126/sciadv.aav6909](https://doi.org/10.1126/sciadv.aav6909).
