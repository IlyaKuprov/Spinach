# examples/dnp_sol/top_dnp/top_parameter_scan.m

- Signature: `top_parameter_scan()`

## Purpose

Maps the final proton (I_z) value after a fixed-contact-time TOP DNP sequence against electron nutation frequency and microwave resonance offset. The source estimates minutes for the powder-averaged scan.

## Model and method

The Q-band model contains one electron and two protons at 80 K. Each TOP simulation uses 10 ns pulses, 14 ns delays, 300 blocks, and a 400-point spherical powder grid. The script calls `topdnp` for each point and records the final contact-curve value.

## Scan and output

The electron nutation frequency spans 10–50 MHz at 30 points; the 120 offsets span −100 to 100 MHz, with an additional −13 MHz reference shift in the sequence offset. A contour plot is updated across the offset loop.

## Reference

[Redrouthu et al., *Science Advances* (2019), DOI: 10.1126/sciadv.aav6909](https://doi.org/10.1126/sciadv.aav6909).
