# examples/dnp_sol/tppm_dnp/tppm_parameter_scan.m

- Signature: `tppm_parameter_scan()`

## Purpose

Maps the final proton (I_z) value after a fixed-contact-time TPPM DNP sequence against electron nutation frequency and microwave resonance offset. The source estimates minutes for the powder-averaged scan.

## Model and method

The Q-band model has one electron and two protons at 80 K. Each sequence uses 48 ns pulses, 150 blocks, a 120° second-pulse phase, and a 400-point spherical powder grid. For each pair of scan values, the script calls `xixdnp` with the TPPM parameters and records the final contact-curve value.

## Scan and output

Electron nutation frequency spans 10–50 MHz at 30 points; the 120 offsets span −100 to 100 MHz, with a −13 MHz reference shift in the sequence offset. The script draws an updating contour plot.

## Reference

[Redrouthu et al., DOI: 10.1063/5.0153053](https://doi.org/10.1063/5.0153053).
