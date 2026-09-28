# examples/dnp_sol/tppm_dnp/tppm_contact_curve.m

- Signature: `tppm_contact_curve()`

## Purpose

Plots the proton (I_z) expectation value during the contact period of a TPPM DNP experiment. The source estimates a runtime of seconds.

## Model and method

The Q-band model has one electron and two protons at 80 K. It uses 33 MHz electron nutation frequency, 16 ns pulses, 250 TPPM blocks, and a second-pulse phase of 120°. The microwave offset is ((-13+4)) MHz relative to the sequence reference. A 400-point spherical powder average calls `xixdnp` with the sequence parameters configured for this TPPM example; the resulting proton signal is plotted against contact time.

## Reference

[Redrouthu et al., DOI: 10.1063/5.0153053](https://doi.org/10.1063/5.0153053).
