# examples/dnp_sol/top_dnp/top_contact_curve.m

- Signature: `top_contact_curve()`

## Purpose

Plots the proton (I_z) expectation value during the contact period of a time-optimised pulsed TOP DNP experiment. The source estimates a runtime of seconds.

## Model and method

The Q-band model (1.2142 T, 80 K) contains one electron and two protons at the specified three-spin coordinates. It uses the Zeeman Hilbert basis and detects proton (I_z). The TOP sequence is evaluated by `topdnp` with 300 blocks, 10 ns pulses, 14 ns delays, 17.8 MHz electron nutation frequency, and a 3200-point spherical powder grid.

The electron offset is set to ((-13+92.5)) MHz relative to the sequence reference. The returned contact curve is plotted against elapsed contact time.

## Reference

[Redrouthu et al., *Science Advances* (2019), DOI: 10.1126/sciadv.aav6909](https://doi.org/10.1126/sciadv.aav6909).
