# examples/nmr_solids/static_powder_gly.m

- Signature: `static_powder_gly()`
- Source: [examples/nmr_solids/static_powder_gly.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/static_powder_gly.m)

## Purpose

A static-powder 13C NMR calculation for glycine using magnetic parameters read from a DFT log. The source assumes proton decoupling and records a calculation time of seconds.

## Spin model and powder average

The parsed glycine model includes 13C and 15N spins; protons are not included, so their decoupling is an assumption rather than an explicitly simulated channel. The field parameter is 14.1 (no unit is stated in the source). The basis retains 15N longitudinal terms with projection +1, and the source sets interaction and proximity cutoffs to 5.0 and 4.0. `powder` performs the static orientation average on `rep_2ang_6400pts_sph`; no rotor or gradient sequence is specified.

## Acquisition and processing

The selected 13C channel starts from and detects `L+`; no spins are listed for decoupling. Acquisition uses sweep 5e4, offset 18000, 128 points, and zero filling to 512. The displayed axis is configured in ppm; units for sweep and offset are not stated. The FID receives exponential apodisation with parameter 6, then a zero-filled Fourier transform; the real spectrum is plotted.
