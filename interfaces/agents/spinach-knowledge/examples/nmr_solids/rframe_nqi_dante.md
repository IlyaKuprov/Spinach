# examples/nmr_solids/rframe_nqi_dante.m

- Signature: `rframe_nqi_dante()`
- Source: [examples/nmr_solids/rframe_nqi_dante.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/rframe_nqi_dante.m)

## Purpose

Simulates a DANTE MAS spectrum of one quadrupolar 14N nucleus using a one-dimensional Fokker–Planck equation and a spherical grid. The source says it is set to reproduce Figure 3d of [the cited paper](https://doi.org/10.1016/j.jmr.2012.05.024), and estimates minutes of calculation time.

## Spin system and rotor sampling

The field parameter is 18.8 and the quadrupolar interaction is constructed as `eeqq2nqi(1.18e6, 0.50, 1, [0 0 0])`; the source does not state units for these interaction arguments. The simulation uses the full spherical-tensor Liouville basis, `singlerot` with the DANTE callback in the lab frame, rate 62.5e3, rotor axis `[1,1,1]`, maximum rank 35, and grid `rep_2ang_200pts_sph`. The selected rotating-frame order is 2. No gradient is configured.

## DANTE acquisition and processing

The initial state is the 14N `Lz` state and the receiver is `L+`. The sequence parameters are pulse duration 1.2e-6, pulse amplitude 88e3, two pulses, and two periods. Acquisition uses sweep 2000000, 1024 points, zero-fill to 4096, and offset 2200; the frequency-axis units are explicitly set to Hz. The FID receives exponential apodisation with parameter 6, is Fourier transformed, and the plotted spectrum is its magnitude.
