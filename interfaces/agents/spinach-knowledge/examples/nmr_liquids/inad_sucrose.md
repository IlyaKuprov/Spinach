# examples/nmr_liquids/inad_sucrose.m

- Signature: `inad_sucrose()`
- Source: [examples/nmr_liquids/inad_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/inad_sucrose.m)

## Purpose

A one-dimensional liquid-state `13C` INADEQUATE example for sucrose. It selects double-quantum coherence from coupled carbon pairs and converts it back for carbon detection; the source estimates calculation time in minutes.

## Implementation

The spin system is imported from `../standard_systems/sucrose.log` using the vacuum-DFT setup, with `1H` and `13C` isotope channels. The code replaces isotropic shifts for spin indices `1:19` and `24:30` with the experimental values `[94.5 73.4 74.9 71.5 74.7 62.4 63.6 106.0 78.7 76.3 83.7 64.7 5.49 3.63 3.83 3.54 3.90 3.90 3.90 3.75 3.75 4.29 4.12 3.96 3.90 3.90]` ppm. The field is `11.7` T. It uses the scalar-coupling IK-1 Liouville basis, with proximity cutoff `4.0` and `options.min_j=5.0`.

The simulated detection channel is `13C`, with `1H` decoupling; the sequence parameters include `J=50` Hz, `offset=10000`, `sweep=8000`, 4096 points, and zero filling to 8192 points. The source does not state units for its offset and sweep fields. It creates pair-labelled `13C` isotopomers, calculates each pair coupling as one third of the trace of its coupling matrix, and simulates only when the magnitude passes the source test `abs(J)>2*pi*1.0`. Each accepted FID receives exponential apodisation with parameter 6; the zero-filled Fourier transforms are accumulated and the real spectrum is plotted in ppm with the axis inverted.
