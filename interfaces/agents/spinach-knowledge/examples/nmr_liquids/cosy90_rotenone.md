# examples/nmr_liquids/cosy90_rotenone.m

- MATLAB implementation: [examples/nmr_liquids/cosy90_rotenone.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/cosy90_rotenone.m)

- Signature: `cosy90_rotenone()`

## Purpose

A liquid-state, homonuclear proton COSY-90 simulation for a 22-site rotenone spin model. The source cites [doi:10.1002/jhet.5570250160](http://dx.doi.org/10.1002/jhet.5570250160) and estimates a calculation time of minutes.

## Spin system and basis

All 22 sites are `1H`. The source assigns chemical shifts (ppm) of 6.72, 6.40, 4.13, 4.56, 4.89, 6.46, 7.79, 3.79, 2.91, 3.27, 5.19, 4.89, 5.03, 1.72, 1.72, 1.72, 3.72, 3.72, 3.72, 3.76, 3.76, and 3.76. Pairwise scalar couplings are specified in the source in Hz, with listed nonzero values from 0.7 to 15.8 Hz; for example, J(3,4)=12.1 Hz, J(9,10)=15.8 Hz, J(10,11)=9.8 Hz, and J(9,11)=8.1 Hz. The field is 5.9 T.

The Liouville-space basis uses the IK-2 approximation, scalar-coupling connectivity, and proximity level 1. The source enables the greedy system-building option and groups sites 14-16, 17-19, and 20-22 as three S3 symmetry groups. Zero track elimination is also explicitly enabled with `zte` in `sys.enable`.

## COSY acquisition and processing

The sequence uses a pi/2 pulse, offset 1200 Hz, and sweep width 2000 Hz. The sampled 2D FID is 512 by 512 points and is zero-filled to 2048 by 2048 for the 2D FFT; the displayed axes use ppm. Cosine apodisation is applied along both dimensions. The plotted array is the real part of the shifted spectrum, with a two-dimensional contour display.

## Interpretation and scope

The output is the simulated real-component COSY spectrum of the specified spin model, not a reported measurement. The example records its model shifts, coupling network, basis approximation, and processing choices; it does not provide a separate experimental data set or a quantitative comparison to one.
