# examples/nmr_liquids/cosy90_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/cosy90_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/cosy90_strychnine.m)

- Signature: `cosy90_strychnine()`

## Purpose

A liquid-state homonuclear proton COSY-90 calculation for strychnine. The source estimates a calculation time of minutes.

## Spin system and basis

Rather than listing shifts and couplings locally, the function imports the proton model with `strychnine({'1H'})`; that helper supplies `sys` and `inter`. The field is set to 5.9 T. The basis is Liouville-space IK-2 with scalar-coupling connectivity and proximity level 1. The greedy option is enabled, with a proximity cutoff of 4.0 (the example does not specify a unit for this cutoff).

## COSY acquisition and processing

The pulse angle is pi/2, the offset is 1200 Hz, and the sweep width is 2200 Hz. The simulation samples 512 by 512 points, then zero-fills to 2048 by 2048 for the 2D FFT; the axes are labelled in ppm. A cosine window is applied in both dimensions. The plotted signal is the real part of the shifted spectrum, displayed with two-dimensional contours.

## Interpretation and scope

The function calls the liquid-state COSY sequence and displays its simulated spectrum for the imported strychnine parameter set. The source does not reproduce the helper's shift or coupling table in this file and does not report an experimental comparison or a numerical peak assignment. The displayed spectrum therefore describes this configured simulation, not an independent measurement.
