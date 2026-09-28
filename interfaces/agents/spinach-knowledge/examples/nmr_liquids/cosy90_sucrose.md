# examples/nmr_liquids/cosy90_sucrose.m

- Signature: `cosy90_sucrose()`

## Purpose

COSY spectrum of sucrose (magnetic parameters computed with DFT). Calculation time: minutes

## Physical / mathematical content

This is a homonuclear proton COSY simulation for sucrose. Its spin-system parameters are generated from the vacuum DFT log at `../standard_systems/sucrose.log` by `g2spinach`, with hydrogen nuclei mapped to `1H`; the example then models the liquid-state COSY response.

## Numerical / algorithmic content

The DFT conversion uses `min_j=2.0` and `no_xyz=1`, with the conversion argument 31.8. The subsequent Spinach setup sets field value 5.9 and uses greedy mode, proximity cutoff 4.0, an IK-2 scalar-coupling Liouville basis at proximity level 1, and a 90-degree COSY angle. It uses offset 800, sweep 1700, 512 points and 2048 zero-fill points in both dimensions, followed by two-dimensional cosine apodisation and a shifted 2D FFT; the plotted data are the real spectrum.

## Implementation structure

The code converts the DFT log into `sys` and `inter`, builds the basis, runs `liquid(...,@cosy,...,'nmr')`, then windows and Fourier-transforms the FID before plotting.
