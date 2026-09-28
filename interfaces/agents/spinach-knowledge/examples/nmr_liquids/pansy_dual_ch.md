# examples/nmr_liquids/pansy_dual_ch.m

- Signature: `pansy_dual_ch()`

## Purpose

Simulates PANSY-COSY spectra of camphor with natural ¹³C abundance. Coordinates, shieldings, and J-couplings are computed with DFT; the source estimates a calculation time of seconds.

## Physical / mathematical content

- Uses DFT-derived molecular coordinates and magnetic parameters, and generates isotopomers to represent the natural ¹³C content.
- Simulates both PANSY-COSY signal channels.

## Numerical / algorithmic content

- Loops over the generated isotopomers in parallel, applies apodisation and two-dimensional Fourier transforms to both signal components, and plots both spectra.

## Implementation structure

- Reads the DFT data, generates isotopomers, preallocates results, then builds and simulates each spin-system basis in a parallel loop.
