# examples/nmr_spen/psyche_strychnine.m

- Signature: `psyche_strychnine()`

## Purpose

PSYCHE pure-shift NMR spectrum of strychnine. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- The source defines the strychnine proton spin system by its isotope, chemical-shift, and scalar-coupling data.
- The imaging simulation calls the `@psyche` sequence and sets sample, 1H state, and saltire-chirp parameters. The diffusion coefficient is set to zero.

## Numerical / algorithmic content

- The sequence produces a two-dimensional imaging FID. The code reconstructs the pure-shift FID from its first acquisition chunk, applies Gaussian apodisation, and Fourier-transforms both the full FID and the pure-shift projection.
- The script plots the magnitude of the 2D spectrum and the imaginary part of the 1D projection.

## Implementation structure

- Builds the spin system and basis, configures sequence parameters, runs imaging with `@psyche`, and performs the FID reconstruction and spectral processing.
