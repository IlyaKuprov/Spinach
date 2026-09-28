# examples/nmr_spen/psyche_rotenone.m

- Signature: `psyche_rotenone()`

## Purpose

PSYCHE pure-shift NMR spectrum of rotenone. Calculation time: hours, faster on a GPU.

## Physical / mathematical content

- The source defines a proton spin system for rotenone at 11.7 T, with chemical shifts and scalar couplings specified explicitly.
- The imaging simulation calls the `@psyche` sequence. It sets a 15 mm sample, a 1H initial/detection state, and saltire-chirp pulse parameters; diffusion is set to zero.

## Numerical / algorithmic content

- The two acquisition dimensions use sweeps of 100 and 5000 Hz, with 32 and 2048 points and zero fills of 128 and 8192 points, respectively.
- The code extracts a pure-shift FID from the imaging output, applies Gaussian apodisation (parameter 6), then computes a 2D Fourier spectrum and a 1D projected spectrum.

## Implementation structure

- Constructs the spin system and basis, configures the sample and pulse sequence, runs imaging, and processes and plots the resulting spectra.
