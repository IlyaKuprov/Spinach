# examples/nmr_liquids/pa_styrene.m

- Signature: `pa_styrene()`

## Purpose

Simulates a ¹H NMR spectrum of styrene in Hilbert space to demonstrate parallel propagation, as described in [the cited paper](https://doi.org/10.1063/1.3679656). The source estimates a calculation time of seconds.

## Physical / mathematical content

- Represents the coupled styrene proton spin system in Zeeman-Hilbert space.

## Numerical / algorithmic content

- Uses parallel propagation, acquires the FID, applies Gaussian apodisation, Fourier-transforms it, and plots the spectrum.

## Implementation structure

- Defines the styrene spin system and acquisition, propagates in Hilbert space, then processes and plots the signal.
