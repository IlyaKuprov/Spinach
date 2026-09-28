# examples/nmr_liquids/cosy90_strychnine.m

- Signature: `cosy90_strychnine()`

## Purpose

COSY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

This example obtains a proton spin system from Spinach's `strychnine({'1H'})` helper and simulates its homonuclear liquid-state COSY spectrum. The sequence is run with a 90-degree angle and detects a two-dimensional signal.

## Numerical / algorithmic content

The model sets field value 5.9 and uses the greedy option and proximity cutoff 4.0, with the `sphten-liouv` formalism, IK-2 approximation, scalar-coupling connectivity and proximity level 1. The sequence uses offset 1200, sweep 2200, 512 points and 2048 zero-fill points on each axis. A cosine window is applied on both dimensions before the shifted 2D FFT; the plot uses the real spectrum.

## Implementation structure

The function loads the strychnine proton parameters, sets the field and basis, and passes the resulting system to `liquid(...,@cosy,...,'nmr')`. It then applies the two-dimensional cosine apodisation, performs the zero-filled 2D Fourier transform and plots the spectrum.
