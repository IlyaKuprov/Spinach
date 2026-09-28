# examples/nmr_liquids/ct_cosy_three_spin.m

- Signature: `ct_cosy_three_spin()`

## Purpose

CT-COSY of three spin system. Calculation time: seconds

## Physical / mathematical content

This constant-time COSY example defines a three-proton system with chemical shifts 2.70, 4.10 and 6.50, and pairwise scalar couplings 10, 8 and 4. The liquid-state `ct_cosy` simulation produces its two-dimensional signal.

## Numerical / algorithmic content

The calculation uses the full `sphten-liouv` basis, greedy mode and proximity cutoff 4.0. It sets field value 14.1, offset 2700, sweep [3500 3500], 256 points and 512 zero-fill points on each axis. Squared-cosine apodisation precedes the shifted 2D FFT; plotting uses spectrum magnitude in positive mode.

## Implementation structure

The function defines the three-site spin system and couplings, builds the Spinach basis, runs `liquid(...,@ct_cosy,...,'nmr')`, then apodises, Fourier-transforms and plots the signal.
