# examples/nmr_liquids/ct_cosy_2spins.m

- Signature: `ct_cosy_2spins()`

## Purpose

CT COSY spectrum for 2 spins. Calculation time: minutes. Source assignment: [doi:10.1002/jhet.5570250160](http://dx.doi.org/10.1002/jhet.5570250160).

## Physical / mathematical content

This is a constant-time COSY simulation of a two-proton spin system. The two sites have shifts 2.00 and 5.00, with a scalar coupling of 7.0; the spectrum is calculated by Spinach's liquid-state `ct_cosy` sequence.

## Numerical / algorithmic content

The source sets field value 5.9 and uses the full `sphten-liouv` basis (no approximation), offset 500, sweep [2000 2000], 512 points and 2048 zero-fill points on each axis. It applies a squared-cosine apodisation to both dimensions, computes a shifted 2D FFT, and plots the spectrum magnitude in positive mode.

## Implementation structure

The function defines the two proton sites and their interaction, builds the Spinach system and basis, then calls `liquid(...,@ct_cosy,...,'nmr')`. The FID is windowed and Fourier-transformed before plotting.
