# examples/nmr_liquids/hsqc_strychnine.m

- Signature: `hsqc_strychnine()`

## Purpose

HSQC spectrum of strychnine with natural content of 13C isotope. Calculation time: minutes

## Physical / mathematical content

The example simulates a strychnine HSQC spectrum with natural 13C abundance by diluting the spin system into 13C isotopomers. The sequence specifies J=140 Hz, observes 13C and 1H, and decouples 1H in F1 and 13C in F2.

## Numerical / algorithmic content

The IK-2 sphten-liouv basis uses scalar-coupling connectivity and proximity level 1. Simulations run over isotopomers in a `parfor` loop; square-cosine apodisation and two Fourier transforms with States combination form the 2D spectrum.

## Implementation structure

The field is 5.9 T. Acquisition uses sweeps [10000 3000] Hz, offsets [4000 1000] Hz, 128 points per dimension and 512-point zero filling. The greedy algorithm uses proximity cutoff 4.0. The real spectrum is plotted with positive polarity.
