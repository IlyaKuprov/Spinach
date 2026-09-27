# examples/nmr_liquids/mqs_propanol.m

- Signature: `mqs_propanol()`

## Purpose

Multiple-quantum NMR experiment for a propanol spin system. Calculation time: seconds.

## Implementation

The model has seven 1H spins at 14.1 T, with the listed 2J and 3J couplings and S2, S2, and S3 permutation symmetry. The sequence selects third-quantum coherence and acquires a two-dimensional spectrum at three delays (0.0333, 0.0710, and 0.5000 s). Each spectrum is sine-apodised, Fourier transformed in both dimensions, and plotted separately.
