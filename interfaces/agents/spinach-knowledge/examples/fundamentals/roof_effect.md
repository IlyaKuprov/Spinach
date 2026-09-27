# examples/fundamentals/roof_effect.m

- Signature: `roof_effect()`

## Purpose

Illustrates the roof effect in the NMR spectrum of a strongly J-coupled two-spin system as the two resonance offsets approach one another.

## Physical / mathematical content

- The model has two `1H` spins at 5.9 T, chemical shifts 0.95 and 1.45, and a 7.0 Hz scalar coupling. It uses the `sphten-liouv` basis without approximation.

## Numerical / algorithmic content

- For each separation parameter 0.2, 0.05, 0.0125, and 0.00625, the code updates the two Zeeman frequencies symmetrically about 1.2 ppm, simulates liquid-state acquisition, applies exponential apodisation (10), and Fourier-transforms the FID.
- Acquisition uses 300 Hz offset and sweep, 1024 points, zero filling to 4096, Hz axis units, and an inverted axis.

## Implementation structure

- Sets `L+` as both the initial state and receiver, runs `liquid` with the NMR assumption, applies apodisation and an FFT, and plots the real spectrum for each separation.
