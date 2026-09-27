# examples/nmr_solids/mqmas_nqi.m

- Signature: `mqmas_nqi()`

## Purpose

Rotor-synchronous MQMAS spectrum of a 87Rb compound, transmitter set to the isotropic chemical shift. Calculation time: minutes

## Physical / mathematical content

- Models the quadrupolar nucleus `87Rb` using an NQI interaction built by `eeqq2nqi(5e6,0.50,3/2,[0 0 0])`.
- Simulates a rotor-synchronous multiple-quantum MAS experiment with MQ order 3 and the transmitter offset at zero, as stated in the example description.

## Numerical / algorithmic content

- Uses `singlerot` to simulate the two-dimensional lab-frame pulse sequence, then applies squared-cosine apodisation along both dimensions and a two-dimensional Fourier transform with zero filling to `[256 256]`.
- The experiment uses a 62.5 kHz rotor rate, rank-7 orientation grid `rep_2ang_1600pts_sph`, two 128-point dimensions, pulse amplitudes `2π × [250e3 250e3]`, and durations of 2 μs and 1 μs. The plotting call passes `20` as an argument.

## Implementation structure

- Creates a 9.4 T `87Rb` spin system in the `sphten-liouv` basis without approximation and disables trajectory-level output.
- Uses `Lz` as the initial state and `L+` as the receiver; runs `mqmas` through `singlerot` in the lab frame with the specified rotor and pulse parameters.
- Applies the two-dimensional apodisation and Fourier transform, then plots the magnitude spectrum.
