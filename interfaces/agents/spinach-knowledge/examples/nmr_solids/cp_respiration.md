# examples/nmr_solids/cp_respiration.m

- Signature: `cp_respiration()`

## Purpose

1H-13C RESPIRATION-CP experiment in the doubly rotating frame. Magic angle spinning simulation using Fokker-Planck formalism. Calculation time: seconds

## Physical / mathematical content

The model is a two-spin 1H/13C system with a 2.00 Å internuclear separation at 11.7 T. It implements RESPIRATION cross-polarisation under magic-angle spinning, with 1H excitation and 13C detection.

## Numerical / algorithmic content

The code calls `singlerot` with the `respiration` pulse sequence, a 20 kHz spinning rate, maximum rank 8, and the `rep_2ang_100pts_sph` grid. The pulse train has 16 loops and θ=π/20; acquisition uses 512 points over a 40 kHz sweep. The FID is exponentially apodised (parameter 5), zero-filled to 16384 points, Fourier transformed, and plotted in kHz.

## Implementation structure

Creates the spin system and full sphten-liouv basis, sets the Fokker–Planck and acquisition parameters, simulates the signal, then applies exponential apodisation and an FFT before plotting the real spectrum.
