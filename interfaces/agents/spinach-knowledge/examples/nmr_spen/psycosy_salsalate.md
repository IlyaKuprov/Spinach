# examples/nmr_spen/psycosy_salsalate.m

- Signature: `psycosy_salsalate()`

## Purpose

PSYCOSY of one salsalate ring. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- The source defines one salsalate-ring spin system with five 1H shifts (8.14, 7.44, 7.71, 7.28, and 7.60 ppm) at 14.1 T, and scalar couplings including 7.9, 7.5, 8.1, 1.6, and 1.2 Hz. The spatial model is 15 mm long with 100 points.
- It simulates `@psycosy` with Spinach's imaging function. The 2D acquisition uses 512 points per dimension, zero-filled to 1024, with 720 Hz sweep and 4620 Hz offset; sequence settings include a 110 ms mixing time, 0.01 T/m gradient, and a 20° saltire chirp (15 ms pulse, 10 kHz sweep, 50 ms gradient duration).

## Numerical / algorithmic content

- The resulting two-dimensional FID is square-sine apodised in each dimension and transformed with a 2D FFT. The script plots the magnitude spectrum.

## Implementation structure

- Sets up the spin system, basis, and sequence parameters; runs imaging and performs the stated apodisation, Fourier transform, and plot.
