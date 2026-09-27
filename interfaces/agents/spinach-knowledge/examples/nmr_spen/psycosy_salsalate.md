# examples/nmr_spen/psycosy_salsalate.m

- Signature: `psycosy_salsalate()`

## Purpose

PSYCOSY of one salsalate ring. Calculation time: minutes on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

- The source defines the spin system for one salsalate ring by its isotope, chemical-shift, and scalar-coupling data.
- It simulates the `@psycosy` sequence with Spinach's imaging function.

## Numerical / algorithmic content

- The resulting two-dimensional FID is square-sine apodised in each dimension and transformed with a 2D FFT. The script plots the magnitude spectrum.

## Implementation structure

- Sets up the spin system, basis, and sequence parameters; runs imaging and performs the stated apodisation, Fourier transform, and plot.
