# examples/nmr_proteins/hnca_simple.m

- Signature: `hnca_simple()`

## Purpose

A minimal example of HNCA pulse sequence simulation. Calculation time: seconds.

## Physical / mathematical content

- Simulates a four-spin system containing `15N`, `13C` (CA), `1H`, and `13C` (C) at a magnetic field of 14.1. Scalar Zeeman values are `[110 60 8 180]`; the specified scalar couplings are N–H 92, N–CA 11, N–C 15, CA–C 55, CA–H 2, and H–C 4.
- Uses the `sphten-liouv` basis with no approximation.

## Numerical / algorithmic content

- Calls `liquid(spin_system,@hnca,parameters,'nmr')` with spins `{'15N','13C','1H'}`, sweep widths `[2800 5000 3000]`, offsets `[-7200 8600 5100]`, 64 points per dimension, 256-point zero filling per dimension, and `ppm` axis units.
- Applies squared-cosine apodisation in all three dimensions to each of the four acquired components: `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg`.
- Performs shifted, zero-filled FFTs along F3, F2, then F1. Conjugate combinations of the phase components form the absorption parts of the F3 and F2 signals before the final transform.

## Implementation structure

- Defines the field, spin system, interactions, basis, and sequence parameters; constructs the Spinach system with `create` and `basis`.
- Simulates the HNCA sequence, processes the four signal components into a three-dimensional spectrum, and plots `-real(spectrum)` using `plot_3d` with positive display.