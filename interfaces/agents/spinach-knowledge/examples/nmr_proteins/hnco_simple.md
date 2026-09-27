# examples/nmr_proteins/hnco_simple.m

- Signature: `hnco_simple()`

## Purpose

A minimal example of HNCO pulse sequence simulation. Calculation time: seconds.

## Physical / mathematical content

- Simulates a four-spin system at a 14.1 T magnetic field: `15N` (N), `13C` (CA), `1H` (H), and `13C` (C).
- Scalar Zeeman values are `[110 55 8.0 180]`. Nonzero scalar couplings are N–H `92`, N–CA `11`, N–C `15`, and CA–C `55`.
- Uses the `sphten-liouv` formalism with no basis approximation.

## Numerical / algorithmic content

- Runs `hnco` through `liquid(spin_system,@hnco,parameters,'nmr')` with spins `{'15N','13C','1H'}`, offsets `[-7200 25000 4800]`, sweeps `[5000 10000 5000]`, points `[63 64 65]`, zero-fill sizes `[255 256 257]`, delays `[2.25e-3 14e-3 4e-3]`, `f1_decouple=1`, and axis units of `ppm`.
- Applies squared-cosine apodisation in all three dimensions to each of the four returned coherence components: `pos_pos`, `pos_neg`, `neg_pos`, and `neg_neg`.
- Performs shifted, zero-filled Fourier transforms successively along dimensions 3, 2, and 1. Between transforms, combines components with their complex conjugates to form the absorption parts of the F3 and F2 signals.

## Implementation structure

1. Define the magnetic field, spin system, scalar interactions, and basis; initialize the Spinach spin system with `create` and `basis`.
2. Set sequence parameters and simulate the HNCO signal with `liquid`.
3. Apodise all four signal components, transform along F3, and combine conjugate components into `f3_pos` and `f3_neg`.
4. Transform along F2, combine into `f3f2`, then transform along F1 to obtain `spectrum`.
5. Create a figure and plot `real(spectrum)` using `plot_3d` with threshold `10`, view bounds `[0.2 0.9 0.2 0.9]`, dimension `2`, and `'positive'` contours.