# examples/nmr_liquids/ct_hsqc_sucrose.m

- Signature: `ct_hsqc_sucrose()`

## Purpose

CT HSQC spectrum of sucrose with natural content of 13C isotope (magnetic parameters computed with DFT). Calculation time: seconds

## Physical / mathematical content

- Two-dimensional constant-time HSQC of sucrose using 13C and 1H spins. Magnetic parameters are initialized from a vacuum DFT log, then selected isotropic shifts are replaced with experimental values.
- The simulation treats 13C isotopomers separately, applies squared-cosine apodisation to the positive and negative FIDs, forms a States signal, and Fourier transforms both dimensions.

## Numerical / algorithmic content

- Spin-system generation uses `g2spinach` with `min_j=3.0` and `no_xyz=1`; the code then sets the isotropic shifts for listed spins. The basis is sphten-liouv / IK-2 with scalar-coupling connectivity and proximity level 1; greedy settings use `prox_cutoff=4.0`.
- Sequence settings are `J=140`, sweep `[3350 950]`, offset `[5000 1100]`, `npoints=[128 128]`, and `zerofill=[512 512]`; F2 13C is decoupled. Isotopomer calculations use `parfor` (no GPU path is present).

## Implementation structure

- Build the sucrose spin system from the vacuum DFT log and replace the listed isotropic shifts with experimental values; set the field to 5.9 T and define the selected basis and CT-HSQC parameters.
- Generate 13C isotopomers, simulate each in parallel, then apodise the FIDs, form the States signal, Fourier transform both dimensions, and plot the real spectrum.
