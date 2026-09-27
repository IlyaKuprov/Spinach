# examples/relaxation_theory/csa_csa_xcorr_2.m

- Signature: `csa_csa_xcorr_2()`

## Purpose

Illustrate CSA–CSA cross-correlation in the `103Rh` subsystem and its effect on the widths of the proton triplet lines. Calculation time: seconds.

## Spin system and relaxation

- Field: `11.75 T`; isotopes: `{'1H','103Rh','103Rh'}`.
- The proton shielding matrix is diagonal with `[6.9 6.9 6.9]`; each rhodium shielding matrix is diagonal with `[7250 8000 7250]`.
- Scalar couplings are 4 Hz between the proton and each rhodium, and 100 Hz between the rhodium nuclei.
- Redfield relaxation uses `tau_c={10e-9}`, zero equilibrium, and secular retention; the basis is `sphten-liouv` with no approximation.

## Acquisition and processing

The liquid-state acquisition uses `1H` with `L+` initial and detection states, no decoupling, offset `6.9*500`, sweep 50, 2048 points, and zero-fill 16384. The axis is in ppm and inverted. The FID is generated with `liquid(...,@acquire,...,'nmr')`, exponentially apodised with 20, Fourier transformed, and plotted.
