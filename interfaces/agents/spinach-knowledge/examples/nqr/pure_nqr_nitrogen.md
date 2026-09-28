# examples/nqr/pure_nqr_nitrogen.m

- Signature: `pure_nqr_nitrogen()`

## Purpose

Powder NQR spectrum of a system with a single (^{14}mathrm{N}) nucleus at zero magnetic field. Calculation time: seconds.

## Physical / mathematical content

- The model contains one (^{14}mathrm{N}) spin and a quadrupolar interaction specified by `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`. With `sys.magnet=0`, the spectrum is generated without a Zeeman field; the transition frequencies are governed by the quadrupolar interaction and its asymmetry parameter.

## Numerical / algorithmic content

- A powder calculation uses `rep_2ang_200pts_sph` and `hp_acquire` with a 5 MHz sweep and 512 points. The acquisition uses (L_+) detection and an (L_x) pulse operator with a (pi/2) pulse angle.
- The code applies laboratory-frame damping at (10^5), sets zero equilibrium and temperature 298, applies exponential apodisation with parameter 6, then plots the imaginary part of the shifted Fourier transform. The frequency axis is labelled in MHz.

## Implementation structure

- Set zero field, isotope, and quadrupolar coupling; select the sphten-liouv basis with no approximation; configure damping and zero equilibrium; construct the Spinach system and basis; set powder acquisition parameters; calculate the powder FID, apodise, Fourier transform, and plot.
