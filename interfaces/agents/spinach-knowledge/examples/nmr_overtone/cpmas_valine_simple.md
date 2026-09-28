# examples/nmr_overtone/cpmas_valine_simple.m

- Signature: `cpmas_valine_simple()`

## Purpose

Simulates proton-to-`14N` overtone cross-polarisation in N-acetylvaline under MAS with the Fokker–Planck formalism. The source attributes the valine quadrupolar tensor data to [this paper](https://doi.org/10.1039/c4cp03994g) and estimates hours of calculation time.

## Model and calculation

The system is `14N` and `1H` at 14.10220742 T. Nitrogen quadrupole parameters are 3.21 MHz, asymmetry 0.27, and spin 1; nitrogen shifts are [57.5, 81.0, 227.0] and proton shifts are zero. Diagonal damping relaxation is used with rate 2000; the basis is `sphten-liouv` without approximation, and the code disables the Krylov and trajectory-level options.

At the magic angle, the average-treatment calculation uses rank 9, rate −19.840 kHz, grid `rep_2ang_6400pts_sph`, a 70–105 kHz sweep, and 256 points with 256-point zero-fill. It sets the `14N` RF frequency to 86.30 kHz and duration to 100 μs; the RF-power pair is 2π × [55.0, 35.1] kHz divided by sin of the magic angle. The spectrum is calculated with `singlerot` and `overtone_cp`.
