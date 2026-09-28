# examples/nmr_solids/fslghetcor_mas_gly.m

- Signature: `fslghetcor_mas_gly()`

## Purpose

FSLG-HETCOR of alpha-glycine powder under MAS. Calculation time: hours on NVidia Tesla A100, much longer on CPU

## Physical / mathematical content

The script uses a 400 MHz spectrometer and constructs an alpha-glycine spin system from a PCM-DFT `glycine.log` input and assigns the isotropic shifts shown in the source for CO, CA, CA protons, and NH protons. It simulates FSLG heteronuclear correlation under MAS, starting from 1H longitudinal magnetisation and detecting in quadrature on 13C.

## Numerical / algorithmic content

The basis is IK-0 with interaction level 4; interactions below an angular-frequency cutoff of 2π×200 rad/s (equivalent to 200 Hz) are ignored. The calculation uses a rank-7 expansion and `rep_2ang_100pts_sph`, with a 10 kHz MAS rate and four FSLG blocks. The source specifies 83 kHz high-power irradiation, 60/50 kHz CP powers, 100 μs CP duration, and a 128×512 acquisition zero-filled to 512×2048. It applies square-cosine apodisation in both dimensions, forms the quadrature FSLG-HETCOR signal, and Fourier-transforms both dimensions; the F1 sweep is calculated from the high-power field and block count.

## Implementation structure

Parses the glycine spin system and chemical shifts, sets up the reduced basis and two-dimensional experiment, simulates with `singlerot` and `fslghetcor`, processes the cosine and sine channels, and plots the real 2D spectrum. GPU enablement is present only as a commented option in the source.
