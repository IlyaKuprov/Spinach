# examples/nmr_proteins/noesyhsqc_ubiquitin_deut.m

- Signature: `noesyhsqc_ubiquitin_deut()`

## Purpose

A 1H–1H–15N NOESY-HSQC spectrum of 15N-labelled ubiquitin at 900 MHz with a 90 ms mixing time. The protein is not 13C-labelled; selected positions are deuterated, and deuterium nuclei are represented explicitly as spin-1 particles. The source estimates a week on 32 cores and 512 GB of RAM.

## Setup and acquisition

The example imports all atoms from `1D3Z.pdb` / `1D3Z.bmrb` and deuterates the listed H* positions; the field is 21.1356 T. It uses inter/proximal cutoffs 2.0/4.0. Relaxation includes Redfield and T1/T2 terms, assigns R1 and R2 rates of 100 to deuterium spins, keeps the kite portion, sets zero equilibrium, and uses `tau_c=1e-8` s. The IK-1 sphten-liouv basis uses scalar-coupling connectivity and levels 4/3; propagation caching and greedy optimization are enabled, and asymptotic Redfield is disabled. The code removes `13C` spins. Sequence parameters are `tmix=0.090` s, `J=90.0`, points [128, 64, 128], zero-fill sizes [512, 256, 512], spins `1H` / `15N` / `1H`, offsets [4250, -10600, 4250], sweeps [10750, 3000, 10750], and ppm axes.

## Processing

The `liquid` simulation uses `@noesyhsqc`. Four phase-cycle FIDs receive squared-cosine apodisation; conjugate components are combined to form absorption-mode signals through the F3 and F2 transforms, followed by F1 transformation and plotting of the negative real 3D spectrum.
