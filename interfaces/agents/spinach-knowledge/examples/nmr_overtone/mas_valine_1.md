# examples/nmr_overtone/mas_valine_1.m

- Signature: `mas_valine_1()`

## Purpose

Simulates the 14N overtone MAS spectrum of N-acetylvaline using Fokker–Planck formalism. The valine quadrupolar tensor data are attributed to [the cited paper](http://dx.doi.org/10.1039/c4cp03994g). The source estimates a calculation time of hours.

## Scientific and numerical content

The single 14N spin is set at 14.102 T. Its quadrupolar parameters are `eeqq2nqi(3.21e6,0.27,1,[0 0 0])`; the anisotropic Zeeman eigenvalues are [57.5, 81.0, 227.0] with Euler angles [-90, -90, -17] degrees. The basis is `sphten-liouv` without approximation, with diagonal damping at rate 2000 and zero equilibrium. Krylov and trajectory-level methods are disabled.

At the magic angle, `singlerot` runs `@overtone_pa` using `qnmr`. The setup uses rank 8, rotor rate -19840, grid `rep_2ang_6400pts_sph`, 75–100 kHz sweep, and 256 points plus 256-point zero filling. Angle-weighted Lz/Lx operators define preparation and detection; the 70 μs pulse uses the source's 55 kHz power factor and 86 kHz offset. The plotted real spectrum is phased by 1.75 rad.
