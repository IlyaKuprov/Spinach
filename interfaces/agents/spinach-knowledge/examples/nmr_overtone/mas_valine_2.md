# examples/nmr_overtone/mas_valine_2.m

- Signature: `mas_valine_2()`

## Purpose

Simulates 14N Z-detected overtone MAS NMR of N-acetylvaline in the Fokker–Planck formalism. The quadrupolar tensor data are attributed to [the cited paper](http://dx.doi.org/10.1039/c4cp03994g). The source estimates a calculation time of hours.

## Scientific and numerical content

The spin and interactions match the companion valine example: 14N at 14.102 T, quadrupolar parameters `eeqq2nqi(3.21e6,0.27,1,[0 0 0])`, and Zeeman eigenvalues [57.5, 81.0, 227.0] with Euler angles [-90, -90, -17] degrees. It uses an unapproximated `sphten-liouv` basis, diagonal damping rate 2000, zero equilibrium, and disables Krylov and trajectory-level methods.

The rank-8 setup uses rotor rate -19840, grid `rep_2ang_6400pts_sph`, sweep 75–100 kHz, and 256 points with 256-point zero filling. The initial state is 14N Lz and the coil detects the same Lz operator. The script calls `singlerot` with `@overtone_a` and `qnmr`; it does not set the pulse and phase parameters used by `mas_valine_1`.
