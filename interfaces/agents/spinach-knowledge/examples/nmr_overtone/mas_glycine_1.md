# examples/nmr_overtone/mas_glycine_1.m

- Signature: `mas_glycine_1()`

## Purpose

Simulates the 14N overtone MAS spectrum of glycine with the Fokker–Planck formalism. The quadrupolar tensor data are attributed to [O'Dell and Ratcliffe](http://dx.doi.org/10.1016/j.cplett.2011.08.030); the simulation parameters reproduce Figure 3b of [the cited paper](http://dx.doi.org/10.1039/C4CP03994G). The source estimates a calculation time of minutes.

## Scientific and numerical content

The model uses a single 14N spin at 14.1 T, with the quadrupolar interaction set by `eeqq2nqi(1.18e6,0.53,1,[0 0 0])` and scalar Zeeman value 32.4. The basis is `sphten-liouv` with no approximation; diagonal damping is specified at rate 300 with zero equilibrium. The script disables Krylov and trajectory-level methods.

At the magic angle, it calls `singlerot` with `@overtone_pa` and the `qnmr` pathway. The setup uses rank 6, rotor rate -19840, the `rep_2ang_6400pts_sph` grid, a 44–52 kHz sweep and 256 points with 256-point zero filling. It prepares and detects using the angle-weighted Lz/Lx combinations, applies a 260 μs pulse (55 kHz power factor and 48 kHz offset), and phases the result by 1.35 rad before plotting its real part.
