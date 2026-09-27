# examples/nmr_overtone/powder_glycine_1.m

- Signature: `powder_glycine_1()`

## Purpose

Simulates the 14N powder overtone NMR spectrum of glycine using Fokker–Planck formalism. The quadrupolar tensor data are attributed to [O'Dell and Ratcliffe](http://dx.doi.org/10.1016/j.cplett.2011.08.030). The source notes that it uses a very short pulse with unphysically large power and estimates a calculation time of seconds.

## Scientific and numerical content

The model uses 14N at 14.1 T, quadrupolar parameters `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`, scalar Zeeman value 32.4, an unapproximated `sphten-liouv` basis and diagonal damping rate 500 with zero equilibrium. Krylov and trajectory-level methods are disabled.

The powder calculation uses the `rep_2ang_6400pts_sph` grid, 0–15 kHz sweep and 256 points with 256-point zero filling. At the magic angle, the initial 14N Lz state and angle-weighted Lz/Lx coil and Lx operator are used. It calls `powder` with `@overtone_pa` and `qnmr`; the pulse is 1 μs with the source's `2*pi*11.3e6/sin(theta)` power expression and 10 kHz offset.
