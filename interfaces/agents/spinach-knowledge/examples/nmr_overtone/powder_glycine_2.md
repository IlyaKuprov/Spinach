# examples/nmr_overtone/powder_glycine_2.m

- Signature: `powder_glycine_2()`

## Purpose

Simulates the 14N powder overtone NMR spectrum of glycine using Fokker–Planck formalism. The quadrupolar tensor data are attributed to [O'Dell and Ratcliffe](http://dx.doi.org/10.1016/j.cplett.2011.08.030). The source states that the T2,-2 coherence gives rise to the overtone signal in a static sample and estimates a calculation time of seconds.

## Scientific and numerical content

The model uses 14N at 14.1 T, quadrupolar parameters `eeqq2nqi(1.18e6,0.53,1,[0 0 0])`, scalar Zeeman value 32.4, an unapproximated `sphten-liouv` basis and diagonal damping rate 500 with zero equilibrium. Krylov and trajectory-level methods are disabled.

The powder setup uses grid `rep_2ang_6400pts_sph`, sweep 0–15 kHz and 256 points with 256-point zero filling. The initial state is explicitly `T2,-2` for 14N; the coil is the magic-angle weighted Lz/Lx combination. It calls `powder` with `@overtone_a` and `qnmr`; unlike `powder_glycine_1`, it has no pulse parameters.
