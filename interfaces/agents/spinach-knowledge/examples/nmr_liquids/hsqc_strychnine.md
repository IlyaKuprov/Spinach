# examples/nmr_liquids/hsqc_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/hsqc_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hsqc_strychnine.m)

- Signature: `hsqc_strychnine()`

## Purpose

Build a natural-abundance 13C HSQC spectrum for strychnine. The source comment estimates minutes for the calculation.

## Spin system and method

`strychnine({'1H','13C'})` selects 22 `1H` and 21 `13C` sites from the molecule definition; its two `15N` sites are excluded by that isotope selection. The script then generates 13C isotopomer contributions with `dilute(spin_system,'13C')`. The basis is sparse Liouville space (`sphten-liouv`, IK-2), scalar-coupling connectivity and proximity level 1. Greedy setup uses proximity cutoff 4.0; the field is 5.9 T.

## HSQC acquisition

The wrapper orders the F1/F2 spins as `{'13C','1H'}`, so 13C is indirect and 1H is the direct detected nucleus. It specifies `J=140 Hz`, sweeps [10000, 3000] Hz, offset values [4000, 1000] (the sources do not state their unit), [128, 128] points and [512, 512] zero-fill points; axes are in ppm. It requests 1H decoupling in F1 and 13C decoupling in F2.

## Processing and scope

Each isotopomer is simulated with `liquid(...,@hsqc,...,'nmr')`. Both States channels receive square-cosine apodisation in both dimensions; the wrapper Fourier-transforms both dimensions, combines the States channels as the source specifies, and accumulates the result. It plots the real 2D spectrum with positive polarity. No explicit relaxation parameters are set by this wrapper. Pulse-train details belong to `experiments/nmr_liquids/hsqc.m`; this script supplies the system, basis, acquisition settings and processing.

## Sources

Strychnine parameters are attributed in the molecule source to Berger and Braun's `200 and more NMR experiments: a practical course`; the source also cites the one-bond C18-H18b coupling and major-conformer coordinates:
http://dx.doi.org/10.1016/j.jmr.2014.02.003
http://dx.doi.org/10.1039/C0CC04114A
HSQC sequence sources: https://doi.org/10.1016/0009-2614(80)80041-8 and https://doi.org/10.1002/cmr.a.10095

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
