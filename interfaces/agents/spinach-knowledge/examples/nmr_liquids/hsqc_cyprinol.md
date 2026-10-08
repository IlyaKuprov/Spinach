# examples/nmr_liquids/hsqc_cyprinol.m

- MATLAB implementation: [examples/nmr_liquids/hsqc_cyprinol.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hsqc_cyprinol.m)

- Signature: `hsqc_cyprinol()`

## Purpose

Build a natural-abundance 13C HSQC spectrum for cyprinol. The source comment estimates seconds for the calculation.

## Spin system and method

`cyprinol()` supplies 42 `1H` and 27 `13C` sites with scalar-coupling and isotropic-shift data. Its source says values are taken from the cited report where available and otherwise estimated. `dilute(spin_system,'13C')` creates the 13C isotopomer contributions. The calculation uses a sparse Liouville-space basis (`sphten-liouv`, IK-1), scalar-coupling connectivity, proximity level 1 and interaction level 3. Greedy system setup uses proximity cutoff 4.0 and interaction cutoff 5.0; the field is 11.7 T.

## HSQC acquisition

The wrapper passes `{'13C','1H'}` as the F1/F2 spins: 13C is indirect and 1H is the direct detected nucleus. It specifies `J=150 Hz`, sweeps [12000, 2500] Hz, offset values [5000, 1250] (the sources do not state their unit), [128, 128] points and [512, 512] zero-fill points; the axes are in ppm. It requests 1H decoupling in F1 and 13C decoupling in F2.

## Processing and scope

Each isotopomer is simulated with `liquid(...,@hsqc,...,'nmr')`. The positive and negative States FIDs receive square-cosine apodisation in both dimensions. The wrapper Fourier-transforms both dimensions, combines the States data as the source specifies, sums isotope contributions, then plots the real 2D spectrum with positive polarity. It sets no explicit relaxation parameters. The wrapper supplies acquisition settings and data processing; pulse-train details are implemented by `experiments/nmr_liquids/hsqc.m`, not here.

## Sources

Cyprinol shift/coupling source: http://dx.doi.org/10.1002/mrc.4782
HSQC sequence sources: https://doi.org/10.1016/0009-2614(80)80041-8 and https://doi.org/10.1002/cmr.a.10095

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
