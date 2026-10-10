# examples/nmr_liquids/hsqc_sucrose.m

- MATLAB implementation: [examples/nmr_liquids/hsqc_sucrose.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hsqc_sucrose.m)

- Signature: `hsqc_sucrose()`

## Purpose

Build a natural-abundance 13C HSQC spectrum for sucrose from the supplied vacuum-DFT output. The source comment estimates seconds for the calculation.

## Spin system and method

The wrapper parses `examples/standard_systems/sucrose.log` and calls `g2spinach` with H mapped to `1H` and C to `13C`; the input contains 22 proton and 12 carbon sites, while oxygen is not imported as a spin. The DFT-derived couplings are filtered with `options.min_j=3.0 Hz`, and `options.no_xyz=1`. The script replaces isotropic shifts at spin indices [1:19, 24:30] with its 26-entry experimental-shift array, then generates the 13C isotopomer contributions. The DFT reference arguments are [31.8, 182.1]; the wrapper does not annotate their units. The basis is sparse Liouville space (`sphten-liouv`, IK-2), scalar-coupling connectivity and proximity level 1. Greedy setup uses proximity cutoff 4.0; the field is 5.9 T.

## HSQC acquisition

The F1/F2 spins are `{'13C','1H'}`, making 13C indirect and 1H directly detected. The wrapper specifies `J=140 Hz`, sweeps [3350, 950] Hz, offset values [5000, 1100] (the sources do not state their unit), [128, 128] points and [512, 512] zero-fill points; axes are in ppm. It requests 1H decoupling in F1 and 13C decoupling in F2.

## Processing and scope

Each isotopomer is simulated with `liquid(...,@hsqc,...,'nmr')`. The positive and negative States FIDs receive cosine apodisation in both dimensions; the script Fourier-transforms both dimensions, combines the States channels as the source specifies, and sums contributions. It plots the real 2D spectrum with positive polarity. No explicit relaxation parameters are set in this wrapper. It selects the molecular parameters and acquisition/processing settings; pulse-train details are in `experiments/nmr_liquids/hsqc.m`. The source does not identify the DFT functional or basis set in this wrapper.

## Sources

HSQC sequence sources: https://doi.org/10.1016/0009-2614(80)80041-8 and https://doi.org/10.1002/cmr.a.10095

Zero track elimination is explicitly enabled with `zte` in `sys.enable`.
