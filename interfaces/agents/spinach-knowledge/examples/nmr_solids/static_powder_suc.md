# examples/nmr_solids/static_powder_suc.m

- Signature: `static_powder_suc()`
- Source: [examples/nmr_solids/static_powder_suc.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/static_powder_suc.m)

## Purpose and model

Calculates a static powder 13C NMR spectrum of sucrose; the source assumes proton decoupling and estimates hours of runtime. The spin system and interactions are generated from `../standard_systems/sucrose.log` with `gparse` and `g2spinach` (the source labels the input as PCM DFT data), selecting C/13C, at field parameter 14.1. The basis is `sphten-liouv` with IK-0 approximation, projection +1, and inter-level 3. The source sets interaction and proximity cutoffs to 5.0 and 4.0 and disables trajectory-level algorithms.

## Powder acquisition and processing

`powder` uses the `rep_2ang_800pts_sph` orientation grid. Acquisition parameters are sweep 5e4, 128 points, 512-point zero-fill, and offset 15000; the frequency axis is labelled ppm and inverted. Both initial and detection states are 13C `L+`. Although the source describes a proton-decoupled spectrum, its `decouple` setting is empty and it specifies no explicit pulse sequence, rotor, or gradient. The FID is exponentially apodised with parameter 6, Fourier transformed, and its real spectrum is plotted.
