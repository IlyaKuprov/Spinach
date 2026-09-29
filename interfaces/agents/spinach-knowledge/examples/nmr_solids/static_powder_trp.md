# examples/nmr_solids/static_powder_trp.m

- Signature: `static_powder_trp()`
- Source: [examples/nmr_solids/static_powder_trp.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_solids/static_powder_trp.m)

## Purpose and model

Calculates a static powder 13C NMR spectrum of tryptophan; the source estimates hours of runtime. It builds the spin system from `../standard_systems/trp_xray.out` with `gparse` and `g2spinach`, selecting 13C and 15N. The source describes coordinates and chemical-shift anisotropies as DFT-derived, and replaces isotropic shifts at sites 2-12 with experimental values: 2: 124.2, 3: 110.1, 4: 118.0, 5: 119.3, 6: 114.7, 7: 107.5, 8: 134.9, 9: 125.0, 10: 26.8, 11: 54.6, and 12: 174.4. The field parameter is 14.1. The `sphten-liouv` basis uses IK-0, 15N longitudinal order, projection +1, and inter-level 3; interaction and proximity cutoffs are 5.0 and 4.0, and trajectory-level algorithms are disabled.

## Powder acquisition and processing

The static powder average uses `rep_2ang_6400pts_sph`. Acquisition selects 13C with sweep 6e4, 128 points, 512-point zero-fill, and offset 18000; the frequency axis is labelled ppm and inverted. Initial and detection states are both 13C `L+`. The source assumes proton decoupling, but sets `decouple` to empty and specifies no explicit pulse sequence, rotor, or gradient. The FID is exponentially apodised with parameter 6, Fourier transformed, and plotted as its real spectrum.
