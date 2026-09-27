# examples/nmr_solids/static_powder_trp.m

- Signature: `static_powder_trp()`

## Purpose

Simulates the 13C NMR spectrum of tryptophan powder. The source takes the coordinates and chemical-shift anisotropies from DFT data, substitutes experimental isotropic shifts, and assumes proton decoupling. Estimated runtime: hours.

## Model and basis

Spin-system data are read from `../standard_systems/trp_xray.out`; the field is 14.1 T. The source sets isotropic shifts (ppm) for sites 2–12 as follows: 2: 124.2, 3: 110.1, 4: 118.0, 5: 119.3, 6: 114.7, 7: 107.5, 8: 134.9, 9: 125.0, 10: 26.8, 11: 54.6, and 12: 174.4. The basis is `sphten-liouv` with IK-0 approximation, 15N longitudinal order, +1 projection, and inter-level 3; the interaction and proximity cutoffs are 5.0 and 4.0.

## Simulation and processing

The 13C powder acquisition uses the `rep_2ang_6400pts_sph` grid, 60 kHz sweep, 128 points, 512-point zero-fill, 18000 offset, and an inverted ppm axis. The 13C `L+` state is both the initial and detection state. The powder FID is apodised exponentially with parameter 6, Fourier transformed, and plotted.
