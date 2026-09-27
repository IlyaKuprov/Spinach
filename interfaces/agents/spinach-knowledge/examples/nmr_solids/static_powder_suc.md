# examples/nmr_solids/static_powder_suc.m

- Signature: `static_powder_suc()`

## Purpose

Simulates the 13C NMR spectrum of static sucrose powder, assuming proton decoupling. The source estimates a calculation time of hours.

## Model and basis

The spin system and interactions are read from `../standard_systems/sucrose.log` with `gparse` and `g2spinach` (PCM DFT data), at 14.1 T. The calculation uses the `sphten-liouv` formalism, IK-0 approximation, +1 projection, and inter-level 3. It disables trajectory-level algorithms and sets the interaction and proximity cutoffs to 5.0 and 4.0, respectively.

## Simulation and processing

The static 13C powder acquisition uses the `rep_2ang_800pts_sph` grid, 50 kHz sweep, 128 points, 512-point zero-fill, and 15000 offset; the axis is in ppm and inverted. The 13C `L+` state is used for both initial state and detection. The FID is apodised exponentially with parameter 6, Fourier transformed, and plotted.
