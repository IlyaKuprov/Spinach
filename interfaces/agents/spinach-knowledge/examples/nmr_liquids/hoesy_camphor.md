# examples/nmr_liquids/hoesy_camphor.m

- MATLAB implementation: [examples/nmr_liquids/hoesy_camphor.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hoesy_camphor.m)

- Signature: `hoesy_camphor()`

## Purpose

Simulates the `13C{1H}` HOESY spectrum of camphor at natural `13C` content. The source header says the coordinates, shielding anisotropies and J couplings come from vacuum DFT and estimates calculation time as minutes; that is source documentation, not a timing measured here.

## Spin system and model

The example reads `../standard_systems/camphor.log` with `gparse` and `g2spinach`, requesting `1H` and `13C` spins with the source arguments `[31.5 189.2]`; `options.min_j=3.0` and `options.no_xyz=0`. The spin model therefore tracks the specified proton and carbon-13 network, not the other molecular nuclei. `dilute(spin_system,'13C')` generates the carbon-13 isotopomer subsystems; the script simulates each and accumulates their spectra.

The field setting is `14.1` T. The basis is spherical-tensor Liouville space (`sphten-liouv`), `IK-2`, scalar-coupling connectivity, and proximity level 3. Relaxation is Redfield with IME equilibrium, `rlx_keep='kite'`, correlation time `50e-12` s and temperature `298` K. The algorithm options are `greedy`, proximity cutoff 5.0 and interaction cutoff 2.0. Zero track elimination is also explicitly enabled with `zte` in `sys.enable`.

## Acquisition and processing

Mixing time is `0.5` s. The source orders the dimensions as `{'1H','13C'}`, with `decouple_f1={'13C'}`; the detected signal is carbon-13, as stated by the source's `13C{1H}` description. Sweeps are `[1800 9000]` Hz, offsets `[900 4500]` Hz, and both dimensions have 128 acquired points and 512-point zero filling. The plotted axes use ppm.

For each isotopomer, `liquid(subsystem,@hoesy,parameters,'nmr')` returns cosine and sine FIDs. Both are apodised with `sqcos` in both dimensions; the script zero-fills and Fourier-transforms F2, forms the States signal `f1_cos-1i*f1_sin`, then Fourier-transforms F1. The accumulated real spectrum is plotted with positive display polarity using `plot_2d`.

## Sequence boundary

This wrapper configures the spin system and parameters, calls the shared `liquid` driver with `@hoesy`, and processes the returned FIDs. It does not contain the HOESY pulse-program implementation; pulse timing and internal coherence-transfer steps are delegated to that sequence function.
