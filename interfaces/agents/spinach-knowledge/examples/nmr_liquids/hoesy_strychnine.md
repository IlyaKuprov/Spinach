# examples/nmr_liquids/hoesy_strychnine.m

- MATLAB implementation: [examples/nmr_liquids/hoesy_strychnine.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hoesy_strychnine.m)

- Signature: `hoesy_strychnine()`

## Purpose

Simulates the `13C{1H}` HOESY spectrum of strychnine at natural `13C` content. The source estimates calculation time as minutes; this is source documentation, not a timing measured here.

## Spin system and model

The example obtains a strychnine spin system from `strychnine({'1H','13C'})`, so the explicit spin labels are proton and carbon-13. `dilute(spin_system,'13C')` generates carbon-13 isotopomer subsystems; the script simulates each and accumulates their spectra. The field setting is `14.1` T. The basis is spherical-tensor Liouville space (`sphten-liouv`), `IK-1`, scalar-coupling connectivity, proximity level 3 and interaction level 4. Relaxation is Redfield with IME equilibrium, `rlx_keep='kite'`, correlation time `50e-12` s and temperature `298` K. Algorithm options are `greedy`, proximity cutoff 5.0 and interaction cutoff 2.0.

## Acquisition and processing

Mixing time is `0.5` s. The source orders dimensions as `{'1H','13C'}`, with `decouple_f1={'13C'}`; carbon-13 is the detected nucleus in the source's `13C{1H}` description. Sweeps are `[6000 18000]` Hz, offsets `[3000 12000]` Hz, and each dimension has 128 acquired points and 512-point zero filling. Axes are labelled in ppm.

Each carbon-13 isotopomer is sent to `liquid(subsystem,@hoesy,parameters,'nmr')`. The wrapper applies `sqcos` apodisation to cosine and sine FIDs in both dimensions, zero-fills and Fourier-transforms F2, forms the States signal `f1_cos-1i*f1_sin`, then Fourier-transforms F1 and adds the result to the accumulated spectrum. The real spectrum is plotted with positive display polarity using `plot_2d`.

## Sequence boundary

This wrapper builds the strychnine spin systems, configures acquisition and processing, and calls the shared `liquid` driver with `@hoesy`. It does not define the HOESY pulse-program internals; pulse timing and internal coherence-transfer steps belong to the sequence function.
