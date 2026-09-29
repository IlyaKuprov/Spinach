# examples/nmr_liquids/hoesy_ftyr_a.m

- MATLAB implementation: [examples/nmr_liquids/hoesy_ftyr_a.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hoesy_ftyr_a.m)

- Signature: `hoesy_ftyr_a()`

## Purpose

Simulates the `1H -> 19F` HOESY experiment for 3-fluorotyrosine. The source selects the transfer direction to minimise the time that aromatic fluorine spends in the transverse plane, noting the short aromatic `19F T2` in proteins.

## Spin system and model

The source reads `../standard_systems/3_fluoro_tyr.log` through `gparse` and `g2spinach`, requesting only `1H` and `19F` spins and passing `[31.82 192.97]` as the absolute isotropic shielding references for `1H` and `19F`, respectively, to place the reference substances at zero ppm. Thus this simulation's explicit spin network is proton plus fluorine-19; it does not request carbon, nitrogen, or oxygen spins. The field setting is `14.1` T. The basis is spherical-tensor Liouville space (`sphten-liouv`), `IK-2`, scalar-coupling connectivity, and proximity level 3. Relaxation is Redfield with IME equilibrium, `rlx_keep='kite'`, correlation time `10e-9` s (commented as a large protein) and temperature `298` K. Algorithm options are `greedy`, proximity cutoff 5.0 and interaction cutoff 2.0.

## Acquisition and processing

The mixing time is `0.5` s (commented as quite long). Dimensions are ordered `{'1H','19F'}`, with `decouple_f1={'19F'}`; fluorine-19 is the detected nucleus in the stated transfer direction. Sweeps are `[4000 2500]` Hz, offsets `[3000 -70000]` Hz, with 128 acquired points and 512-point zero filling in each dimension. The axes are labelled in ppm.

The wrapper calls `liquid(spin_system,@hoesy,parameters,'nmr')`. It applies `sqcos` apodisation to both cosine and sine FIDs in both dimensions, zero-fills and Fourier-transforms F2, forms `f1_cos-1i*f1_sin`, and Fourier-transforms F1. The real spectrum is plotted with negative display polarity using `plot_2d`.

## Sequence boundary

This file supplies the spin-system, acquisition and processing settings to the shared `liquid` driver and passes `@hoesy` as the sequence. The HOESY pulse-program internals are not implemented in this wrapper, so no pulse timings or internal transfer steps are specified here.
