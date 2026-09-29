# examples/nmr_liquids/hoesy_ftyr_b.m

- MATLAB implementation: [examples/nmr_liquids/hoesy_ftyr_b.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/hoesy_ftyr_b.m)

- Signature: `hoesy_ftyr_b()`

## Purpose

Simulates the reverse-direction `19F -> 1H` HOESY experiment for 3-fluorotyrosine. The source cautions that this is not the preferred direction for proteins because aromatic fluorine has short `19F T2`, while noting that fluorine is phase-encoded in this experiment.

## Spin system and model

The source reads `../standard_systems/3_fluoro_tyr.log` through `gparse` and `g2spinach`, requesting only `1H` and `19F` spins and passing `[31.82 192.97]` as the absolute isotropic shielding references for `1H` and `19F`, respectively, to place the reference substances at zero ppm. The explicit spin network is therefore proton plus fluorine-19, with no carbon, nitrogen, or oxygen spins requested. The field setting is `14.1` T. The basis is spherical-tensor Liouville space (`sphten-liouv`), `IK-2`, scalar-coupling connectivity, and proximity level 3. Relaxation is Redfield with IME equilibrium, `rlx_keep='kite'`, correlation time `10e-9` s (commented as a large protein) and temperature `298` K. Algorithm options are `greedy`, proximity cutoff 5.0 and interaction cutoff 2.0.

## Acquisition and processing

The mixing time is `0.5` s (commented as quite long). Dimensions are ordered `{'19F','1H'}`, with `decouple_f1={'1H'}`; proton is the detected nucleus in this transfer direction. Sweeps are `[2500 4000]` Hz, offsets `[-70000 3000]` Hz, and each dimension has 128 acquired points and 512-point zero filling. Axes are labelled in ppm.

The wrapper calls `liquid(spin_system,@hoesy,parameters,'nmr')`, applies `sqcos` apodisation to cosine and sine FIDs in both dimensions, zero-fills and Fourier-transforms F2, forms `f1_cos-1i*f1_sin`, and Fourier-transforms F1. The real spectrum is plotted with negative display polarity using `plot_2d`.

## Sequence boundary

This file supplies the spin-system, acquisition and processing settings to the shared `liquid` driver and passes `@hoesy` as the sequence. The HOESY pulse-program internals are not implemented in this wrapper, so no pulse timings or internal transfer steps are specified here.
