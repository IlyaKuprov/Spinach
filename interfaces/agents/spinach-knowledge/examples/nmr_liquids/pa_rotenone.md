# examples/nmr_liquids/pa_rotenone.m

- Signature: `pa_rotenone()`
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pa_rotenone.m)

## Purpose

Simulates a liquid-state pulse-acquire 1H NMR FID for rotenone. The source uses a T1/T2 relaxation model but does not implement an inversion-recovery sequence: acquisition is through `liquid(...,@acquire,...,'nmr')`. It is not an INADEQUATE, NOE, or NOESY example.

## System and relaxation

The source defines 22 1H spins at `sys.magnet=5.9` and lists chemical shifts and scalar couplings. Relaxation is `inter.relaxation={'t1_t2'}`, retaining diagonal terms; the code assigns `r1_rates` as 1.0 for each spin and `r2_rates` as 3.0 for each spin, with `equilibrium='zero'`. These are the source's rate values; this file does not state their units. The basis is `sphten-liouv` / `IK-2`, scalar-coupling connected, proximal level 1, with S3 symmetry on three specified three-spin groups.

## Acquisition and processing

The observed spins, initial state, and receiver are all 1H, using `L+` for state and receiver and no decoupling. Settings are `offset=1200`, `sweep=2000`, `npoints=4096`, and `zerofill=16536`; the axis is ppm and inverted. Offset and sweep units are not specified in this source. The FID receives Gaussian apodisation parameter 10, followed by a shifted zero-filled Fourier transform and plotting of its real part.

## Source limits

The magnetic parameters are cited to [the reported rotenone study](https://doi.org/10.1002/jhet.5570250160), and the source estimates seconds of computation. This is a simulation recipe, not a supplied measured spectrum or a report of experimental relaxation-rate measurements.
