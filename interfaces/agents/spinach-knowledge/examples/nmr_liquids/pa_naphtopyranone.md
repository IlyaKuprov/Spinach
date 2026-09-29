# examples/nmr_liquids/pa_naphtopyranone.m

- Signature: `pa_naphtopyranone()`
- Source: [MATLAB example](https://github.com/IlyaKuprov/Spinach/blob/main/examples/nmr_liquids/pa_naphtopyranone.m)

## Purpose

Simulates a liquid-state pulse-acquire 1H NMR spectrum for 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one. It is a proton acquisition through `liquid(...,@acquire,...,'nmr')`, not an INADEQUATE, inversion-recovery, NOE, or NOESY sequence.

## System and interactions

The source specifies twelve 1H spins at `sys.magnet=14.095` T, their scalar shifts, and pairwise scalar couplings. The coded shift list is `[8.345, 7.741, 8.097, 8.354, 7.784, 8.330, 7.059, 7.941, 7.466, 7.326, 7.466, 7.941]`; the displayed chemical-shift axis is ppm. Representative coded couplings are `J(1,2)=7.8`, `J(1,3)=0.9`, and `J(4,5)=8.4`. It uses a scalar-coupling-connected `sphten-liouv` basis with `IK-2` approximation and proximal level 1.

## Acquisition and processing

The source uses 1H `L+` for both initial state and receiver, with no decoupling. It sets `offset=4600`, `sweep=1200`, `npoints=4096`, and `zerofill=32768`, with ppm units and inverted axis. The source does not label offset or sweep units. After liquid-state acquisition it applies exponential apodisation parameter 10, computes the shifted zero-filled Fourier transform, and plots the real spectrum.

## Source limits

Magnetic parameters are cited to [the original report](https://doi.org/10.1016/j.saa.2010.11.015); the example comments estimate seconds of computation. No relaxation model or measured/computed peak values are supplied by this source file.
