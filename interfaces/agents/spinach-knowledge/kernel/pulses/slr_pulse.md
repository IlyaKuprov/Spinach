# kernel/pulses/slr_pulse.m

- Signature: `[Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)`

## Purpose

Designs a Shinnar-Le Roux (SLR) linear-phase selective excitation pulse and returns its X/Y controls, slice durations, RF amplitudes, and phases.

## Design

The beta polynomial is obtained by continuous weighted least squares in a linear-phase cosine basis. The complementary minimum-phase alpha polynomial and RF waveform are obtained by the inverse SLR transform. The ripple arguments enter the excitation-pulse transform and the transition-width estimate of Pauly et al.; they are design targets, not guaranteed minimax error bounds. For flip angles below `pi/2`, they do not specify angle-independent magnetisation error bounds. The output controls are calibrated for Spinach propagation under `exp(-1i*H*t)` and may be passed directly to `shaped_pulse_xy()`.

## Parameters / inputs

- `npts` - even number of piecewise-constant pulse slices
- `dur` - total pulse duration, seconds
- `tbw` - time-bandwidth product, defined as pulse duration times the nominal full passband width
- `flip_angle` - on-resonance flip angle between zero and `pi/2`, radians
- `pass_rip` - dimensionless 90-degree excitation passband ripple target used in the prototype design
- `stop_rip` - dimensionless 90-degree excitation stopband ripple target used in the prototype design

## Outputs

- `Cx` - X control amplitudes, rad/s, 1 x npts row vector
- `Cy` - Y control amplitudes, rad/s, 1 x npts row vector
- `durs` - pulse slice durations, seconds, 1 x npts row vector
- `amps` - RF amplitudes, rad/s, 1 x npts row vector
- `phis` - RF phases, radians, 1 x npts row vector

## Reference

J. Pauly, P. Le Roux, D. Nishimura, and A. Macovski, *IEEE Transactions on Medical Imaging* 10(1), 53-65 (1991), https://doi.org/10.1109/42.75611

Source Wiki page: https://spindynamics.org/wiki/index.php?title=slr_pulse.m
