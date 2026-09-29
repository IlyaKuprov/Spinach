# kernel/pulses/slr_pulse.m

MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/slr_pulse.m
Source Wiki page: https://spindynamics.org/wiki/index.php?title=slr_pulse.m

- Signature: `[Cx,Cy,durs,amps,phis]=slr_pulse(npts,dur,tbw,flip_angle,pass_rip,stop_rip)`

## Purpose

Designs a Shinnar–Le Roux (SLR) linear-phase selective excitation pulse as `npts` piecewise-constant slices. It returns Cartesian controls, slice durations, and the equivalent RF amplitudes and phases.

## Design

The implementation converts the 90-degree excitation-profile ripple targets to beta targets, estimates a transition width using the Pauly relation, and solves a continuous weighted least-squares problem in a linear-phase cosine basis. The stop-band weight is the ratio of the beta passband and stopband targets. The beta coefficients are scaled by `sin(flip_angle/2)`; the source then checks the sampled beta response and requires a positive complementary spectrum before cepstral factorisation of the minimum-phase alpha polynomial. An inverse SLR recursion recovers one slice rotation at a time.

The ripple inputs are design targets, not guaranteed minimax error bounds. For flip angles below `pi/2`, they do not specify angle-independent magnetisation error bounds.

## Parameters / inputs

- `npts` - even integer of at least 2; number of equal-duration, piecewise-constant slices
- `dur` - finite positive total pulse duration, seconds
- `tbw` - finite positive time-bandwidth product: pulse duration times the nominal full passband width
- `flip_angle` - finite on-resonance flip angle in radians, strictly greater than zero and no greater than `pi/2`
- `pass_rip` - finite dimensionless 90-degree excitation passband ripple target, strictly between zero and one
- `stop_rip` - finite dimensionless 90-degree excitation stopband ripple target, strictly between zero and one

The derived passband and stopband edges must satisfy `0 < pass_edge < stop_edge < 1`; otherwise the function rejects the specification. Each returned slice has duration `dur/npts`.

## Outputs

- `Cx`, `Cy` - X and Y control components, respectively, in rad/s; each is a `1 x npts` row vector
- `durs` - slice durations in seconds, a `1 x npts` row vector
- `amps` - RF amplitudes in rad/s, a `1 x npts` row vector
- `phis` - RF phases in radians, a `1 x npts` row vector

## Reference

J. Pauly, P. Le Roux, D. Nishimura, and A. Macovski, *IEEE Transactions on Medical Imaging* 10(1), 53–65 (1991), https://doi.org/10.1109/42.75611
