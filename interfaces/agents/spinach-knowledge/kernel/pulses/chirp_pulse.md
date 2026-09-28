# kernel/pulses/chirp_pulse.m

- Signature: `[Cx,Cy,durs,ints,amps,phis,frqs]=chirp_pulse(npts,dur,bwidth,smp,type)`

## Purpose

Generates a frequency-swept chirp pulse with a WURST or smoothed amplitude envelope, or a saltire pulse formed from a smoothed chirp. The waveform is calibrated to produce an inversion pulse.

## Inputs

- `npts` — positive integer number of waveform points.
- `dur` — finite positive pulse duration in seconds.
- `bwidth` — finite positive sweep bandwidth around zero frequency, in Hz.
- `type` — one of `wurst`, `smoothed`, or `saltire`; append `-adaptive` for adaptive sampling.
- `smp` — envelope parameter. For WURST, it is the power in `1-abs(sin(pi*time_grid).^smp)` and must exceed 1. For smoothed and saltire pulses it sets the edge-fade duration as a percentage from 0 to 50; 0 gives a square envelope and 50 a sine-bell envelope.

## Outputs

- `Cx`, `Cy` — real and imaginary Cartesian waveform components in rad/s.
- `durs` — time-slice durations for piecewise-constant propagation, in seconds.
- `ints` — interval durations for piecewise-linear propagation, in seconds.
- `amps`, `phis`, `frqs` — waveform amplitude (rad/s), phase (rad), and instantaneous frequency (Hz). The saltire branch sets `Cy` to zero and returns before assigning `frqs`; requesting more than six outputs for a saltire pulse is rejected.

The adaptive mode uses a nonuniform time grid; otherwise the time grid is uniform. The phase is quadratic in the normalized time coordinate and the instantaneous frequency sweeps linearly across the specified bandwidth. The amplitude envelope is scaled by `2*pi*sqrt(bwidth/dur)`.

## Reference

[Spin Dynamics Wiki: `chirp_pulse.m`](https://spindynamics.org/wiki/index.php?title=chirp_pulse.m)
