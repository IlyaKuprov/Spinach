# kernel/pulses/chirp_pulse.m

[Source: `kernel/pulses/chirp_pulse.m`](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/chirp_pulse.m)

- Signature: `[Cx,Cy,durs,ints,amps,phis,frqs]=chirp_pulse(npts,dur,bwidth,smp,type)`

## Purpose

Builds a frequency-swept RF waveform calibrated for an inversion pulse. The sweep is centred on zero and linear in time; its phase is quadratic in normalised time. The supported families are WURST, smoothed, and saltire chirps, with an optional `-adaptive` suffix for a nonlinear sample grid.

## Inputs and discretisation

- `npts` — finite positive integer number of waveform points.
- `dur` — finite positive pulse duration in seconds.
- `bwidth` — finite positive sweep bandwidth in Hz, centred on zero.
- `type` — `'wurst'`, `'smoothed'`, or `'saltire'`; append `'-adaptive'` to use the nonlinear normalised time grid instead of the uniform grid.
- `smp` — WURST edge power for the WURST family (the source accepts values of at least 1); for smoothed and saltire it is the percentage of duration affected by the quarter-sine edge ramps, from 0 (square envelope) through 50 (sine-bell envelope).

With a uniform grid the function returns N point samples, N piecewise-constant slice durations in `durs` summing to `dur`, and N−1 piecewise-linear interval durations in `ints`. The adaptive option uses a nonlinear normalised grid and returns the corresponding nonuniform durations. The phase is `pi*dur*bwidth*t.^2` and frequency is `bwidth*t` on normalised time `t`; the amplitude envelope is calibrated by `2*pi*sqrt(bwidth/dur)`.

## Outputs and checks

- `Cx`, `Cy` — real and imaginary RF components in rad/s. For saltire, `Cy` is identically zero; the real waveform's sign is represented by phases of 0 or pi.
- `durs`, `ints` — piecewise-constant slice durations and piecewise-linear interval durations, respectively, in seconds.
- `amps`, `phis`, `frqs` — waveform amplitude in rad/s, phase in radians, and instantaneous frequency in Hz.

The routine checks the scalar/range constraints on its parameters and rejects inadequate phase sampling when a sample-to-sample phase jump exceeds pi and fewer than seven outputs are requested. Requesting the seventh output, `frqs`, bypasses that particular error check. The function returns waveforms and grids only; it does not write files or configure hardware.

## Reference

[Spin Dynamics Wiki: `chirp_pulse.m`](https://spindynamics.org/wiki/index.php?title=chirp_pulse.m)
