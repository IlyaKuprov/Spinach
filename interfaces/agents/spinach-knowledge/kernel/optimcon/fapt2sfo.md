# kernel/optimcon/fapt2sfo.m

- Signature: `[wave,dt,time_grid]=fapt2sfo(fapt,time_grid)`

## Purpose

Converts frequency-amplitude-phase-time pulse events into a two-row X/Y waveform for GRAPE. Each event contributes between its start and end times, inclusive; overlapping events add.

## Parameters / inputs

- `fapt` — cell array of real five-element vectors: [frequency (Hz), amplitude (rad/s), phase at t=0 (radians), start time (seconds), end time (seconds)]. Amplitudes must be nonnegative, and each end time must exceed its start time.
- `time_grid` — optional real row vector of time ticks. If supplied, it is used as given and `dt=[]`.

## Outputs

- `wave` — 2-by-N waveform; row 1 is X and row 2 is Y. For an event with amplitude A, frequency f and phase phi, its contribution is `X=A*cos(2*pi*f*t+phi)` and `Y=A*sin(2*pi*f*t+phi)` on the event's time interval.
- `dt` — time step for an automatically generated grid; empty when `time_grid` is supplied.
- `time_grid` — row vector of time ticks used to construct the waveform.

## Sampling and rotation convention

Without an explicit grid, the routine samples from zero to the latest event end time with spacing no greater than `1/(4*max(abs(frequency)))`; an all-zero frequency list therefore requires an explicit grid. The positive sign in the Y row gives anticlockwise rotation: in a drift offset by `2*pi*f*Lz`, an event at frequency f is on resonance. Reversing the Y sign gives the opposite rotation sense and is off resonance by 2*f for nonzero f; at f=0 it reflects the nutation axis to -phi.

[Source documentation](https://spindynamics.org/wiki/index.php?title=fapt2sfo.m)