# kernel/optimcon/fapt2sfo.m

- Signature: `[wave,dt,time_grid]=fapt2sfo(fapt,time_grid)`

## Purpose

Converts a freq-ampl-phase-time specification of a pulse sequ- uence into the corresponding single frequency origin waveform that is compatible with GRAPE optimisations. Syntax: [wave,dt,time_grid]=fapt2sfo(fapt,time_grid)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

## Parameters / inputs

- fapt -a cell array of 5-element row vectors with of
- the following structure: [frequency (Hz), amp-
- litude (rad/s), phase at t=0 (radians), start
- time (seconds), end time (seconds)]
- time_grid -optional vector of time grid ticks; when
- not provided, the grid is made at twice
- the Nyquist-Shannon minimum sampling ra-
- te of the highest frequency present
- Output:
- wave -pulse sequence as a single waveform, a matrix
- with two rows, corresponding to X and Y compo-
- nents
- dt -step duration of the time grid, seconds
- time_grid -row vector of time grid ticks

## Implementation structure

- Converts a freq-ampl-phase-time specification of a pulse sequ-
- uence into the corresponding single frequency origin waveform
- that is compatible with GRAPE optimisations. Syntax:
- [wave,dt,time_grid]=fapt2sfo(fapt,time_grid)
- fapt -a cell array of 5-element row vectors with of
- the following structure: [frequency (Hz), amp-
- litude (rad/s), phase at t=0 (radians), start
- time (seconds), end time (seconds)]
- time_grid -optional vector of time grid ticks; when
- not provided, the grid is made at twice
- the Nyquist-Shannon minimum sampling ra-
- te of the highest frequency present
