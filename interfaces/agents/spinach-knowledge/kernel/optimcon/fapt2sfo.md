# kernel/optimcon/fapt2sfo.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/optimcon/fapt2sfo.m`
- Signature: `[wave,dt,time_grid]=fapt2sfo(fapt,time_grid)`
- Total lines: 111

## Purpose

Converts a freq-ampl-phase-time specification of a pulse sequ- uence into the corresponding single frequency origin waveform that is compatible with GRAPE optimisations. Syntax: [wave,dt,time_grid]=fapt2sfo(fapt,time_grid)

## Physical / mathematical content

- Optimal-control core routines. These files implement GRAPE-style objective evaluation, quasi-Newton search, line search, regularisation, distortion models, and waveform parameterisations.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(fapt,exist('time_grid','var'))`.
- Lines 38-39: Make time grid if not provided; implemented by `if ~exist('time_grid','var')`.
- Lines 41-42: Find the highest frequency in the sequence; implemented by `max_freq=max(abs(cellfun(@(x)x(1),fapt)))`.
- Lines 44-45: Find the end time of the sequence; implemented by `end_time=max(cellfun(@(x)x(5),fapt))`.
- Lines 47-48: Get the time grid with dt < half Nyquist; implemented by `dt_hN=1/(4*max_freq); npts=ceil(end_time/dt_hN)+1`.
- Lines 53-56: Make sure time grid is a real row vector; implemented by `if (~isnumeric(time_grid))|| (~isreal(time_grid))|| (~isrow(time_grid))`.
- Lines 60-61: Accept what is provided; implemented by `npts=numel(time_grid); dt=[]`.
- Lines 65-66: Preallocate waveform; implemented by `wave=zeros(2,npts)`.
- Lines 68-69: Loop over events; implemented by `for n=1:numel(fapt)`.
- Lines 71-72: Pull the parameters out; implemented by `freq=fapt{n}(1); ampl=fapt{n}(2); phi=fapt{n}(3)`.
- Lines 75-76: Get the time interval mask; implemented by `time_mask=(time_grid>=start_time)&(time_grid<=end_time)`.
- Lines 78-79: Add event to the waveform; implemented by `wave(1,:)=wave(1,:)+ampl*cos(2*pi*freq*time_grid+phi).*time_mask`.

### Control flow inferred from the code

- Line 39: conditional branch on `~exist('time_grid','var')`.
- Line 54: conditional branch on `(~isnumeric(time_grid))||`.
- Line 69: `for` loop over `n=1:numel(fapt)`.

### Key state/data transformations

- Lines 42: computes `max_freq` using `max_freq=max(abs(cellfun(@(x)x(1),fapt)))`.
- Lines 48: computes `dt_hN` using `dt_hN=1/(4*max_freq); npts=ceil(end_time/dt_hN)+1`.
- Lines 49: computes `time_grid` using `time_grid=linspace(0,end_time,npts); dt=time_grid(2)`.
- Lines 61: computes `npts` using `npts=numel(time_grid); dt=[]`.
- Lines 66: computes `wave` using `wave=zeros(2,npts)`.
- Lines 72: computes `freq` using `freq=fapt{n}(1); ampl=fapt{n}(2); phi=fapt{n}(3)`.
- Lines 73: computes `start_time` using `start_time=fapt{n}(4); end_time=fapt{n}(5)`.
- Lines 79: computes `wave(1,:)` using `wave(1,:)=wave(1,:)+ampl*cos(2*pi*freq*time_grid+phi).*time_mask`.
- Lines 80: computes `wave(2,:)` using `wave(2,:)=wave(2,:)-ampl*sin(2*pi*freq*time_grid+phi).*time_mask`.

### Local helper functions

- Line 87: `grumble()` — `function grumble(fapt,have_grid)`.
  - Representative operation: `if ~iscell(fapt)`.
  - Representative operation: `error('fapt must be a cell array of five-element vectors.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `exist()`, `cellfun()`, `time_grid()`, `isrow()`, `wave()`, `iscell()`.
