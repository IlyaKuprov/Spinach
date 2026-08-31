# kernel/pulses/triwave.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/triwave.m`
- Signature: `waveform=triwave(amplitude,frequency,time_grid)`
- Total lines: 49

## Purpose

Returns a triangular waveform. Syntax: waveform=triwave(amplitude,frequency,time_grid)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 24-25: Check consistency; implemented by `grumble(amplitude,frequency,time_grid)`.
- Lines 27-28: Compute the waveform; implemented by `waveform=abs(sawtooth(amplitude,frequency,time_grid))`.

### Key state/data transformations

- Lines 28: computes `waveform` using `waveform=abs(sawtooth(amplitude,frequency,time_grid))`.

### Local helper functions

- Line 33: `grumble()` — `function grumble(amplitude,frequency,time_grid)`.
  - Representative operation: `if (numel(amplitude)~=1)||(~isnumeric(amplitude))||(~isreal(amplitude))`.
  - Representative operation: `error('amplitude parameter must be a real number.')`.

## Parameters / inputs

- amplitude -amplitude at the tooth top
- frequency -waveform frequency, Hz
- time_grid -vector of time points, seconds

## Outputs

- waveform -waveform array of the same shape
- as the time_grid input

## Implementation structure

- Returns a triangular waveform. Syntax:
- waveform=triwave(amplitude,frequency,time_grid)
- amplitude - amplitude at the tooth top
- frequency - waveform frequency, Hz
- time_grid - vector of time points, seconds
- waveform - waveform array of the same shape
- as the time_grid input
- Check consistency
- Compute the waveform
- Consistency enforcement
- I am here to determine whether what you had just done is simple
- incompetence or deliberate sabotage.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `sawtooth()`.
