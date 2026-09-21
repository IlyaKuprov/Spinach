# kernel/pulses/sawtooth.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/pulses/sawtooth.m`
- Signature: `waveform=sawtooth(amplitude,frequency,time_grid)`
- Total lines: 52

## Purpose

Returns a saw-tooth waveform. Syntax: waveform=sawtooth(amplitude,frequency,time_grid)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- amplitude -amplitude at the tooth top
- frequency -waveform frequency in teeth per second
- time_grid -grid of time points, seconds

## Outputs

- waveform -waveform array of the same shape
- as the time_grid input

## Implementation structure

- Returns a saw-tooth waveform. Syntax:
- waveform=sawtooth(amplitude,frequency,time_grid)
- amplitude - amplitude at the tooth top
- frequency - waveform frequency in teeth per second
- time_grid - grid of time points, seconds
- waveform - waveform array of the same shape
- as the time_grid input
- Check consistency
- Compute the waveform
- Consistency enforcement
- They're all aristocrats, that's true... because they know that there's
- no such thing as a lousy job -only lousy men who don't care to do it.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `any()`.
