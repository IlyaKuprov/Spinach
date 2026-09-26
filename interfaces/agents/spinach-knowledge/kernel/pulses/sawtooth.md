# kernel/pulses/sawtooth.m

- Signature: `waveform=sawtooth(amplitude,frequency,time_grid)`

## Purpose

Returns a saw-tooth waveform. Syntax: waveform=sawtooth(amplitude,frequency,time_grid)

## Physical / mathematical content

- Pulse and waveform utilities. These files encode shaped RF pulses, gradient events, rotating-frame transformations, resonator response, and Lie-group integration of time-dependent driven dynamics.

## Numerical / algorithmic content

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
