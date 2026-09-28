# kernel/pulses/triwave.m

- Signature: `waveform=triwave(amplitude,frequency,time_grid)`

## Purpose

Returns a triangular waveform.

## Numerical / algorithmic content

The function computes the absolute value of the sawtooth waveform: `abs(sawtooth(amplitude,frequency,time_grid))`.

## Parameters / inputs

- `amplitude` - amplitude at the tooth top
- `frequency` - waveform frequency, Hz
- `time_grid` - vector of time points, seconds

## Outputs

- `waveform` - waveform array with the same shape as `time_grid`

Source Wiki page: https://spindynamics.org/wiki/index.php?title=triwave.m
