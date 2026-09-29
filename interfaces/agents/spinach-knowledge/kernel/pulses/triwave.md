# kernel/pulses/triwave.m

MATLAB source: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/triwave.m
Source Wiki page: https://spindynamics.org/wiki/index.php?title=triwave.m

- Signature: `waveform=triwave(amplitude,frequency,time_grid)`

## Purpose

Returns a triangular waveform evaluated at the supplied time points.

## Waveform and inputs

The function takes the absolute value of the corresponding sawtooth waveform. The sawtooth helper forms `amplitude*(2*frequency*mod(time_grid,1/frequency)-1)`, so the triangle repeats with period `1/frequency`; for positive amplitude it reaches the tooth top at period boundaries and zero halfway through each period. The waveform inherits the shape of `time_grid`.

- `amplitude` - finite real scalar giving the tooth-top amplitude
- `frequency` - positive finite real scalar, in Hz (teeth per second)
- `time_grid` - finite real array of time points, in seconds; the source documentation describes this input as a vector

The helper's finite-value and positive-frequency checks apply when `triwave` calls it.

## Output

- `waveform` - triangular waveform array with the same shape as `time_grid`
