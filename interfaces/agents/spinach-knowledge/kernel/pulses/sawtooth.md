# kernel/pulses/sawtooth.m

- MATLAB source: [kernel/pulses/sawtooth.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/pulses/sawtooth.m)
- Spinach wiki: [sawtooth.m](https://spindynamics.org/wiki/index.php?title=sawtooth.m)
- Signature: `waveform=sawtooth(amplitude,frequency,time_grid)`

## Purpose

Evaluates a sawtooth waveform directly at the supplied time points.

## Waveform

For each time value, the source computes `amplitude*(2*frequency*mod(time_grid,1/frequency)-1)`. The period is `1/frequency`. For positive amplitude the waveform ramps from `-amplitude` to values just below `amplitude`, then resets to `-amplitude` at each period boundary; the output has the same shape as `time_grid`.

## Inputs and output

- `amplitude` — finite real numeric scalar; the source does not require it to be positive.
- `frequency` — positive finite real numeric scalar, in teeth per second.
- `time_grid` — finite real numeric array of time points, in seconds. Values are evaluated elementwise; a sorted or uniformly spaced grid is not required by the implementation.
- `waveform` — real-valued samples with the same shape as `time_grid`.

The function evaluates the stated formula without an additional window or filtering step.
