# kernel/pulses/sawtooth.m

- Signature: `waveform=sawtooth(amplitude,frequency,time_grid)`

## Purpose

Returns a sawtooth waveform evaluated at the supplied time points.

## Algorithm

For each time `t`, the function evaluates `amplitude*(2*frequency*mod(t,1/frequency)-1)`. The waveform rises linearly from `-amplitude` to just below `amplitude` during each period `1/frequency`, then resets; the result has the same shape as `time_grid`.

## Parameters / inputs

- `amplitude` — finite real scalar setting the magnitude at the tooth top.
- `frequency` — positive finite real frequency, in teeth per second.
- `time_grid` — finite real array of time points, in seconds.

## Output

- `waveform` — sawtooth values with the same shape as `time_grid`.

[Spinach wiki page](https://spindynamics.org/wiki/index.php?title=sawtooth.m)
