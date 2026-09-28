# kernel/pulses/restrans.m

- Signature: `[X,Y,dt]=restrans(X_user,Y_user,dt_user,omega,Q,model,up_factor)`

## Purpose

Models the RLC circuit response of a probe, converting an ideal in-phase and out-of-phase pulse waveform into the waveform after the circuit response.

## Algorithm

The input is interpolated onto a finer time grid: `pwc` treats the supplied samples as slice midpoints and uses nearest-neighbour interpolation; `pwl` and `pwl_tsc` treat them as slice edges and use linear interpolation. The function forms the carrier-modulated input, applies a second-order RLC transfer function, demodulates and low-pass filters the output, then downsamples according to `up_factor`. For `pwl_tsc`, it also shifts the output time grid by `2*Q/omega`. With no output arguments, it produces diagnostic plots.

## Parameters / inputs

- `X_user`, `Y_user` — real column vectors for the in-phase and out-of-phase rotating-frame waveform components; they must have the same length.
- `dt_user` — finite positive input slice duration, in seconds; it must not be less than `pi/omega`.
- `omega` — finite positive scalar RLC resonance frequency, in radians per second.
- `Q` — finite positive scalar RLC quality factor.
- `model` — `'pwc'` (piecewise-constant), `'pwl'` (piecewise-linear), or `'pwl_tsc'` (piecewise-linear with time-shift compensation).
- `up_factor` — finite positive integer controlling output waveform discretisation relative to the input; the source describes about 100 as a safe guess.

## Outputs

- `X`, `Y` — in-phase and out-of-phase rotating-frame components after the RLC response.
- `dt` — output slice duration, in seconds.

[Spinach wiki page](https://spindynamics.org/wiki/index.php?title=restrans.m)
