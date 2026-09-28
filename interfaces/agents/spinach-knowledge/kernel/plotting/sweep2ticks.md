# kernel/plotting/sweep2ticks.m

- Signature: `axis_hz=sweep2ticks(offs,sweep,npoints)`

## Purpose

Converts an offset, sweep width, and point count into a frequency axis in Hz, suitable for use with MATLAB plotting functions such as `plot()`.

## Physical / mathematical content

The axis is centered on `offs` and runs from `offs+sweep/2` to `offs-sweep/2`, with `npoints` evenly spaced ticks.

## Numerical / algorithmic content

The column vector is constructed as `axis_hz=-linspace(-sweep/2,sweep/2,npoints)'+offs`.

## Parameters / inputs

- `offs` — offset from carrier frequency, Hz; must be a real numeric scalar.
- `sweep` — sweep width, Hz; must be a real numeric scalar.
- `npoints` — number of points in the spectrum; must be a real numeric scalar integer of at least 1.

## Outputs

- `axis_hz` — column vector of frequency-axis ticks, Hz.

## Implementation structure

- Calls `grumble(offs,sweep,npoints)` to check the inputs.
- Builds the axis using `linspace`, transposes it into a column vector, reverses its direction, and adds the offset.

## Reference

- [Spinach documentation: sweep2ticks.m](https://spindynamics.org/wiki/index.php?title=sweep2ticks.m)