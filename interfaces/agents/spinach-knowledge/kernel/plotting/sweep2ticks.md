# kernel/plotting/sweep2ticks.m

- Source: [kernel/plotting/sweep2ticks.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/sweep2ticks.m)
- Wiki: [sweep2ticks.m](https://spindynamics.org/wiki/index.php?title=sweep2ticks.m)

## Purpose

Convert an offset, sweep width, and point count into a column vector of frequency ticks in Hz, for example as an axis supplied to a plotting function.

## Inputs and discretisation

- `offs` is a real scalar offset from the carrier frequency in Hz.
- `sweep` is a real scalar sweep width in Hz.
- `npoints` is a real integer scalar of at least 1.

The implementation is `-linspace(-sweep/2,sweep/2,npoints)' + offs`. For a positive sweep and more than one point, this yields evenly spaced samples from `offs+sweep/2` down to `offs-sweep/2`, including both endpoints. It returns the axis only: it does not plot, change axes, or write a file. The input checks require real scalar offset and sweep values and an integer point count; they do not impose a positivity check on `sweep`.
