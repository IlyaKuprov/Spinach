# kernel/plotting/stack_2d.m

- Signature: `stack_2d(spin_system,spectrum,parameters,stack_dim,alpha_fun)`

## Purpose

Plot a 2D NMR spectrum as a stack of lines in the current figure.

## Parameters / inputs

- `spin_system` — spin-system data used for axis conversion and reporting.
- `spectrum` — numeric 2D spectrum. If it has a nonzero imaginary component, the real and imaginary components are plotted separately.
- `parameters.sweep` — one or two sweep widths in Hz.
- `parameters.spins` — cell array of one or two character strings specifying the working spins.
- `parameters.offset` — one or two transmitter offsets in Hz; defaults to zero offsets if omitted.
- `parameters.axis_units` — `ppm`, `Gauss`, `Hz`, `kHz`, `MHz`, or `points`; defaults to `ppm` if omitted.
- `stack_dim` — stacking dimension, `1` or `2`.
- `alpha_fun` — optional function handle applied to each spectral slice to determine stack-line opacity. Defaults to `@(x)sqrt(norm(x,2))`.

## Outputs

The function updates the current figure; it does not return an output argument.

## Numerical / algorithmic content

The function constructs F1 and F2 axes from the offsets, sweep widths, and spectrum dimensions, then converts and labels them in the requested units. A single spin, offset, or sweep value is duplicated for both dimensions. For `stack_dim=1` or `stack_dim=2`, it draws lines from slices in the corresponding direction. Slice opacity values are divided by their maximum and values below `0.01` are set to `0.01`. The plot uses a perspective projection, reverses both horizontal axes, and labels F1 and F2.

## Reference

- [Spinach `stack_2d.m` documentation](https://spindynamics.org/wiki/index.php?title=stack_2d.m)