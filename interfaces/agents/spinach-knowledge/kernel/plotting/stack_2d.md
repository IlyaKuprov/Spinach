# kernel/plotting/stack_2d.m

- Source: [kernel/plotting/stack_2d.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/stack_2d.m)
- Wiki: [stack_2d.m](https://spindynamics.org/wiki/index.php?title=stack_2d.m)

## Purpose

Draw a 2D NMR spectrum as stacked line traces in the current axes. Spectrum values are used as trace heights; the function builds frequency or point-index coordinates and does not rescale the spectrum amplitudes.

## Inputs

- `spin_system` supplies the spin and field information used for axis conversion and reporting.
- `spectrum` is a numeric 2D array. Its columns provide the F1 samples and its rows the F2 samples; F2 is horizontal and F1 vertical in the plot.
- `parameters.sweep` contains one or two sweep widths in Hz. A single width is duplicated for both dimensions.
- `parameters.offset` contains one or two transmitter offsets in Hz and defaults to zero when omitted.
- `parameters.spins` is a cell array of one or two spin labels. `parameters.axis_units` defaults to `ppm`; supported values are `ppm`, `Gauss`, `Hz`, `kHz`, `MHz`, and `points`.
- `stack_dim` selects 1 or 2. With 1, each trace runs along F1 and traces are taken at successive F2 positions; with 2, each trace runs along F2 and traces are taken at successive F1 positions.
- `alpha_fun` optionally maps each slice to its line opacity. The default is `@(x)sqrt(norm(x,2))`.

## Plot construction and effects

Frequency coordinates are made with `ft_axis` using offsets, sweep widths, and the matching spectrum dimension. `ppm` and `Gauss` axes use the spin and field data; `Hz`, `kHz`, and `MHz` rescale the frequency axis, while `points` uses sample indices. If any imaginary data are nonzero, the function recursively plots the real part and then the imaginary part on the same axes.
Each slice opacity is evaluated by `alpha_fun`, normalised by the largest opacity, and floored at `0.01`. The function draws patch-line traces, tightens the horizontal limits, enables the grid and perspective projection, orbits the camera, reverses both horizontal axis directions, and labels F1/F2. It updates the current figure and returns no data or file.
