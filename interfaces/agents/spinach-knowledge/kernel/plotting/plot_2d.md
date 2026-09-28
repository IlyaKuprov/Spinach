# kernel/plotting/plot_2d.m

- Signature: `[axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,...`

## Purpose

Plots 2D NMR spectra using non-linear adaptive contour spacing to show small cross-peaks alongside large diagonal peaks. Syntax: see the signature above.

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- `spectrum` — real matrix containing the 2D NMR spectrum.
- `parameters.sweep` — one or two sweep widths, Hz.
- `parameters.spins` — cell array with one or two character strings specifying the working spins.
- `parameters.offset` — one or two transmitter offsets, Hz.
- `parameters.axis_units` — `ppm`, `Gauss`, `Hz`, `kHz`, `MHz`, or `points`.
- `ncont` — number of contours (20 is a reasonable value).
- `delta` — minimum and maximum contour elevations as fractions of total intensity; the first pair applies to positive contours and the second to negative contours. A suggested starting value is `[0.02 0.2 0.02 0.2]`.
- `k` — controls contour-spacing curvature: 1 is linear; values above 1 increase sampling near the baseline. A reasonable value is 2.
- `ncol` — number of colors in the color map (around 256).
- `m` — color-map curvature: 1 gives a linear red ramp for positive contours and blue for negative contours; 6 is a reasonable high-contrast value.
- `signs` — `positive`, `negative`, or `both`, selecting which contours to plot.

## Outputs

- Returns `axis_f1` and `axis_f2` ticks and the 2D `spectrum` for external plotting utilities.
- Complex spectra are shown as real and imaginary panels side by side. If a panel is all zero, it has no contours; the function draws empty axes with the correct ranges and labels the panel “all-zero spectrum”.

## Implementation structure

- Applies defaults and validates inputs; for complex data, recursively plots the real and imaginary parts side by side.
- Computes positive and negative contour levels with `contspacing`, constructs axes from sweep widths and offsets, and converts axes to the requested units.
- Transposes the spectrum for plotting, draws its contours (or empty axes for an all-zero panel), applies the color map, and reverses both axes.
- Draws the color bar unless `colorbar` is listed in `spin_system.sys.disable`.

[Source reference](https://spindynamics.org/wiki/index.php?title=plot_2d.m)
