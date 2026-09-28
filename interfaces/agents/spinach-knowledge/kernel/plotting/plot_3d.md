# kernel/plotting/plot_3d.m

- Signature: `plot_3d(spin_system,spectrum,parameters,nsurf,delta,k,signs)`

## Purpose

Volume isosurface plotting utility with non-linear adaptive surface spacing. The function plots the 3D volume and the three projections onto the coordinate planes. Syntax: plot_3d(spin_system,spectrum,parameters,nsurf,delta,k,signs)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- `spectrum` — real cube containing the 3D NMR spectrum.
- `parameters.sweep` — three sweep widths, Hz.
- `parameters.spins` — cell array of three character strings specifying the working spins.
- `parameters.offset` — three transmitter offsets, Hz.
- `parameters.axis_units` — `ppm`, `Hz`, or `Gauss`.
- `parameters.npoints` and `parameters.zerofill` — each a three-element vector of finite positive integers; spectrum dimensions must match `parameters.zerofill`.
- `nsurf` — number of surfaces (20 is a reasonable value).
- `delta` — minimum and maximum surface elevations as fractions of total intensity; the first pair applies to positive surfaces and the second to negative ones. A suggested starting value is `[0.02 0.2 0.02 0.2]`.
- `k` — controls surface-spacing curvature: 1 is linear; values above 1 increase sampling near the baseline. A reasonable value is 2.
- `signs` — `positive`, `negative`, or `both`, selecting which surfaces to plot.

## Outputs

- Creates a figure showing the 3D isosurfaces and projections onto the coordinate planes.

## Implementation structure

- Validates the real spectrum cube and plotting parameters, then derives surface levels with `contspacing`.
- Builds and converts frequency axes, draws the isosurfaces, and labels the 3D plot.
- Produces F3–F2, F2–F1, and F3–F1 projections by summing the spectrum along dimensions 1, 3, and 2, respectively, and plots each projection with `plot_2d` using 20 contour levels, independently of `nsurf`.

[Source reference](https://spindynamics.org/wiki/index.php?title=plot_3d.m)
