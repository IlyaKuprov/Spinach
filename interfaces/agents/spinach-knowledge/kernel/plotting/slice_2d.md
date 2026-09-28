# kernel/plotting/slice_2d.m

- Signature: `slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`

## Purpose

Plot a 2D NMR spectrum with non-linear adaptive contour spacing and use mouse input to extract and display 1D slices.

## Parameters / inputs

- `spin_system` — spin system used by the plotting routines; working spins must be isotopes present in the system.
- `spectrum` — numeric matrix containing the 2D NMR spectrum, described as real in the source documentation.
- `parameters.sweep` — one or two sweep widths, in Hz.
- `parameters.spins` — cell array of one or two character strings specifying the working spins.
- `parameters.offset` — one or two transmitter offsets, in Hz. If absent, zero offsets are assumed.
- `parameters.zerofill` — one or two positive integer point counts in F1 and F2.
- `parameters.axis_units` — character string specifying axis units (`'ppm'`, `'Hz'`, or `'Gauss'` in the source documentation). If absent, `'ppm'` is assumed.
- `ncont` — positive integer number of contours; 20 is suggested.
- `delta` — four real values between 0 and 1 specifying contour elevations relative to the spectrum intensity. The first pair applies to positive contours and the second to negative contours; `[0.02 0.2 0.02 0.2]` is a suggested starting value.
- `k` — positive integer controlling contour-spacing curvature. `k=1` gives linear spacing; values above 1 increase sampling density near the baseline. A suggested value is 2.
- `ncol` — positive integer number of colours in the colour map; about 256 is suggested.
- `m` — positive integer controlling colour-map curvature. `m=1` gives a linear ramp toward red for positive contours and blue for negative contours; 6 is suggested for high-contrast plotting.
- `signs` — character string specifying `'positive'`, `'negative'`, or `'both'` contours.

## Physical / mathematical content

The source documents these contour-level expressions, where `smin` and `smax` are computed from the spectrum matrix:

- `cont_levs_pos=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);`
- `cont_levs_neg=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);`

## Numerical / algorithmic content

The function applies defaults and checks its inputs, then duplicates single spin, offset, and sweep entries for homonuclear 2D sequences. It calls `plot_2d` for the contour plot. For each mouse-selected point, it constructs F1 and F2 axis grids, uses spline `griddedInterpolant` interpolation to compute two traces, and displays them with `plot_1d` as F1 and F2 slices. Both slice plots use the minimum and maximum values of `spectrum` as their vertical limits. The mouse-selection loop continues indefinitely.

## Outputs

The function creates a figure containing the 2D spectrum and two 1D slice plots; it has no return value.

## Source link

<https://spindynamics.org/wiki/index.php?title=slice_2d.m>