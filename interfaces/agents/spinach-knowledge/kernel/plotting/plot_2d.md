# kernel/plotting/plot_2d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/plot_2d.m`
- Signature: `[axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,...`
- Total lines: 262

## Purpose

Contour plotting utility with non-linear adaptive contour spacing. The function is useful for NMR data where small cross-peaks must be adequa- tely contoured next to large diagonal peaks. Syntax: [axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,... parameters,ncont,delta,... k,ncol,m,signs)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `defaults()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep -one or two sweep widths, Hz
- parameters.spins -cell array with one ot two character
- strings specifying the working spins
- parameters.offset -one or two transmitter offsets, Hz
- parameters.axis_units -axis units ('ppm','Hz','Gauss')
- ncont -the number of contours, a reasonable value is 20
- delta -minimum and maximum elevation (as a fraction of the
- total intensity) of the contours above the baseline.
- A good starting value is [0.02 0.2 0.02 0.2]. The
- first pair of numbers refers to the positive conto-
- urs and the second pair to the negative ones.
- k -a coefficient that controls the curvature of the contour
- spacing function: k=1 corresponds to linear spacing and
- k>1 bends the spacing curve to increase the sampling den-
- sity near the baseline. A reasonable value is 2.
- ncol -number of colours in the colour map; around 256 is fine
- m -the curvature of the colour map: m=1 corresponds to a li-
- near colour ramp into the red for positive contours, and
- into the blue for negative contours. A reasonable value
- for high-contrast plotting is 6.
- signs -can be set to 'positive', 'negative' or 'both' -this
- will cause the corresponding contours to be plotted.

## Outputs

- a figure is drawn and the following parameters returned
- axis_f1, axis_f2 -F1 and F2 axis ticks for external
- plotting utilities
- spectrum -2D spectrum array for external
- plotting utilities
- Note: the following functions are used to compute contour levels:
- cont_levs_pos=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
- cont_levs_neg=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);
- where smin and smax are computed from the spectrum matrix.

## Implementation structure

- Contour plotting utility with non-linear adaptive contour spacing. The
- function is useful for NMR data where small cross-peaks must be adequa-
- tely contoured next to large diagonal peaks. Syntax:
- [axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,...
- parameters,ncont,delta,...
- k,ncol,m,signs)
- spectrum -a real matrix containing the 2D NMR spectrum
- parameters.sweep - one or two sweep widths, Hz
- parameters.spins - cell array with one ot two character
- strings specifying the working spins
- parameters.offset - one or two transmitter offsets, Hz
- parameters.axis_units - axis units ('ppm','Hz','Gauss')

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `defaults()`, `grumble()`, `report()`, `nnz()`, `subplot()`, `ktitle()`, `transpose()`, `contspacing()`, `isscalar()`, `ft_axis()`, `spin()`, `set()`, `kxlabel()`, `kylabel()`, `any()`, `colormap()`.
