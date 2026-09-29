# kernel/plotting/slice_2d.m

- Signature: `slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)`

## Purpose

Display a 2D spectrum and interactively plot two 1D traces through a selected point. This is a rendering and sampling helper, not a propagation or spectrum-calculation routine.

## Inputs and contour settings

`spectrum` is a numeric 2D matrix. The routine applies defaults for a missing `parameters.offset` (zeros) and `parameters.axis_units` (`'ppm'`); a single spin, offset, or sweep entry is duplicated for the two dimensions. It passes the spectrum and parameters to `plot_2d`, which supplies the plotted `f1`, `f2`, and contour matrix `S`; this function does not independently derive those frequency axes.

The source suggests `ncont=20`. Its documented `delta` example is `[0.02 0.2 0.02 0.2]`; the first pair describes positive contour elevation bounds and the second pair negative bounds, as fractions of the total intensity. `k=1` gives linear contour spacing, while `k>1` adds sampling near the baseline; `k=2` is a suggested value. `ncol` controls the number of colormap entries (around 256 is suggested). The source describes `m` as colormap curvature: at `m=1` the ramp is linear into red for positive contours and blue for negative contours; `m=6` is suggested for high contrast. These values are passed to `plot_2d`; `slice_2d` itself does not define the colormap.

## Slice construction and display effects

The figure is divided into three side-by-side subplots and its width is set to three times the root default figure width. In the first panel, `ginput(1)` records a selected point `(x,y)`. The code builds `[F1,F2]=ndgrid(f1,f2)`, constructs a spline `griddedInterpolant(F1,F2,transpose(S))`, then samples one trace at fixed `f1=x` across `f2` and another at fixed `f2=y` across `f1`. The traces are passed to `plot_1d` in the other two panels; both 1D axes receive `YLim=[min(spectrum(:)) max(spectrum(:))]`. The loop is `while true`, so selection continues indefinitely unless interrupted outside this routine. It also disables MATLAB warning `MATLAB:griddedInterpolant:MeshgridEval2DWarnId` and does not restore that warning state in the function.

## Checks and references

The local guards require a numeric matrix; one or two spin entries; one or two offsets, sweep widths, and zero-fill values; a character-string axis_units value (the default is ppm; the plotting helper handles unit conversion); positive integer values for `ncont`, `k`, `ncol`, and `m`; four real `delta` values from 0 to 1; and a character-string `signs` value. Missing offsets and units are filled before these checks; other required plotting inputs are not defaulted. The routine creates/updates a figure and returns no value.

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/slice_2d.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=slice_2d.m)
