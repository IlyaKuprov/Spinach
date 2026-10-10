# kernel/plotting/int_2d.m

- Signature: `int_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs,filename)`

## Purpose

Calls `plot_2d` to draw the contour display, then integrates selected rectangles either from mouse input or from a saved interval file. This is a display-and-integration utility, not a spin-dynamics calculation.

## Contours, axes and interpolation

- `spectrum` is a real 2-D matrix. `parameters.sweep` and `parameters.offset` accept one or two values in Hz; `parameters.axis_units` may be `'ppm'`, `'Hz'` or `'Gauss'`. `parameters.spins` is also required by the delegated plotter: a cell array with one or two isotope-label character strings; one label is reused on both axes. The actual contour rendering and axis construction are delegated to `plot_2d`; this routine uses its returned `f2`, `f1` and `S` coordinates/data.
- The source documents `ncont=20` as a reasonable contour count. `delta=[0.02 0.2 0.02 0.2]` is a suggested positive/negative fractional level range; `k=2` is a reasonable curvature exponent, `ncol=256` is a suitable number of colours, and `m=6` is suggested for higher contrast (with `m=1` giving a linear red/blue ramp). These are usage suggestions, not defaults assigned by `int_2d`.
- With spectrum extrema `smax` and `smin`, levels follow `(delta(2)-delta(1))*smax*linspace(0,1,ncont).^k+smax*delta(1)` and `(delta(4)-delta(3))*smin*linspace(0,1,ncont).^k+smin*delta(3)`. `signs` selects positive, negative or both contour sets.
- For integration it makes `[F1,F2]=ndgrid(f1,f2)` and a spline `griddedInterpolant(F1,F2,transpose(S),'spline')`, then integrates each selected/file-specified rectangle using `integral2` with relative and absolute tolerances `1e-3`. The coordinates and units are those of the plotted axes.

## Bounds, files and side effects

There are no bounds automatically supplied by this wrapper. If `filename` does not exist, `ginput(2)` obtains each rectangle's opposite corners, and each current `ranges` value is saved to that file. If the file exists, `ranges` is loaded and integrated automatically. Integrals and selected ranges are reported through `report`; the function has no return arguments. It turns off MATLAB warning `MATLAB:griddedInterpolant:MeshgridEval2DWarnId` and does not restore the previous warning state.

The wrapper itself checks that `filename` is a character string; plot inputs are checked by the delegated `plot_2d` routine. It creates no extra colormap or axes configuration directly beyond calling `plot_2d`.

## Links

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/int_2d.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=int_2d.m)
