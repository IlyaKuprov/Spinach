# kernel/plotting/plot_3d.m

- Signature: `plot_3d(spin_system,spectrum,parameters,nsurf,delta,k,signs)`

## Purpose

Render an input real 3D spectrum as isosurfaces, with its three coordinate-plane projections. This is a plotting routine; it does not propagate a spin system or calculate a spectrum.

## Inputs and axes

- `spectrum` is a real numeric cube with dimensions `N1 × N2 × N3`. The three sweep widths and offsets are in Hz; `parameters.spins` has three entries.
- Each frequency vector is built by `ft_axis(offset(i),sweep(i),size(spectrum,i))`. That helper starts with `N+1` evenly spaced samples from `-sweep/2` to `sweep/2`; for odd `N` it drops the first and shifts the remaining values by half a bin, while for even `N` it drops the last, then adds the offset. [See `ft_axis`](ft_axis.md).
- `axis_units='Hz'` leaves those vectors unchanged. For `ppm`, each is converted as `1e6*(2*pi*frequency)/(spin(spins{i})*spin_system.inter.magnet)`. For frequency-swept Gauss, each becomes `1e4*(B0-2*pi*frequency/spin('E'))`, where `B0=spin_system.inter.magnet`. These conversions change tick coordinates, not the spectrum values.

## Surface levels and projections

The routine obtains `xmax` and `xmin` from the cube and passes them with `delta`, `k`, `signs`, and `nsurf` to `contspacing` to make the isosurface levels. The source documents `nsurf=20` as a reasonable value. `delta` contains four fractions: the first pair sets positive-surface minimum/maximum elevations and the second pair sets the negative ones; its example is `[0.02 0.2 0.02 0.2]`. `k=1` gives linear spacing; `k>1` concentrates levels nearer the baseline, with `k=2` suggested. `signs` selects `'positive'`, `'negative'`, or `'both'`.

The 3D panel uses `meshgrid(axis_f2,axis_f1,axis_f3)`, so its X, Y, and Z coordinates are F2, F1, and F3 and align with the three spectrum dimensions. Isosurfaces are red with no edges; explicit axis limits are the minimum and maximum of each converted axis vector. The other panels call `plot_2d` on `squeeze(sum(spectrum,1))` (F3–F2), `squeeze(sum(spectrum,3))` (F2–F1), and `squeeze(sum(spectrum,2))` (F3–F1), each with 20 contours, `delta`, `k=2`, 256 colormap entries, `m=6`, and `signs`.

## Rendering and checks

The function uses a 2-by-2 subplot layout; the first panel is created with `'replace'`. It makes the 3D axes square, boxed, gridded, reverses X/Y/Z directions, sets a camera position from the axis extents, and labels the axes. The projection panels are square and gridded with reversed X/Y directions. Finally, it sets the current figure position to `[100 100 2*default_width 2*default_height]`. It returns no value.

The local guards require a real numeric 3D cube; three-element offsets, sweeps, point counts, zero-fill sizes, and spin entries; finite positive integer point-count and zero-fill vectors; and a spectrum shape matching `parameters.zerofill`. Units must be a character string recognised as `'ppm'`, `'Hz'`, or `'Gauss'`; `nsurf` and `k` must be positive integers and `delta` must contain four real values from 0 to 1.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/plot_3d.m)
- [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=plot_3d.m)
