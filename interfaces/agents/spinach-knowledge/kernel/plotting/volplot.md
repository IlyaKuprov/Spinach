# kernel/plotting/volplot.m

- Source: [kernel/plotting/volplot.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/plotting/volplot.m)
- Wiki: [volplot.m](https://spindynamics.org/wiki/index.php?title=volplot.m)
- Paper cited in the source: [DOI 10.1039/C4CP03106G](http://dx.doi.org/10.1039/C4CP03106G)

## Purpose

Render a real scalar field as a volumetric 3D plot: sign determines colour and magnitude determines opacity. The source recommends displaying a colour bar, but the function does not add one.

## Inputs

- `data_cube` is a real numeric 3D array ordered as `[X Y Z]`, with at least two samples along each dimension.
- `axis_ranges` is `[xmin xmax ymin ymax zmin zmax]`; each lower bound must be less than its upper bound. It defaults to `[-1 1 -1 1 -1 1]`.
- `clip_ranges` gives the positive and negative clipping fractions, respectively. It defaults to `[1 1]`; both values must lie in `(0,1]`.

## Scaling and rendering

Positive and negative data are normalised independently: the positive maximum maps to `1`, and the magnitude of the negative minimum maps to `-1`. A clipping fraction below 1 caps that sign at the requested fraction and remaps the retained range to the full sign interval. The source permutes the cube for surface plotting, draws orthogonal stacks of `surf` planes, and uses absolute scaled values for opacity. Values with magnitude below `1/64` in a plotted plane are changed to `NaN` so that they are not rendered. The blue-white-red colour map represents sign; the colour scale is fixed at `[-1 1]`.

## Figure effects

The function clears the current figure with `clf`, draws into its axes, and changes the figure colour map and alpha map. It sets the requested axis extents independently of plotted geometry, equal data aspect, perspective projection, grid, and X/Y/Z labels, then turns hold off. It returns no value and writes no file; it does not create the recommended colour bar.
