# kernel/plotting/volplot.m

- Signature: `volplot(data_cube,axis_ranges,clip_ranges)`

## Purpose

Volumetric 3D plot of a scalar field. Sign is mapped to colour and amplitude to opacity, with separate scaling for positive and negative values. Displaying the colour bar is recommended.

## Parameters / inputs

- `data_cube` — real 3D numeric array ordered as `[X Y Z]`, with at least two points along each dimension.
- `axis_ranges` — real six-element vector `[xmin xmax ymin ymax zmin zmax]`; each minimum must be less than its maximum. Defaults to `[-1 1 -1 1 -1 1]` if omitted.
- `clip_ranges` — optional real two-element vector of positive and negative clipping fractions, respectively. Each element must be in `(0,1]`; the default is `[1 1]`. Clipping can be useful for steep functions.

## Outputs

The function produces a figure; it has no return value.

## Numerical / algorithmic content

Positive values are divided by their positive maximum and mapped into `[0,1]`; negative values are scaled by the magnitude of their negative minimum and mapped into `[-1,0]`. If a corresponding clipping fraction is less than one, values beyond that fraction are clipped and the remaining values remapped to the same interval. The function reports these scaling and clipping operations in the command window.

The array is permuted to match the `surf`/`meshgrid` convention. The figure is drawn from surfaces parallel to the XY, XZ, and YZ planes. Values with magnitude below `1/64` are omitted from each plane. Surface colour represents the signed value, while alpha data uses its absolute value.

The plot uses a blue–white–red colormap and fixes colour limits at `[-1 1]`. Its default alpha map is interpolated, divided by five, and set to zero where the resulting opacity is below `0.01`. Axis extents come from `axis_ranges`, independently of the plotted geometry; the plot uses equal data aspect ratios, perspective projection, a grid, and X, Y, and Z labels.

## Links

- [volplot.m documentation](https://spindynamics.org/wiki/index.php?title=volplot.m)
- [Pseudocontact-shift paper linked in the source](http://dx.doi.org/10.1039/C4CP03106G)