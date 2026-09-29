# experiments/pseudocon/centroid.m

Source: [canonical MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/experiments/pseudocon/centroid.m) · [Spinach wiki](https://spindynamics.org/wiki/index.php?title=centroid.m)

Signature: `[x,y,z]=centroid(probden,ranges)`.

## Behaviour

This helper finds the coordinate centroid of a three-dimensional real numeric array. It creates coordinate grids with `ndgrid` from three inclusive `linspace` axes spanning `[xmin,xmax]`, `[ymin,ymax]`, and `[zmin,zmax]`, with each axis length taken from the corresponding dimension of `probden`. It then returns three ratios of nested trapezoidal integrals: each coordinate multiplied by `probden`, divided by the nested trapezoidal integral of `probden`. The resulting scalar coordinates use the units of the supplied ranges; the source specifies no physical coordinate unit.

Despite its pseudocontext folder, this function does not calculate a pseudocontact shift or tensor. It accepts no susceptibility tensor, magnetic-field orientation, nucleus, or spin-system argument; it returns only the centroid of the input array.

## Inputs, outputs, and limits

`probden` is described as a probability-density cube ordered `[X Y Z]`. The code checks that it is a numeric, real, three-dimensional array; it does not check non-negativity or unit normalisation. `ranges` must be a real numeric six-element vector ordered `[xmin xmax ymin ymax zmin zmax]`, with each lower bound strictly less than its upper bound. Outputs `x`, `y`, and `z` are scalar coordinates.

The integrals use nested default `trapz` calls without explicit spacing arguments. The source has no special handling for a zero normalisation, so a nonzero integrated density is required for finite centroid coordinates. No runtime result or accuracy assessment is implied.
