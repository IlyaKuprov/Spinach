# interfaces/orca/ocparse.m

- Signature: `[density,ext,dx,dy,dz]=ocparse(filename,pad_factor)`

## Purpose

Parses an ORCA cube file containing spin density in “3D simple format” and returns the normalised probability density and associated grid metrics.

## Inputs

- `filename`: path to an existing cube file.
- `pad_factor`: nonnegative integer scalar controlling zero-padding on each side of the grid.

## Processing and outputs

The parser takes the absolute density values, arranges the density as `[X Y Z]`, normalises it by numerical integration, and zero-pads the grid according to `pad_factor`. It returns:

- `density`: padded, normalised density.
- `ext`: updated extents `[xmin xmax ymin ymax zmin zmax]` in Angstrom.
- `dx`, `dy`, `dz`: grid steps in Angstrom.

## Source

Ilya Kuprov, Elizaveta Suturina, and Petra Pikulova. [ocparse.m on the Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=ocparse.m)
