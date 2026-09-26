# kernel/plotting/zoom_3d.m

- Signature: `[density,ext]=zoom_3d(density,ext,zoom_ranges)`

## Purpose

Zooms a 3D data cube to the fractional limits specified by the user. Syntax: [density,ext]=zoom_3d(density,ext,zoom_ranges)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]
- zoom_ranges -zoom ranges along each axis as fractions,
- ordered as [xmin xmax ymin ymax zmin zmax],
- e.g. [0.3 0.6 0.1 0.2 0.5 0.8]

## Outputs

- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]

## Implementation structure

- Zooms a 3D data cube to the fractional limits specified
- by the user. Syntax:
- [density,ext]=zoom_3d(density,ext,zoom_ranges)
- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]
- zoom_ranges -zoom ranges along each axis as fractions,
- ordered as [xmin xmax ymin ymax zmin zmax],
- e.g. [0.3 0.6 0.1 0.2 0.5 0.8]
- Check consistency
- Generate axis ticks
