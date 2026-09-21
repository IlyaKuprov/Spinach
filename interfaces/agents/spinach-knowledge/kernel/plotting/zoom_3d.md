# kernel/plotting/zoom_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/plotting/zoom_3d.m`
- Signature: `[density,ext]=zoom_3d(density,ext,zoom_ranges)`
- Total lines: 87

## Purpose

Zooms a 3D data cube to the fractional limits specified by the user. Syntax: [density,ext]=zoom_3d(density,ext,zoom_ranges)

## Physical / mathematical content

- This file belongs to the `kernel` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ext()`, `zoom_ranges()`, `density()`, `oldxvals()`, `oldyvals()`, `oldzvals()`, `ndims()`, `any()`.
