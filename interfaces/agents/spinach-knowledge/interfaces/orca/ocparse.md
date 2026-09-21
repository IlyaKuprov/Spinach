# interfaces/orca/ocparse.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/orca/ocparse.m`
- Signature: `[density,ext,dx,dy,dz]=ocparse(filename,pad_factor)`
- Total lines: 95

## Purpose

ORCA cube file parser. Extracts the normalised probability density and the associated metric information from ORCA spin density in "3D simple format" (see ORCA manual). Syntax: [density,ext,dx,dy,dz]=ocparse(filename,pad_factor)

## Physical / mathematical content

- ORCA interfaces. They recover quantum-chemistry tensors and metadata and convert them to Spinach conventions.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- filename -character string specifying the file to load
- pad_factor -padding factor specifying how many multiples
- the array dimension in zeros to add on each
- side of the cube

## Outputs

- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]
- dx,dy,dz -grid steps in the three directions, Angstrom

## Implementation structure

- ORCA cube file parser. Extracts the normalised probability density and
- the associated metric information from ORCA spin density in "3D simple
- format" (see ORCA manual). Syntax:
- [density,ext,dx,dy,dz]=ocparse(filename,pad_factor)
- filename -character string specifying the file to load
- pad_factor -padding factor specifying how many multiples
- the array dimension in zeros to add on each
- side of the cube
- density -probability density cube with dimensions
- ordered as [X Y Z]
- ext -grid extents in Angstrom, ordered as
- [xmin xmax ymin ymax zmin zmax]

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `importdata()`, `str2num()`, `npts()`, `dxdydz()`, `trapz()`, `corner_xyz()`, `padarray()`, `ext()`, `ischar()`, `exist()`, `isscalar()`.
