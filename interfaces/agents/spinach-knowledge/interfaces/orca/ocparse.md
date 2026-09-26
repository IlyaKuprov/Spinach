# interfaces/orca/ocparse.m

- Signature: `[density,ext,dx,dy,dz]=ocparse(filename,pad_factor)`

## Purpose

ORCA cube file parser. Extracts the normalised probability density and the associated metric information from ORCA spin density in "3D simple format" (see ORCA manual). Syntax: [density,ext,dx,dy,dz]=ocparse(filename,pad_factor)

## Physical / mathematical content

- ORCA interfaces. They recover quantum-chemistry tensors and metadata and convert them to Spinach conventions.

## Numerical / algorithmic content

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
