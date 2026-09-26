# experiments/pseudocon/interpmat.m

- Signature: `P=interpmat(cube_dims,ranges,xyz)`

## Purpose

Returns a matrix that acts on a stretched pseudocontact shift density cube and projects out the values of the PCS at the Cartesian coordinates given. Tricubic interpolation is used. Syntax: P=interpmat(cube_dims,ranges,xyz)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Parameters / inputs

- cube_dims -pseudocontact shift cube grid sizes, a vector of
- three integers ordered as [X Y Z]
- ranges -cartesian axis extents for the pseudocontact shift
- cube as [xmin xmax ymin ymax zmin zmax] in Angstroms.
- xyz -nuclear coordinates as [x y z] with multiple rows) at
- which PCS is to be evaluated, in Angstroms.
- Output:
- P -matrix projecting out PCS values at the specified
- nuclear positions from the stretched PCS cube.
- Note: this function is a part of the PCS inverse problem solver module; it
- should not normally be called directly by the user.

## Implementation structure

- Returns a matrix that acts on a stretched pseudocontact shift density cube
- and projects out the values of the PCS at the Cartesian coordinates given.
- Tricubic interpolation is used. Syntax:
- P=interpmat(cube_dims,ranges,xyz)
- cube_dims -pseudocontact shift cube grid sizes, a vector of
- three integers ordered as [X Y Z]
- ranges -cartesian axis extents for the pseudocontact shift
- cube as [xmin xmax ymin ymax zmin zmax] in Angstroms.
- xyz -nuclear coordinates as [x y z] with multiple rows) at
- which PCS is to be evaluated, in Angstroms.
- Output:
- P -matrix projecting out PCS values at the specified
