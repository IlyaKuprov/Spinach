# experiments/pseudocon/interpmat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/pseudocon/interpmat.m`
- Signature: `P=interpmat(cube_dims,ranges,xyz)`
- Total lines: 129

## Purpose

Returns a matrix that acts on a stretched pseudocontact shift density cube and projects out the values of the PCS at the Cartesian coordinates given. Tricubic interpolation is used. Syntax: P=interpmat(cube_dims,ranges,xyz)

## Physical / mathematical content

- Paramagnetic-pseudocontact inference routines. The mathematics includes inverse problems, tensor parameterisation, interpolation, and regularisation.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.
- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `ranges()`, `cube_dims()`, `xyz()`, `x_grid()`, `y_grid()`, `z_grid()`, `spalloc()`, `x_intvec()`, `fdweights()`, `y_intvec()`, `z_intvec()`, `cell2mat()`, `any()`.
