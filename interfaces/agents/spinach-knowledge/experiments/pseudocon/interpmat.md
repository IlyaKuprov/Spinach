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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(cube_dims,ranges,xyz)`.
- Lines 35-36: Extract the bounds; implemented by `xmin=ranges(1); xmax=ranges(2)`.
- Lines 40-41: Compute grids; implemented by `x_grid=linspace(xmin,xmax,cube_dims(1))`.
- Lines 45-46: Compute grid intervals; implemented by `x_grid_interval=(xmax-xmin)/(cube_dims(1)-1)`.
- Lines 50-51: Preallocate return arrays; implemented by `rows=cell(size(xyz,1),1)`.
- Lines 55-56: Loop over points; implemented by `parfor n=1:size(xyz,1)`.
- Lines 58-59: Move into fractional coordinates; implemented by `x=(xyz(n,1)-xmin)/x_grid_interval`.
- Lines 63-67: Decide stencil points; implemented by `x_stencil=[min([max([1 (ceil(x)-1)]) (cube_dims(1)-3)]) min([max([2 (ceil(x)-0)]) (cube_dims(1)-2)]) min([max([3 (ceil(x)+1)]) (cube_dims(1)-1)]) min([max([4 (ceil(x)+2)…`.
- Lines 77-78: Extract subgrid values; implemented by `x_subgrid=x_grid(x_stencil); y_subgrid=y_grid(y_stencil); z_subgrid=z_grid(z_stencil)`.
- Lines 80-81: Compute interpolation vectors; implemented by `x_intvec=spalloc(1,cube_dims(1),4); x_intvec(x_stencil)=fdweights(xyz(n,1),x_subgrid,0)`.
- Lines 85-87: Store non-zeroes; implemented by `[~,cols{n},vals{n}]=find(kron(z_intvec, kron(y_intvec,x_intvec)))`.
- Lines 92-93: Merge cells; implemented by `rows=cell2mat(rows); cols=cell2mat(cols); vals=cell2mat(vals)`.
- Lines 95-96: Form the matrix; implemented by `P=sparse(rows,cols,vals,size(xyz,1),prod(cube_dims))`.

### Control flow inferred from the code

- Line 56: `parfor` loop over `n=1:size(xyz,1)`.

### Key state/data transformations

- Lines 36: computes `xmin` using `xmin=ranges(1); xmax=ranges(2)`.
- Lines 37: computes `ymin` using `ymin=ranges(3); ymax=ranges(4)`.
- Lines 38: computes `zmin` using `zmin=ranges(5); zmax=ranges(6)`.
- Lines 41: computes `x_grid` using `x_grid=linspace(xmin,xmax,cube_dims(1))`.
- Lines 42: computes `y_grid` using `y_grid=linspace(ymin,ymax,cube_dims(2))`.
- Lines 43: computes `z_grid` using `z_grid=linspace(zmin,zmax,cube_dims(3))`.
- Lines 46: computes `x_grid_interval` using `x_grid_interval=(xmax-xmin)/(cube_dims(1)-1)`.
- Lines 47: computes `y_grid_interval` using `y_grid_interval=(ymax-ymin)/(cube_dims(2)-1)`.
- Lines 48: computes `z_grid_interval` using `z_grid_interval=(zmax-zmin)/(cube_dims(3)-1)`.
- Lines 51: computes `rows` using `rows=cell(size(xyz,1),1)`.
- Lines 52: computes `cols` using `cols=cell(size(xyz,1),1)`.
- Lines 53: computes `vals` using `vals=cell(size(xyz,1),1)`.
- Lines 59: computes `x` using `x=(xyz(n,1)-xmin)/x_grid_interval`.
- Lines 60: computes `y` using `y=(xyz(n,2)-ymin)/y_grid_interval`.
- Lines 61: computes `z` using `z=(xyz(n,3)-zmin)/z_grid_interval`.
- Lines 64-67: computes `x_stencil` using `x_stencil=[min([max([1 (ceil(x)-1)]) (cube_dims(1)-3)]) min([max([2 (ceil(x)-0)]) (cube_dims(1)-2)]) min([max([3 (ceil(x)+1)]) (cube_dims(1)-1)]) min([max([4 (ceil(x)+2)…`.
- Lines 68-71: computes `y_stencil` using `y_stencil=[min([max([1 (ceil(y)-1)]) (cube_dims(2)-3)]) min([max([2 (ceil(y)-0)]) (cube_dims(2)-2)]) min([max([3 (ceil(y)+1)]) (cube_dims(2)-1)]) min([max([4 (ceil(y)+2)…`.
- Lines 72-75: computes `z_stencil` using `z_stencil=[min([max([1 (ceil(z)-1)]) (cube_dims(3)-3)]) min([max([2 (ceil(z)-0)]) (cube_dims(3)-2)]) min([max([3 (ceil(z)+1)]) (cube_dims(3)-1)]) min([max([4 (ceil(z)+2)…`.

### Local helper functions

- Line 101: `grumble()` — `function grumble(cube_dims,ranges,xyz)`.
  - Representative operation: `if (~isnumeric(cube_dims))||(~isreal(cube_dims))|| (numel(cube_dims)~=3)||any(mod(cube_dims,1)~=0)`.
  - Representative operation: `(numel(cube_dims)~=3)||any(mod(cube_dims,1)~=0)`.

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
