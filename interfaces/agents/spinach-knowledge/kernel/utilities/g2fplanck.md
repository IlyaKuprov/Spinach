# kernel/utilities/g2fplanck.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/g2fplanck.m`
- Signature: `G=g2fplanck(spin_system,parameters)`
- Total lines: 146

## Purpose

Returns gradient operators within the Fokker-Planck formalism used in the imaging module of Spinach. Syntax: G=g2fplanck(spin_system,parameters)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 40-41: Check consistency; implemented by `grumble(spin_system,parameters)`.
- Lines 43-44: Call Spinach to build magnet Zeeman Hamiltonian (per Tesla); implemented by `H=hamiltonian(assume(spin_system,'labframe','zeeman'))/spin_system.inter.magnet`.
- Lines 46-47: Empty arrays as default; implemented by `Gx=[]; Gy=[]; Gz=[]`.
- Lines 49-50: Build gradient operators; implemented by `switch numel(parameters.npts)`.
- Lines 54-55: Generate normalized gradient operators in 1D; implemented by `Gx=linspace(-0.5,0.5,parameters.npts(1))`.
- Lines 61-62: Generate normalized gradient operators in 2D; implemented by `Gx=parameters.dims(1)*linspace(-0.5,0.5,parameters.npts(1))`.
- Lines 69-70: Rotate gradients if necessary; implemented by `if isfield(parameters,'grad_angles')`.
- Lines 80-81: Generate normalized gradient operators; implemented by `Gx=parameters.dims(1)*linspace(-0.5,0.5,parameters.npts(1))`.
- Lines 102-103: Complain and bomb out; implemented by `error('incorrect number of spatial dimensions.')`.
- Lines 107-108: Pack the operators; implemented by `G={Gx,Gy,Gz}`.

### Control flow inferred from the code

- Line 50: dispatches on `numel(parameters.npts)`; cases `1`, `2`, `3`.
- Line 70: conditional branch on `isfield(parameters,'grad_angles')`.
- Line 92: conditional branch on `isfield(parameters,'grad_angles')`.

### Key state/data transformations

- Lines 44: computes `H` using `H=hamiltonian(assume(spin_system,'labframe','zeeman'))/spin_system.inter.magnet`.
- Lines 47: computes `Gx` using `Gx=[]; Gy=[]; Gz=[]`.
- Lines 63: computes `Gy` using `Gy=parameters.dims(2)*linspace(-0.5,0.5,parameters.npts(2))`.
- Lines 71: computes `R` using `R=[ cos(parameters.grad_angles) sin(parameters.grad_angles)`.
- Lines 73: computes `Gx_new` using `Gx_new=Gx*R(1,1)+Gy*R(1,2)`.
- Lines 74: computes `Gy_new` using `Gy_new=Gx*R(2,1)+Gy*R(2,2)`.
- Lines 83: computes `Gz` using `Gz=parameters.dims(3)*linspace(-0.5,0.5,parameters.npts(3))`.
- Lines 96: computes `Gz_new` using `Gz_new=Gx*R(3,1)+Gy*R(3,2)+Gz*R(3,3)`.
- Lines 108: computes `G` using `G={Gx,Gy,Gz}`.

### Local helper functions

- Line 113: `grumble()` — `function grumble(spin_system,parameters)`.
  - Representative operation: `if spin_system.inter.magnet==0`.
  - Representative operation: `error('the primary magnet field must be non-zero.')`.

## Parameters / inputs

- parameters.dims -a vector with one, two or three
- elements giving the dimensions
- of the box, metres
- parameters.npts -a vector with one, two or three
- elements giving number of points
- in each dimension of the box

## Outputs

- G -a cell array with the three gradient operators
- ordered as {Gx,Gy,Gz}, normalised to 1 T/m, empty
- matrices for non-exitent dimensions
- Note: gradients are assumed to be linear and centered on the
- middle of the sample.
- Note: the direct product order is Z(x)Y(x)X(x)Spin, this cor-
- responds to a column-wise vectorization of a 3D array
- with dimensions ordered as [X Y Z].
- Note: polyadic objects are returned, use inflate() to get the
- corresponding sparse matrix.

## Implementation structure

- Returns gradient operators within the Fokker-Planck formalism
- used in the imaging module of Spinach. Syntax:
- G=g2fplanck(spin_system,parameters)
- parameters.dims -a vector with one, two or three
- elements giving the dimensions
- of the box, metres
- parameters.npts -a vector with one, two or three
- elements giving number of points
- in each dimension of the box
- G -a cell array with the three gradient operators
- ordered as {Gx,Gy,Gz}, normalised to 1 T/m, empty
- matrices for non-exitent dimensions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `hamiltonian()`, `assume()`, `spdiags()`, `polyadic()`, `opium()`, `isfield()`, `euler2dcm()`, `any()`.
