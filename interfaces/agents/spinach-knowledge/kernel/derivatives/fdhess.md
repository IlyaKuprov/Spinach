# kernel/derivatives/fdhess.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fdhess.m`
- Signature: `H=fdhess(A,nstenc)`
- Total lines: 72

## Purpose

Returns the finite-difference Hessian of a 3D array using a finite difference scheme with a user-specified number of stencil points and a unit grid spacing. The dimensions of the 3D array are assumed to be ordered as [X Y Z]. Syntax: H=fdhess(A,nstenc)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(A,nstenc)`.
- Lines 33-34: Compute derivatives; implemented by `d2A_dzdz=reshape(kron(kron(fdmat(size(A,3),nstenc,2),speye(size(A,2))),speye(size(A,1)))*A(:),size(A))`.
- Lines 44-45: Form the Hessian array; implemented by `H={d2A_dxdx d2A_dxdy d2A_dxdz`.

### Key state/data transformations

- Lines 34: computes `d2A_dzdz` using `d2A_dzdz=reshape(kron(kron(fdmat(size(A,3),nstenc,2),speye(size(A,2))),speye(size(A,1)))*A(:),size(A))`.
- Lines 35: computes `d2A_dzdy` using `d2A_dzdy=reshape(kron(kron(fdmat(size(A,3),nstenc,1),fdmat(size(A,2),nstenc,1)),speye(size(A,1)))*A(:),size(A))`.
- Lines 36: computes `d2A_dzdx` using `d2A_dzdx=reshape(kron(kron(fdmat(size(A,3),nstenc,1),speye(size(A,2))),fdmat(size(A,1),nstenc,1))*A(:),size(A))`.
- Lines 37: computes `d2A_dydz` using `d2A_dydz=reshape(kron(kron(fdmat(size(A,3),nstenc,1),fdmat(size(A,2),nstenc,1)),speye(size(A,1)))*A(:),size(A))`.
- Lines 38: computes `d2A_dydy` using `d2A_dydy=reshape(kron(kron(speye(size(A,3)),fdmat(size(A,2),nstenc,2)),speye(size(A,1)))*A(:),size(A))`.
- Lines 39: computes `d2A_dydx` using `d2A_dydx=reshape(kron(kron(speye(size(A,3)),fdmat(size(A,2),nstenc,1)),fdmat(size(A,1),nstenc,1))*A(:),size(A))`.
- Lines 40: computes `d2A_dxdz` using `d2A_dxdz=reshape(kron(kron(fdmat(size(A,3),nstenc,1),speye(size(A,2))),fdmat(size(A,1),nstenc,1))*A(:),size(A))`.
- Lines 41: computes `d2A_dxdy` using `d2A_dxdy=reshape(kron(kron(speye(size(A,3)),fdmat(size(A,2),nstenc,1)),fdmat(size(A,1),nstenc,1))*A(:),size(A))`.
- Lines 42: computes `d2A_dxdx` using `d2A_dxdx=reshape(kron(kron(speye(size(A,3)),speye(size(A,2))),fdmat(size(A,1),nstenc,2))*A(:),size(A))`.
- Lines 45: computes `H` using `H={d2A_dxdx d2A_dxdy d2A_dxdz`.

### Local helper functions

- Line 52: `grumble()` — `function grumble(A,npoints)`.
  - Representative operation: `if (~isnumeric(A))||(ndims(A)~=3)`.
  - Representative operation: `error('A must be a three-dimensional numeric array.')`.

## Parameters / inputs

- A -a 3D array with dimensions ordered as [X Y Z]
- nstenc -number of points inthe finite difference stencil,
- periodic boundary conditions are used

## Outputs

- H -a 3x3 cell array of 3D arrays ordered in the
- following way:
- {d2A_dxdx d2A_dxdy d2A_dxdz
- d2A_dydx d2A_dydy d2A_dydz
- d2A_dzdx d2A_dzdy d2A_dzdz}

## Implementation structure

- Returns the finite-difference Hessian of a 3D array using a finite
- difference scheme with a user-specified number of stencil points and
- a unit grid spacing. The dimensions of the 3D array are assumed to
- be ordered as [X Y Z]. Syntax:
- H=fdhess(A,nstenc)
- A -a 3D array with dimensions ordered as [X Y Z]
- nstenc -number of points inthe finite difference stencil,
- periodic boundary conditions are used
- H -a 3x3 cell array of 3D arrays ordered in the
- following way:
- {d2A_dxdx d2A_dxdy d2A_dxdz
- d2A_dydx d2A_dydy d2A_dydz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fdmat()`, `speye()`, `ndims()`, `any()`.
