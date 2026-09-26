# kernel/derivatives/fdhess.m

- Signature: `H=fdhess(A,nstenc)`

## Purpose

Returns the finite-difference Hessian of a 3D array using a finite difference scheme with a user-specified number of stencil points and a unit grid spacing. The dimensions of the 3D array are assumed to be ordered as [X Y Z]. Syntax: H=fdhess(A,nstenc)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

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
