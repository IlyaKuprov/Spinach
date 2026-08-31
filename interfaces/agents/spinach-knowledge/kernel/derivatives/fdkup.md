# kernel/derivatives/fdkup.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fdkup.m`
- Signature: `K=fdkup(npoints,extents,chi,nstenc)`
- Total lines: 98

## Purpose

Returns a finite difference representation of the Kuprov operator: K[rho]=-(1/3)*Trace(Hessian[rho]*chi) with the number of stencil points in the finite difference approxi- mation specified by user. The resulting operator is a sparse matrix designed to act on the vectorisation of rho. The dimensions of rho are assumed to be ordered as [X Y Z]. For further information, see K=fdkup(npoints,extents,chi,nstenc)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 44-45: Check consistency; implemented by `grumble(npoints,extents,chi,nstenc)`.
- Lines 47-48: Build second derivative operators; implemented by `d2_dzdz=kron(kron(fdmat(npoints(3),nstenc,2),speye(npoints(2))),speye(npoints(1)))`.
- Lines 56-57: Normalise second derivative operators; implemented by `d2_dxdx=(npoints(1)/extents(1))*(npoints(1)/extents(1))*d2_dxdx`.
- Lines 67-70: Form the Kuprov operator; implemented by `K=-(1/3)*(chi(1,1)*d2_dxdx+chi(1,2)*d2_dxdy+chi(1,3)*d2_dxdz+ chi(2,1)*d2_dydx+chi(2,2)*d2_dydy+chi(2,3)*d2_dydz+ chi(3,1)*d2_dzdx+chi(3,2)*d2_dzdy+chi(3,3)*d2_dzdz)`.

### Key state/data transformations

- Lines 48: computes `d2_dzdz` using `d2_dzdz=kron(kron(fdmat(npoints(3),nstenc,2),speye(npoints(2))),speye(npoints(1)))`.
- Lines 49: computes `d2_dydy` using `d2_dydy=kron(kron(speye(npoints(3)),fdmat(npoints(2),nstenc,2)),speye(npoints(1)))`.
- Lines 50: computes `d2_dxdx` using `d2_dxdx=kron(kron(speye(npoints(3)),speye(npoints(2))),fdmat(npoints(1),nstenc,2))`.
- Lines 51: computes `d2_dydx` using `d2_dydx=kron(kron(speye(npoints(3)),fdmat(npoints(2),nstenc,1)),fdmat(npoints(1),nstenc,1))`.
- Lines 52: computes `d2_dzdx` using `d2_dzdx=kron(kron(fdmat(npoints(3),nstenc,1),speye(npoints(2))),fdmat(npoints(1),nstenc,1))`.
- Lines 53: computes `d2_dzdy` using `d2_dzdy=kron(kron(fdmat(npoints(3),nstenc,1),fdmat(npoints(2),nstenc,1)),speye(npoints(1)))`.
- Lines 54: computes `d2_dxdz` using `d2_dxdz=d2_dzdx; d2_dydz=d2_dzdy; d2_dxdy=d2_dydx`.
- Lines 58: computes `d2_dxdy` using `d2_dxdy=(npoints(1)/extents(1))*(npoints(2)/extents(2))*d2_dxdy`.
- Lines 62: computes `d2_dydz` using `d2_dydz=(npoints(2)/extents(2))*(npoints(3)/extents(3))*d2_dydz`.
- Lines 68-70: computes `K` using `K=-(1/3)*(chi(1,1)*d2_dxdx+chi(1,2)*d2_dxdy+chi(1,3)*d2_dxdz+ chi(2,1)*d2_dydx+chi(2,2)*d2_dydy+chi(2,3)*d2_dydz+ chi(3,1)*d2_dzdx+chi(3,2)*d2_dzdy+chi(3,3)*d2_dzdz)`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(npoints,extents,chi,nstenc)`.
  - Representative operation: `if (~isnumeric(npoints))||(numel(npoints)~=3)||(~isreal(npoints))|| (any(npoints<1))||any(mod(npoints,1)~=0)`.
  - Representative operation: `(any(npoints<1))||any(mod(npoints,1)~=0)`.

## Parameters / inputs

- npoints -a three-element vector specifying the dimensions
- of the 3D cube of data that the operator will be
- acting on, in Angstroms. The dimensions are assu-
- med to be ordered as [X Y Z].
- chi -the electron magnetic susceptibility tensor in
- cubic Angstroms, a symmetric 3x3 matrix.
- extents -a three-element vector specifying axis extents
- in Angstroms. The dimensions are assumed to be
- ordered as [X Y Z].
- nstenc -number of finite-difference stencil points for
- the finite-difference approximation. Periodic
- boundary conditions are used.

## Outputs

- K -a sparse matrix designed to act on the vectori-
- zation of the array. The dimensions are assumed
- to be ordered as [X Y Z].

## Implementation structure

- Returns a finite difference representation of the Kuprov operator:
- K[rho]=-(1/3)*Trace(Hessian[rho]*chi)
- with the number of stencil points in the finite difference approxi-
- mation specified by user. The resulting operator is a sparse matrix
- designed to act on the vectorisation of rho. The dimensions of rho
- are assumed to be ordered as [X Y Z]. For further information, see
- K=fdkup(npoints,extents,chi,nstenc)
- npoints - a three-element vector specifying the dimensions
- of the 3D cube of data that the operator will be
- acting on, in Angstroms. The dimensions are assu-
- med to be ordered as [X Y Z].
- chi - the electron magnetic susceptibility tensor in

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fdmat()`, `npoints()`, `speye()`, `extents()`, `chi()`, `any()`.
