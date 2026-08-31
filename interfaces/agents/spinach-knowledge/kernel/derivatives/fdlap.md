# kernel/derivatives/fdlap.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fdlap.m`
- Signature: `L=fdlap(dims,extents,nstenc)`
- Total lines: 112

## Purpose

Returns a finite-difference representation of the Laplacian for an array with a user-specified finite difference stencil size. The re- sulting operator is a sparse matrix designed to act on the vectori- sation of the array. The dimensions of the array are assumed to be ordered as [X Y Z]. Syntax: L=fdlap(npoints,extents,nstenc)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check consistency; implemented by `grumble(dims,extents,nstenc)`.
- Lines 44-45: Get differentiation matrices; implemented by `Dxx=fdmat(dims(1),nstenc,2)`.
- Lines 47-48: Normalize differentiation matrices; implemented by `Dxx=(dims(1)/extents(1))^2*Dxx`.
- Lines 50-51: Compute the Laplacian; implemented by `L=Dxx`.
- Lines 63-65: Compute the Laplacian; implemented by `L=kron(Dyy,speye(dims(1)))+ kron(speye(dims(2)),Dxx)`.
- Lines 79-82: Compute the Laplacian; implemented by `L=kron(kron(Dzz,speye(dims(2))),speye(dims(1)))+ kron(kron(speye(dims(3)),Dyy),speye(dims(1)))+ kron(kron(speye(dims(3)),speye(dims(2))),Dxx)`.
- Lines 86-87: Complain and bomb out; implemented by `error('incorrect number of spatial dimensions.')`.

### Control flow inferred from the code

- Line 40: dispatches on `numel(dims)`; cases `1`, `2`, `3`.

### Key state/data transformations

- Lines 45: computes `Dxx` using `Dxx=fdmat(dims(1),nstenc,2)`.
- Lines 51: computes `L` using `L=Dxx`.
- Lines 57: computes `Dyy` using `Dyy=fdmat(dims(2),nstenc,2)`.
- Lines 72: computes `Dzz` using `Dzz=fdmat(dims(3),nstenc,2)`.

### Local helper functions

- Line 94: `grumble()` — `function grumble(dims,extents,nstenc)`.
  - Representative operation: `if (~isnumeric(dims))||(~isreal(dims))||(any(dims<1))||any(mod(dims,1)~=0)`.
  - Representative operation: `error('npoints must be a three-element vector of positive integers.')`.

## Parameters / inputs

- dims -a one-element, two-element, or three-element
- vector specifying the number of discretisation
- points in each dimension of the 1D, 2D, or 3D
- array of data that the operator will be acting
- on, ordered as [X Y Z].
- extents -a one-element, two-element, or three-element
- vector specifying the size of each dimension
- of the array, ordered as [X Y Z].
- nstenc -number of finite-difference stencil points for
- the finite-difference approximation; periodic
- boundary conditions are used

## Outputs

- L -a sparse matrix designed to act on the vectori-
- zation of the array. The dimensions are assumed
- to be ordered as [X Y Z].

## Implementation structure

- Returns a finite-difference representation of the Laplacian for an
- array with a user-specified finite difference stencil size. The re-
- sulting operator is a sparse matrix designed to act on the vectori-
- sation of the array. The dimensions of the array are assumed to be
- ordered as [X Y Z]. Syntax:
- L=fdlap(npoints,extents,nstenc)
- dims - a one-element, two-element, or three-element
- vector specifying the number of discretisation
- points in each dimension of the 1D, 2D, or 3D
- array of data that the operator will be acting
- on, ordered as [X Y Z].
- extents - a one-element, two-element, or three-element

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fdmat()`, `dims()`, `extents()`, `speye()`, `any()`.
