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
