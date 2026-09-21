# kernel/derivatives/fdmat.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fdmat.m`
- Signature: `D=fdmat(dim,nstenc,order,boundary)`
- Total lines: 100

## Purpose

Returns arbitrary-order central finite-difference differentiation matrices (sparse) with unit grid point spacing. Syntax: D=fdmat(dim,nstenc,order,boundary)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- dim -dimension of the column vector to be
- differentiated
- nstenc -number of points in the finite diffe-
- rence stencil
- order -order of the derivative required
- boundary -'wall' fills the edges with sided
- finite difference schemes, 'pbc'
- assumes periodic boundaries. The
- default is 'pbc'.

## Outputs

- D -finite difference differentiation matrix

## Implementation structure

- Returns arbitrary-order central finite-difference differentiation
- matrices (sparse) with unit grid point spacing. Syntax:
- D=fdmat(dim,nstenc,order,boundary)
- dim -dimension of the column vector to be
- differentiated
- nstenc -number of points in the finite diffe-
- rence stencil
- order -order of the derivative required
- boundary -'wall' fills the edges with sided
- finite difference schemes, 'pbc'
- assumes periodic boundaries. The
- default is 'pbc'.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `grumble()`, `spalloc()`, `fdweights()`, `ischar()`.
