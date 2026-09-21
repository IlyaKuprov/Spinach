# kernel/derivatives/fdvec.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fdvec.m`
- Signature: `dx=fdvec(x,npoints,order)`
- Total lines: 85

## Purpose

Performs arbitrary-order finite-difference differentiation of a user-supplied row or column vector. Uses central finite-differe- nce stencils in the middle and sided stencils of the same order of accuracy on the sides. Syntax: dx=fdvec(x,npoints,order)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- x -column or row vector to be differentiated
- npoints -number of points in the finite difference
- stencil
- order -order of the derivative required

## Outputs

- dx -column or row vector with the derivative

## Implementation structure

- Performs arbitrary-order finite-difference differentiation of a
- user-supplied row or column vector. Uses central finite-differe-
- nce stencils in the middle and sided stencils of the same order
- of accuracy on the sides. Syntax:
- dx=fdvec(x,npoints,order)
- x -column or row vector to be differentiated
- npoints -number of points in the finite difference
- stencil
- order -order of the derivative required
- dx -column or row vector with the derivative
- Check consistency
- Preallocate the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fdweights()`, `isvector()`.
