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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 27-28: Check consistency; implemented by `grumble(x,npoints,order)`.
- Lines 30-31: Preallocate the answer; implemented by `dx=zeros(size(x)); x=x(:)`.
- Lines 33-34: Compute edges with sided schemes; implemented by `for n=1:(npoints-1)/2`.
- Lines 40-41: Fill in the middle with centered schemes; implemented by `stencil=((1-npoints)/2):((npoints-1)/2)`.

### Control flow inferred from the code

- Line 34: `for` loop over `n=1:(npoints-1)/2`.
- Line 43: `for` loop over `n=((npoints-1)/2+1):(numel(x)-(npoints-1)/2)`.

### Key state/data transformations

- Lines 31: computes `dx` using `dx=zeros(size(x)); x=x(:)`.
- Lines 35: computes `w` using `w=fdweights(n,1:npoints,order)`.
- Lines 36: computes `dx(n)` using `dx(n)=w(end,:)*x(1:npoints)`.
- Lines 37: computes `dx(end-n+1)` using `dx(end-n+1)=((-1)^order)*w(end,end:-1:1)*x((end-npoints+1):end)`.
- Lines 41: computes `stencil` using `stencil=((1-npoints)/2):((npoints-1)/2)`.

### Local helper functions

- Line 50: `grumble()` — `function grumble(x,npoints,order)`.
  - Representative operation: `if (~isnumeric(x))||(~isvector(x))`.
  - Representative operation: `error('x argument must be a vector.')`.

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
