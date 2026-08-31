# kernel/derivatives/fdweights.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/derivatives/fdweights.m`
- Signature: `w=fdweights(target_point,grid_points,max_order)`
- Total lines: 81

## Purpose

Calculates finite difference weights for numerical derivatives, including order 0, which amounts to interpolation. Syntax: w=fdweights(target_point,grid_points,max_order)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check consistency; implemented by `grumble(target_point,grid_points,max_order)`.
- Lines 32-33: Compute the weights; implemented by `n=length(grid_points); w=zeros(max_order+1,n)`.

### Control flow inferred from the code

- Line 35: `for` loop over `i=2:n`.
- Line 37: `for` loop over `j=1:(i-1)`.
- Line 39: conditional branch on `j==(i-1)`.

### Key state/data transformations

- Lines 33: computes `n` using `n=length(grid_points); w=zeros(max_order+1,n)`.
- Lines 34: computes `w1` using `w1=1; w4=grid_points(1)-target_point; w(1,1)=1`.
- Lines 36: computes `mn` using `mn=min(i,max_order+1); w2=1; w5=w4; w4=grid_points(i)-target_point`.
- Lines 38: computes `w3` using `w3=grid_points(i)-grid_points(j); w2=w2*w3`.
- Lines 40: computes `w(2:mn,i)` using `w(2:mn,i)=w1*((1:(mn-1))'.*w(1:(mn-1),i-1)-w5*w(2:mn,i-1))/w2`.
- Lines 41: computes `w(1,i)` using `w(1,i)=-w1*w5*w(1,i-1)/w2`.
- Lines 43: computes `w(2:mn,j)` using `w(2:mn,j)=(w4*w(2:mn,j)-(1:(mn-1))'.*w(1:(mn-1),j))/w3`.
- Lines 44: computes `w(1,j)` using `w(1,j)=w4*w(1,j)/w3`.

### Local helper functions

- Line 52: `grumble()` — `function grumble(target_point,grid_points,max_order)`.
  - Representative operation: `if (~isnumeric(target_point))||(~isnumeric(grid_points))||(~isnumeric(max_order))|| (~isreal(target_point))||(~isreal(grid_points))||(~isreal(max_order))`.
  - Representative operation: `(~isreal(target_point))||(~isreal(grid_points))||(~isreal(max_order))`.

## Parameters / inputs

- target_point -the point at which the derivative
- is required
- grid_points -the points at which the function
- is given
- max_order -maximum derivative order

## Outputs

- w -finite difference coefficient array
- with the coefficients for the succes-
- sive derivatives in rows

## Implementation structure

- Calculates finite difference weights for numerical derivatives,
- including order 0, which amounts to interpolation. Syntax:
- w=fdweights(target_point,grid_points,max_order)
- target_point -the point at which the derivative
- is required
- grid_points -the points at which the function
- is given
- max_order -maximum derivative order
- w -finite difference coefficient array
- with the coefficients for the succes-
- sive derivatives in rows
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `grid_points()`, `isvector()`, `issorted()`.
