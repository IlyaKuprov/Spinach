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
