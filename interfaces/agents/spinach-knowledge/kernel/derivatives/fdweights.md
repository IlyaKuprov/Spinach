# kernel/derivatives/fdweights.m

- Signature: `w=fdweights(target_point,grid_points,max_order)`

## Purpose

Calculates finite difference weights for numerical derivatives, including order 0, which amounts to interpolation. Syntax: w=fdweights(target_point,grid_points,max_order)

## Physical / mathematical content

- Derivative utilities. These routines compute finite-difference, analytical, or optimisation-oriented derivatives needed for sensitivity analysis, fitting, and optimal control.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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
