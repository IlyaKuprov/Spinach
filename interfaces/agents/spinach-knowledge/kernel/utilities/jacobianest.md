# kernel/utilities/jacobianest.m

- Signature: `[jac,err] = jacobianest(fun,x0)`

## Purpose

Estimate of the Jacobian matrix of a vector valued function of n variables. Syntax: [jac,err] = jacobianest(fun,x0)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- Finite-difference discretisation appears in the implementation, so numerical accuracy depends on stencil order, boundary handling, and the balance between resolution and conditioning.

## Parameters / inputs

- fun -(vector valued) analytical function to differentiate.
- fun must be a function of the vector or array x0.
- x0 -vector location at which to differentiate fun
- If x0 is an nxm array, then fun is assumed to be
- a function of n*m variables.

## Outputs

- jac -array of first partial derivatives of fun.
- Assuming that x0 is a vector of length p
- and fun returns a vector of length n, then
- jac will be an array of size (n,p)
- err -vector of error estimates corresponding to
- each partial derivative in jac.
- John D'Errico

## Implementation structure

- Estimate of the Jacobian matrix of a vector valued
- function of n variables. Syntax:
- [jac,err] = jacobianest(fun,x0)
- fun -(vector valued) analytical function to differentiate.
- fun must be a function of the vector or array x0.
- x0 -vector location at which to differentiate fun
- If x0 is an nxm array, then fun is assumed to be
- a function of n*m variables.
- jac -array of first partial derivatives of fun.
- Assuming that x0 is a vector of length p
- and fun returns a vector of length n, then
- jac will be an array of size (n,p)
