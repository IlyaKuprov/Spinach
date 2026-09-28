# kernel/utilities/jacobianest.m

- Signature: `[jac,err] = jacobianest(fun,x0)`

## Purpose

Estimate the Jacobian of a vector-valued function at `x0`, together with an entry-wise error estimate.

## Physical / mathematical content

A general numerical differentiation utility; it does not encode a spin-system model.

## Numerical / algorithmic content

For each element of `x0`, the routine evaluates centered finite differences over a geometrically decreasing sequence of 26 step sizes. Romberg extrapolation cancels the leading second- and fourth-order error terms; after trimming the three largest- and smallest-step estimates, the estimate with the smallest predicted error is selected.

## Parameters / inputs

- `fun` - function handle accepting the vector or array `x0` and returning a vector-valued result.
- `x0` - numeric vector or array; each of its `numel(x0)` elements is treated as an independent variable.

## Outputs

- `jac` - Jacobian array with one row per element of `fun(x0)` after linearisation and one column per element of `x0`.
- `err` - estimated error for each corresponding Jacobian entry; it has the same size as `jac`.

## Implementation structure

The result at the centre point determines the output dimension. Each input coordinate is perturbed in both directions, the finite-difference estimates are extrapolated, and the selected derivative and error estimate are stored in the corresponding Jacobian column.
