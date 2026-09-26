# kernel/utilities/tikhol1n.m

- Signature: `[x,err,reg]=tikhol1n(A,y,nnzt)`

## Purpose

L1 norm Tikhonov regularised solver for A*x=y where A is an ill-conditioned matrix. The error functional is norm(A*x-y,2)^2+lambda*norm(x,1), it is minimised using the FISTA algorithm. The user specifies the de- sired number of non-zeroes, lambda parameter is then found by bracketing / bisection. Syntax: [x,err,reg]=tikhol1n(A,y,nnzt)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Parameters / inputs

- A -a real or complex matrix
- y -a real or complex column vector
- nnzt -the target for the number of
- non-zeroes in the solution

## Outputs

- x -a real or complex vector
- err -squared 2-norm of the fitting
- error divided by the squared
- 2-norm of the solution
- reg -1-norm of the solution

## Implementation structure

- L1 norm Tikhonov regularised solver for A*x=y where
- A is an ill-conditioned matrix. The error functional
- is norm(A*x-y,2)^2+lambda*norm(x,1), it is minimised
- using the FISTA algorithm. The user specifies the de-
- sired number of non-zeroes, lambda parameter is then
- found by bracketing / bisection. Syntax:
- [x,err,reg]=tikhol1n(A,y,nnzt)
- A -a real or complex matrix
- y -a real or complex column vector
- nnzt -the target for the number of
- non-zeroes in the solution
- x -a real or complex vector
