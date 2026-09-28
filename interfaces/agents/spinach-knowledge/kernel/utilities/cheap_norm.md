# kernel/utilities/cheap_norm.m

- Signature: `n=cheap_norm(A,t,itmax)`

## Purpose

Returns the least expensive supported norm for the representation of `A`: the infinity norm for GPU arrays, the 1-norm for other non-polyadic arrays, and a lower-bound 1-norm estimate for polyadic objects.

## Physical / mathematical content

The function selects a norm calculation appropriate to the array representation. For polyadic objects, it estimates the 1-norm using probe vectors and the iterative method in Algorithm 2.4 of Higham and Tisseur's paper ([DOI](https://doi.org/10.1137/S0895479899356080)).

## Numerical / algorithmic content

GPU arrays return `norm(A,inf)`; non-polyadic CPU arrays return `norm(A,1)`. For polyadic inputs, the estimator uses matrix-vector and adjoint products to update a lower bound on the 1-norm. The number of probe columns `t` is limited by the number of columns of `A`.

## Parameters / inputs

- `A` — a matrix or polyadic representation
- `t` — optional number of probe columns for the polyadic estimator; defaults to 1
- `itmax` — optional maximum number of estimator iterations; defaults to 5

## Outputs

- `n` — infinity norm for GPU arrays, 1-norm for other non-polyadic arrays, or a lower-bound 1-norm estimate for polyadic objects
- The function uses the least expensive of these supported norm calculations.
