# kernel/overloads/@ttclass/amensum.m

- Signature: `y=amensum(x,tol,opts)`

## Purpose

Sums buffered tensor trains in a single tensor train using AMEn algorithm. Syntax: y=amensum(x,tol,opts)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- x -ttclass with buffered rank-one tensors
- tol -relative tolerance parameter, e.g. 1e-10
- The opts field is optional:
- opts.max_swp -maximum number of iterations
- opts.init_guess_rank -rank of the initial guess
- opts.enrichment_rank -rank of the enrichment
- opts.verb -verbosity switch

## Outputs

- y -ttclass with a single tensor train, such
- that |x-y|<tol*|x| in Frobenius norm

## Implementation structure

- Sums buffered tensor trains in a single tensor train using
- AMEn algorithm. Syntax:
- y=amensum(x,tol,opts)
- x -ttclass with buffered rank-one tensors
- tol -relative tolerance parameter, e.g. 1e-10
- The opts field is optional:
- opts.max_swp -maximum number of iterations
- opts.init_guess_rank -rank of the initial guess
- opts.enrichment_rank -rank of the enrichment
- opts.verb -verbosity switch
- y -ttclass with a single tensor train, such
- that |x-y|<tol*|x| in Frobenius norm
