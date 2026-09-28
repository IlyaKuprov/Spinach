# examples/fundamentals/tensor_structures/kronm_test_1.m

- Signature: `kronm_test_1()`

## Purpose

Checks that `kronm` applies a list of Kronecker factors to a matrix of input columns consistently with forming the full Kronecker-product matrix and multiplying it by the same input. The test covers randomly sized real and complex factors and inputs.

## Physical / mathematical content

- For factors `Q_terms{1},...,Q_terms{nmats}`, the reference matrix is their Kronecker product `Q`; the tested identity is `kronm(Q_terms,x) = Q*x`.
- Each case uses the same number of rows in `x` as the product of the factor dimensions, with a randomly selected number of columns.

## Numerical / algorithmic content

- The number of factors is randomly selected from 3 to 6; each factor dimension is 2 to 4, and the input has 1 to 20 columns.
- Separate real and complex comparisons use the one-norm difference threshold `1e-6`. Both `kronm` and dense `Q*x` timings are printed.

## Implementation structure

- Generate the dimensions, build real factors and their full Kronecker product, and compare `kronm(Q_terms,x)` with `Q*x`.
- Repeat with complex factors and a complex input; each comparison reports a pass or raises an error.
