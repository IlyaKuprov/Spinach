# kernel/utilities/kronm.m

- Signature: `x=kronm(Q,x)`

## Purpose

Apply the Kronecker product of the matrices in `Q` to `x` without constructing the full Kronecker-product matrix.

## Physical / mathematical content

This is a matrix-free implementation of a tensor-product linear transformation.

## Numerical / algorithmic content

The input is reshaped into dimensions matching the column dimensions of the factors, with any columns of `x` carried as a trailing dimension. The routine applies each factor along its corresponding dimension, using reshape and permutation operations as needed, then flattens the row dimensions. The factor order corresponds to `Q{1} kron Q{2} kron ... kron Q{n}` acting on each column of `x`.

## Parameters / inputs

- `Q` - cell array of matrix factors.
- `x` - numeric vector or matrix with a row dimension compatible with the Kronecker-product factors.

## Outputs

- `x` - full vector or matrix after applying the Kronecker-product operator; one output column is returned for each input column.

## Implementation structure

The factors' row and column sizes define the tensor dimensions. Each factor is multiplied against its assigned dimension, with an explicit permutation for non-leading dimensions; the final tensor is reshaped into the output matrix.
