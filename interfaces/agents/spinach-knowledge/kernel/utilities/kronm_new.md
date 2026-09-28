# kernel/utilities/kronm_new.m

- Signature: `M=kronm_new(Q,M)`

## Purpose

Apply the Kronecker product of the matrices in `Q` to `M` without constructing the full Kronecker-product matrix.

## Physical / mathematical content

This is a matrix-free implementation of a tensor-product linear transformation.

## Numerical / algorithmic content

The input matrix is reshaped into the column dimensions of the factors, with its column count as a trailing dimension. The routine contracts each factor with the corresponding tensor dimension using `tensorprod`, then flattens the resulting row dimensions. The factor order corresponds to `Q{1} kron Q{2} kron ... kron Q{n}` acting on each column of `M`.

## Parameters / inputs

- `Q` - cell array of matrix factors.
- `M` - numeric vector or matrix with a row dimension compatible with the Kronecker-product factors.

## Outputs

- `M` - full vector or matrix after applying the Kronecker-product operator; one output column is returned for each input column.

## Implementation structure

The factors' row and column sizes define the tensor dimensions. Each factor is contracted with its assigned dimension, and the resulting tensor is reshaped to the output matrix.
