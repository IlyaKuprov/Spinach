# kernel/utilities/keep_rank.m

- Signature: `A=keep_rank(A,nsvk)`

## Purpose

Truncate a matrix to a requested singular-value rank and return the reconstructed matrix.

## Physical / mathematical content

This is a numerical low-rank approximation based on the singular value decomposition; it does not assume a particular physical model.

## Numerical / algorithmic content

The input is converted to a full matrix and factorised with `svd`. The routine retains the first `nsvk` singular components and forms `U(:,1:nsvk)*S(1:nsvk,1:nsvk)*V(:,1:nsvk)'`.

## Parameters / inputs

- `A` - numeric matrix with more than one row and more than one column; sparse input is converted to full.
- `nsvk` - positive real integer no greater than the smaller matrix dimension.

## Outputs

- `A` - full matrix reconstructed from the retained singular components.

## Implementation structure

Consistency checks precede the full SVD. The requested leading singular-vector columns and matching diagonal block of `S` are multiplied to produce the truncated matrix.
