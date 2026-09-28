# kernel/overloads/@rcv/ctranspose.m

- Signature: `A=ctranspose(A)`

## Purpose

Returns the conjugate transpose of an RCV sparse matrix while retaining its RCV representation.

## Physical / mathematical content

For a complex matrix, the conjugate transpose swaps row and column indices and complex-conjugates each stored value.

## Numerical / algorithmic content

The operation updates the stored coordinate arrays and dimension metadata directly; it does not first materialize a full matrix.

## Parameters / inputs

- A -an RCV sparse matrix

## Outputs

- A -an RCV sparse matrix

## Implementation structure

- Requires A to be an RCV object.
- Swaps row and column coordinate arrays and swaps numRows with numCols.
- Replaces the value array with its complex conjugate.
