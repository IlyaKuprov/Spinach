# kernel/overloads/@rcv/sparse.m

- Signature: `A=sparse(A)`

## Purpose

Converts an RCV sparse-matrix object to a MATLAB sparse matrix of the same dimensions.

## Mathematical content

The conversion preserves the matrix represented by the stored row indices, column indices, and values. An RCV object with no stored entries is converted to an empty sparse matrix with its specified dimensions.

## Numerical / algorithmic content

The function requires an RCV input. For a nonempty object, it passes the stored entry arrays and dimensions to MATLAB’s `sparse` constructor; for an empty object, it uses `spalloc` with zero allocated entries.

## Parameters / inputs

- A -RCV sparse matrix

## Outputs

- A -Matlab sparse matrix

## Implementation structure

- Check that the input is RCV.
- If it has no stored entries, create a zero-allocation sparse matrix with the recorded dimensions.
- Otherwise, construct the MATLAB sparse matrix from the stored row, column, and value arrays and the recorded dimensions.
