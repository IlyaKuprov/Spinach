# kernel/overloads/@rcv/kron.m

- Signature: `C=kron(A,B)`

## Purpose

Forms the Kronecker product of two RCV sparse matrices.

## Mathematical content

For each stored entry pair `A(i,j)` and `B(k,l)`, the result stores their product at row `(i-1)*B.numRows+k` and column `(j-1)*B.numCols+l`. Its dimensions are `A.numRows*B.numRows` by `A.numCols*B.numCols`.

## Numerical / algorithmic content

The function checks that both inputs are RCV matrices. If either input is marked as GPU-resident, it converts both to GPU arrays. It forms all pairs of stored-entry indices with `ndgrid`, computes the corresponding output indices and value products, and constructs the result with `rcv`; the output GPU flag is the logical OR of the input flags.

## Parameters / inputs

- A -left RCV sparse matrix
- B -right RCV sparse matrix

## Outputs

- C -RCV sparse matrix

## Implementation structure

- Check that both inputs are RCV matrices.
- Compute the product dimensions and, when needed, move both inputs to GPU arrays.
- Form all stored-entry pairs, map their indices into the Kronecker-product matrix, and multiply their values.
- Assemble the output RCV matrix and set its GPU flag from the inputs.
