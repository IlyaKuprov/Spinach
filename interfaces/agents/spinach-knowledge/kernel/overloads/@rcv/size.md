# kernel/overloads/@rcv/size.m

- Signature: `[s,ncols]=size(A,dim)`

## Purpose

Returns the dimensions of an RCV sparse matrix in MATLAB-style one-output, dimension-query, or two-output form.

## Mathematical content

With one output and no dimension index, `size` returns `[numRows numCols]`. With a dimension index, it returns the corresponding row or column count for dimensions 1 and 2, and `1` for any other positive integer dimension. With two outputs and no dimension index, it returns the row and column counts separately.

## Numerical / algorithmic content

The function requires an RCV input and, when supplied, a positive integer scalar dimension index. The two-output form cannot be combined with a dimension query.

## Parameters / inputs

- A -RCV sparse matrix
- dim -optional dimension index

## Outputs

- s -size vector, dimension length, or number
- of rows in the two-output form
- ncols -number of columns in the two-output form

## Implementation structure

- Validate the RCV input and optional dimension index.
- Reject a dimension query when two outputs are requested.
- Return `[numRows numCols]`, the selected dimension length (or `1` beyond dimension 2), or the separate row and column counts according to the requested form.
