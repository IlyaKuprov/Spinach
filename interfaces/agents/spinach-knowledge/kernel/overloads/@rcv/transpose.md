# kernel/overloads/@rcv/transpose.m

- Signature: `A=transpose(A)`

## Purpose

Return the transpose of an RCV sparse matrix.

## Physical / mathematical content

Transposition exchanges the row and column indices and swaps the matrix dimensions.

## Parameters / inputs

- `A` - RCV sparse matrix.

## Outputs

- `A` - transposed RCV sparse matrix.

## Implementation structure

After checking that the input is an `rcv` object, the function swaps its `row` and `col` arrays and its `numRows` and `numCols` metadata. Stored values are unchanged.
