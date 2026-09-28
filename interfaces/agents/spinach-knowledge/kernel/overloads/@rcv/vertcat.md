# kernel/overloads/@rcv/vertcat.m

- Signature: `A=vertcat(A,B)`

## Purpose

Vertically concatenate RCV sparse matrices in top-to-bottom order.

## Physical / mathematical content

The inputs must have matching column counts. Concatenation offsets each input's row indices by the number of rows already appended.

## Parameters / inputs

- One or more RCV sparse matrices, supplied in top-to-bottom order; all must have the same number of columns.

## Outputs

- `A` - RCV sparse matrix containing the input rows in order.

## Implementation structure

The function checks the input types and column counts. If any input is on the GPU, it converts all inputs to GPU arrays. It then offsets and concatenates the row, column, and value arrays and sets the result's row count to the sum of the input row counts.
