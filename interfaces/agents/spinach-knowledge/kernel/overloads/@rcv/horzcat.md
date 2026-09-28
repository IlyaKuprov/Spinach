# kernel/overloads/@rcv/horzcat.m

- Signature: `A=horzcat(A,B)`

## Purpose

Horizontally concatenates RCV sparse matrices in the order supplied, producing a matrix with the common row count and the sum of their column counts.

## Physical / mathematical content

The output places each input matrix in a consecutive column block; row indices are preserved and column indices are offset by the widths of preceding blocks.

## Numerical / algorithmic content

The variadic implementation accepts the input sequence as varargin. If any input is GPU-resident, all inputs are moved to GPU before their coordinate and value arrays are concatenated.

## Parameters / inputs

- A -left RCV sparse matrix
- B -right RCV sparse matrix

## Outputs

- A -RCV sparse matrix

## Implementation structure

- Requires every input to be an RCV object and all row counts to match.
- Moves all operands to GPU if at least one input is GPU-resident.
- Offsets each input's column indices by the cumulative column count.
- Concatenates the row, adjusted column, and value arrays; sets numCols to the sum of input column counts.
