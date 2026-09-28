# kernel/overloads/@rcv/full.m

- Signature: `A=full(A)`

## Purpose

Converts an RCV sparse matrix to a full MATLAB matrix.

## Physical / mathematical content

This changes the matrix storage representation, not its represented entries: the output is the dense matrix corresponding to the input's stored coordinates and values.

## Numerical / algorithmic content

After checking that A is an RCV object, the implementation constructs a MATLAB sparse matrix and delegates conversion to MATLAB's `full` function.

## Parameters / inputs

- A -an RCV sparse matrix

## Outputs

- A -a full Matlab matrix

## Implementation structure

- Requires A to be an RCV object.
- Converts A to MATLAB sparse storage, then calls full on that sparse matrix.
