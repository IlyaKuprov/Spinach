# kernel/overloads/@rcv/spy.m

- Signature: `spy(A)`

## Purpose

Display the nonzero pattern of an RCV sparse matrix using MATLAB's `spy`.

## Physical / mathematical content

The plotted pattern corresponds to the row and column indices stored in the RCV matrix.

## Parameters / inputs

- `A` - RCV sparse matrix.

## Outputs

A MATLAB sparsity plot; the function does not return a matrix.

## Implementation structure

The function checks that `A` is an `rcv` object, converts it to MATLAB sparse form with `sparse(A)`, and passes that matrix to MATLAB's `spy`.
