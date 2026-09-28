# kernel/utilities/killcross.m

- Signature: `M=killcross(M,f1idx,f2idx)`

## Purpose

Set selected columns and rows of a matrix to zero.

## Physical / mathematical content

A matrix utility; it makes no assumptions about the matrix's physical interpretation.

## Numerical / algorithmic content

The entries in columns `f1idx` and rows `f2idx` are assigned zero. Column and row selections are applied directly to the input matrix.

## Parameters / inputs

- `M` - numeric matrix.
- `f1idx` - unique positive integer indices of columns to zero.
- `f2idx` - unique positive integer indices of rows to zero.

## Outputs

- `M` - the matrix with the selected rows and columns set to zero.

## Implementation structure

After checking that the matrix is two-dimensional and that each index list contains unique in-range positive integers, the routine zeros `M(f2idx,:)` and `M(:,f1idx)`.
