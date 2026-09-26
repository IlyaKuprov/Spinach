# kernel/overloads/@rcv/ctranspose.m

- Signature: `A=ctranspose(A)`

## Purpose

Returns the conjugate transpose of an RCV sparse matrix. Syntax: A=ctranspose(A)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A -an RCV sparse matrix

## Outputs

- A -an RCV sparse matrix

## Implementation structure

- Returns the conjugate transpose of an RCV sparse matrix. Syntax:
- A=ctranspose(A)
- A -an RCV sparse matrix
- Check consistency
- Efficiently swap rows and columns
- Update row and column dimension information
- Conjugate values
- Consistency enforcement
- At the Tower of London we remembered in
- our prayers the elephant kept there by
- James I, which, poor creature, was never
- given anything to drink but wine.
