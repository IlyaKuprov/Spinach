# kernel/overloads/@rcv/mtimes.m

- Signature: `C=mtimes(A,B)`

## Purpose

Multiplication for RCV sparse matrices. Syntax: C=mtimes(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- C -product A*B as a Matlab sparse matrix if
- both operands are RCV or Matlab matrices

## Implementation structure

- Multiplication for RCV sparse matrices. Syntax:
- C=mtimes(A,B)
- A -left operand
- B -right operand
- C -product A*B as a Matlab sparse matrix if
- both operands are RCV or Matlab matrices
- Check consistency
- RCV sparse by a scalar
- Scalar by RCV sparse
- RCV sparse by RCV sparse
- Check dimension consistency
- Result is Matlab sparse
