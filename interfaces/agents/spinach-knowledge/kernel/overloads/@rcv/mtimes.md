# kernel/overloads/@rcv/mtimes.m

- Signature: `C=mtimes(A,B)`

## Purpose

Implements scalar scaling and matrix multiplication when at least one operand is an RCV sparse matrix.

## Mathematical content

A numeric scalar scales the stored values of an RCV operand. For matrix operands, the function computes the usual product `A*B` after requiring the inner dimensions to agree.

## Numerical / algorithmic content

RCV-by-RCV and mixed RCV/MATLAB-sparse matrix products are converted to MATLAB sparse matrices for multiplication and return a MATLAB sparse result. Scalar scaling returns an RCV object. The input check accepts an RCV operand paired with another RCV matrix, a MATLAB sparse matrix, or a numeric scalar.

## Parameters / inputs

- A -left operand
- B -right operand

## Outputs

- C -product A*B as a Matlab sparse matrix if
- both operands are RCV or Matlab matrices

## Implementation structure

- Check that at least one operand is RCV and that the other operand is an RCV matrix, MATLAB sparse matrix, or numeric scalar.
- For RCV-by-scalar or scalar-by-RCV, scale the RCV `val` array and return the RCV operand.
- For RCV-by-RCV, RCV-by-MATLAB-sparse, or MATLAB-sparse-by-RCV, check inner dimensions, convert the RCV operand(s) to MATLAB sparse form, and multiply.
