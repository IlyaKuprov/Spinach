# kernel/overloads/@rcv/times.m

- Signature: `C=times(A,B)`

## Purpose

Multiplies an RCV sparse matrix by a numeric scalar, in either operand order. Syntax: C=times(A,B)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A,B -an RCV sparse matrix and a numeric
- scalar, in either order

## Outputs

- C -RCV sparse matrix

## Implementation structure

- Multiplies an RCV sparse matrix by a numeric scalar,
- in either operand order. Syntax:
- C=times(A,B)
- A,B -an RCV sparse matrix and a numeric
- scalar, in either order
- C -RCV sparse matrix
- Check consistency
- RCV sparse by a scalar
- Scalar by RCV sparse
- Consistency enforcement
- They say that the fish that gets away
- looks bigger than it really is.
