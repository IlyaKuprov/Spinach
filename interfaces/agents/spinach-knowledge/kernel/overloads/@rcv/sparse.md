# kernel/overloads/@rcv/sparse.m

- Signature: `A=sparse(A)`

## Purpose

Converts an RCV sparse matrix into a Matlab sparse matrix. Syntax: A=sparse(A)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A -RCV sparse matrix

## Outputs

- A -Matlab sparse matrix

## Implementation structure

- Converts an RCV sparse matrix into a Matlab sparse
- matrix. Syntax:
- A=sparse(A)
- A -RCV sparse matrix
- A -Matlab sparse matrix
- Check consistency
- Check if empty
- Empty matrix of a specified size
- Call Matlab's sparse matrix constructor
- Consistency enforcement
- Working 16 hours a day, 7 days a week, 52 weeks
- in a year, and people still calling me lucky.
