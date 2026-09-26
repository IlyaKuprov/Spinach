# kernel/overloads/@rcv/full.m

- Signature: `A=full(A)`

## Purpose

Converts an RCV sparse matrix into a full matrix. Syntax: A=full(A)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A -an RCV sparse matrix

## Outputs

- A -a full Matlab matrix

## Implementation structure

- Converts an RCV sparse matrix into a full matrix. Syntax:
- A=full(A)
- A -an RCV sparse matrix
- A -a full Matlab matrix
- Check consistency
- Delegate to Matlab
- Consistency enforcement
- Whenever you find yourself on the side of the
- majority, it is time to pause and reflect.
- Mark Twain
