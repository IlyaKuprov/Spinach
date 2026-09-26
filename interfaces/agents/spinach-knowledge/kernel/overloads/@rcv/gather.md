# kernel/overloads/@rcv/gather.m

- Signature: `A=gather(A)`

## Purpose

Gathers an RCV sparse matrix from GPU. Syntax: A=gather(A)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

## Parameters / inputs

- A -an RCV sparse matrix

## Outputs

- A -the same matrix with data stored on the CPU

## Implementation structure

- Gathers an RCV sparse matrix from GPU. Syntax:
- A=gather(A)
- A -an RCV sparse matrix
- A -the same matrix with data stored on the CPU
- Check consistency
- Gather to CPU
- Consistency enforcement
- Aerie, I've noticed the unfortunate fact that you live
- by one of the great lessons of history that nothing is
- often a good thing to do and a clever thing to say.
- Edwin Odesseiron, in Baldur's Gate 2
