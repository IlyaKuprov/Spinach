# kernel/overloads/@rcv/gpuArray.m

- Signature: `obj=gpuArray(obj)`

## Purpose

Transfers an RCV sparse matrix to the GPU. Syntax: obj=gpuArray(obj)

## Physical / mathematical content

- RCV sparse-matrix storage utilities. The focus is data structure design for sparse linear algebra and low-overhead composition of large matrices.

## Numerical / algorithmic content

- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Parameters / inputs

- obj -an RCV sparse matrix

## Outputs

- obj -the same matrix with data stored on GPU

## Implementation structure

- Transfers an RCV sparse matrix to the GPU. Syntax:
- obj=gpuArray(obj)
- obj -an RCV sparse matrix
- obj -the same matrix with data stored on GPU
- Check consistency
- Upload to GPU
- Consistency enforcement
- Then it got worse. The book is very, very good. If
- someone's going to beat you to the punch with a great
- book idea, the least they can do is write something
- crap. Not Andrew. Which shouldn't really come as a
- surprise, since the little bastard is prodigiously
