# kernel/overloads/@ttclass/hdot.m

- Signature: `c=hdot(a,b)`

## Purpose

Hadamard dot product between two tensor train matrices. Syntax: c=hdot(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- a,b -tensor train objects representing numerical
- arrays of the same dimensions and having
- the same internal topology

## Outputs

- c -Hadamard product of a and b, a scalar

## Implementation structure

- Hadamard dot product between two tensor train matrices. Syntax:
- c=hdot(a,b)
- a,b -tensor train objects representing numerical
- arrays of the same dimensions and having
- the same internal topology
- c -Hadamard product of a and b, a scalar
- Check consistency
- Read topology and initialize the answer
- Loop over TT buffers
- Multiply coefficients
- Loop over TT cores and compute dot product
- Add to the total
