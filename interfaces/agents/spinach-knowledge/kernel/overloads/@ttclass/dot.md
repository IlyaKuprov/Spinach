# kernel/overloads/@ttclass/dot.m

- Signature: `c=dot(a,b)`

## Purpose

Dot product of TT representations of matrices. Syntax: c=dot(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- a,b -tensor train objects representing numerical
- arrays of consistent dimensions and having
- the same internal topology

## Outputs

- c -inner product of a and b

## Implementation structure

- Dot product of TT representations of matrices. Syntax:
- c=dot(a,b)
- a,b -tensor train objects representing numerical
- arrays of consistent dimensions and having
- the same internal topology
- c -inner product of a and b
- Check consistency
- Compute the product
- Consistency enforcement
- Any product that needs a manual to work is broken.
- Elon Musk
