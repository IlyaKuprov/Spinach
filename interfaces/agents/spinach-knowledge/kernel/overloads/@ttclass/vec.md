# kernel/overloads/@ttclass/vec.m

- Signature: `A=vec(A)`

## Purpose

Stretches arrays into vectors -useful for situations when the stand- ard (:) syntax is not available. Syntax: A=vec(A)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- A -numeric or ttclass array

## Outputs

- A -numeric or ttclass array
- WARNING: for tensor trains this operation proceeds by stretching eve-
- ry core of the tensor train. the result is not the same as
- column-wise matrix stretching (it is an element permutation
- away from it), but the resulting order of elements is consi-
- stent with tensor train Kronecker product operation output.

## Implementation structure

- Stretches arrays into vectors -useful for situations when the stand-
- ard (:) syntax is not available. Syntax:
- A=vec(A)
- A -numeric or ttclass array
- WARNING: for tensor trains this operation proceeds by stretching eve-
- ry core of the tensor train. the result is not the same as
- column-wise matrix stretching (it is an element permutation
- away from it), but the resulting order of elements is consi-
- stent with tensor train Kronecker product operation output.
- Decide how to proceed
- Read tensor train sizes and ranks
- Reshape the cores
