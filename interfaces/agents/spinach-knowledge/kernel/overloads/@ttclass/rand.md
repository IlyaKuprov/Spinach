# kernel/overloads/@ttclass/rand.m

- Signature: `tt=rand(tt,ttrank)`

## Purpose

Generates a tensor train representation of a matrix with random tensor train cores, same physical index topology as the tensor train supplied, and user-specified bond ranks. Syntax: tt=rand(tt,ttrank)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -a tensor train object
- ttrank -bond rank, a positive integer

## Outputs

- tt -a tensor train object

## Implementation structure

- Generates a tensor train representation of a matrix with random
- tensor train cores, same physical index topology as the tensor
- train supplied, and user-specified bond ranks. Syntax:
- tt=rand(tt,ttrank)
- tt -a tensor train object
- ttrank -bond rank, a positive integer
- Check consistency
- Read tensor train sizes
- Reallocate cores
- Fill the cores with random elements
- Unit coefficient and zero tolerance
- Consistency enforcement
