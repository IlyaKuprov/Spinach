# kernel/overloads/@ttclass/pack.m

- Signature: `ttout=pack(tt)`

## Purpose

This subroutine packs all trains from the addition buffer into a single tensor train, but does not perform the recom- pression. Normally you should not call it directly, use ttclass/shrink.m instead. Syntax: ttout=pack(tt)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- tt -tensor train object with unprocessed additions

## Outputs

- ttout -tensor train with additions buffer absorbed
- into the cores of the tensor, but not re-
- compressed

## Implementation structure

- This subroutine packs all trains from the addition buffer
- into a single tensor train, but does not perform the recom-
- pression. Normally you should not call it directly, use
- ttclass/shrink.m instead. Syntax:
- ttout=pack(tt)
- tt - tensor train object with unprocessed additions
- ttout -tensor train with additions buffer absorbed
- into the cores of the tensor, but not re-
- compressed
- Read tensor ranks and dimensions
- Fast return if possible
- Total rank of all summands
