# kernel/overloads/@ttclass/mtimes.m

- Signature: `c=mtimes(a,b)`

## Purpose

Performs tensor train multiplication followed by a shrink. Syntax: c=mtimes(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- a -a scalar or a tensor train
- b -a scalar, a tensor train, or a full matrix

## Outputs

- c -a tensor train object

## Implementation structure

- Performs tensor train multiplication followed by a shrink. Syntax:
- c=mtimes(a,b)
- a -a scalar or a tensor train
- b -a scalar, a tensor train, or a full matrix
- c -a tensor train object
- Decide the type combination
- Multiply tensor train by a scalar from the right
- Read sizes and ranks of the operands
- Check consistency
- Preallocate result
- Loop over the buffers of the operands
- Set current vector as the right-hand side
