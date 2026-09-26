# kernel/overloads/@ttclass/shrink.m

- Signature: `ttrain=shrink(ttrain)`

## Purpose

Approximates a given tensor train with lower TT-ranks. Syntax: ttrain=shrink(ttrain)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- ttrain -a tensor train object

## Outputs

- ttrain -compressed tensor train object with
- right-to-left orthogonalisation

## Implementation structure

- Approximates a given tensor train with lower TT-ranks. Syntax:
- ttrain=shrink(ttrain)
- ttrain -a tensor train object
- ttrain -compressed tensor train object with
- right-to-left orthogonalisation
- Read train sizes
- Summation
- Left-to-right orthogonalisation
- Check the norm and escape if the object is zero
- Truncation
- Convert to a scalar if appropriate
- It's a tough life, being small and delicious.
