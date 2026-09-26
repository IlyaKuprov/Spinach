# kernel/overloads/@ttclass/plus.m

- Signature: `a=plus(a,b)`

## Purpose

Tensor train addition operation. Does not perform the actual addition, but instead concatenates the operands until such time as recompression becomes absolutely necessary. Syntax: c=plus(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- a -a tensor train object
- b -a tensor train object

## Outputs

- c -a tensor train object

## Implementation structure

- Tensor train addition operation. Does not perform the actual addition,
- but instead concatenates the operands until such time as recompression
- becomes absolutely necessary. Syntax:
- c=plus(a,b)
- a -a tensor train object
- b -a tensor train object
- c -a tensor train object
- Validate the input
- Write the sum object
- Filter out zero coeff
- Twinkle, twinkle, little star.
- I don't wonder what you are.
