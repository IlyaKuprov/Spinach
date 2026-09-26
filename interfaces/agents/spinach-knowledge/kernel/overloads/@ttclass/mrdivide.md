# kernel/overloads/@ttclass/mrdivide.m

- Signature: `a=mrdivide(a,b)`

## Purpose

Divides a tensor train object by a scalar. Syntax: c=mrdivide(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- a -tensor train object
- b -a scalar

## Outputs

- c -tensor train object

## Implementation structure

- Divides a tensor train object by a scalar. Syntax:
- c=mrdivide(a,b)
- a -tensor train object
- b -a scalar
- c -tensor train object
- Division of tensor train by a scalar
- Divide the coefficients and update the tolerances
- Complain and bomb out
- It is dangerous to be right in matters on which the established
- authorities are wrong.
- Voltaire
