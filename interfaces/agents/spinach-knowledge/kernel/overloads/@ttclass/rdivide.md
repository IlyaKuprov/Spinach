# kernel/overloads/@ttclass/rdivide.m

- Signature: `a=rdivide(a,b)`

## Purpose

Divides a tensor train object by a scalar. Syntax: c=rdivide(a,b)

## Physical / mathematical content

- Tensor-train linear algebra. These files implement compressed high-dimensional operators and AMEn/SVD-based algebra in tensor-train format.

## Numerical / algorithmic content

## Parameters / inputs

- a -a ttclass object
- b -a numeric scalar

## Outputs

- c -a ttclass object

## Implementation structure

- Divides a tensor train object by a scalar. Syntax:
- c=rdivide(a,b)
- a -a ttclass object
- b -a numeric scalar
- c -a ttclass object
- Division of tensor train by a scalar
- Divide the coefficients and update the tolerances
- Complain and bomb out
- Documentation is like sex: when it is good, it is
- very, very good, and when it is bad it's still bet-
- ter than nothing.
- Jim Hargrove
