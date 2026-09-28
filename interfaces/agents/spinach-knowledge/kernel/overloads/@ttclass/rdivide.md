# kernel/overloads/@ttclass/rdivide.m

- Signature: `a=rdivide(a,b)`

## Purpose

Divides a tensor train object by a scalar by dividing its coefficients by the scalar and scaling its tolerances by the scalar's absolute value.

## Physical / mathematical content

The operation is accepted only when the first argument is a `ttclass` object and the second is scalar; otherwise it raises an error. The core tensors are unchanged.

## Numerical / algorithmic content

The implementation performs `a.coeff=a.coeff/b` and `a.tolerance=a.tolerance/abs(b)`.

## Parameters / inputs

- a - a ttclass object
- b - a numeric scalar

## Outputs

- c - a ttclass object

## Implementation structure

- Divide the coefficients by `b` and tolerances by `abs(b)`.
- Raise an error unless the first input is a tensor train and the second input is scalar.
