# kernel/overloads/@ttclass/mrdivide.m

- Signature: `a=mrdivide(a,b)`

## Purpose

Divides a tensor-train object by a scalar.

## Numerical / algorithmic content

For a `ttclass` first argument and scalar second argument, the function divides the stored coefficients by the scalar and divides the tolerances by its absolute value. The core arrays are unchanged. Other input combinations raise an error.

## Parameters / inputs

- a -tensor train object
- b -a scalar

## Outputs

- a -the tensor train object after division
