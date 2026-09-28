# kernel/overloads/@ttclass/mtimes.m

- Signature: `c=mtimes(a,b)`

## Purpose

Implements scalar scaling, tensor-train multiplication, and multiplication of a tensor train by a full matrix.

## Numerical / algorithmic content

A scalar on either side scales a tensor train's coefficient and tolerance. For a tensor train multiplied by a full matrix, the matrix dimensions must match the product of the train's column-mode sizes; the code contracts the train against the matrix and returns a dense result. For two tensor trains, core counts and contracted mode sizes must agree. Their cores are contracted over the shared modes, coefficients are multiplied, zero-coefficient components are removed, and the result is compressed with `shrink`. Unsupported operand combinations raise an error.

## Parameters / inputs

- a -a scalar or a tensor train
- b -a scalar, a tensor train, or a full matrix

## Outputs

- c -a tensor train for scalar or tensor-train multiplication; a dense matrix when a tensor train is multiplied by a full matrix
