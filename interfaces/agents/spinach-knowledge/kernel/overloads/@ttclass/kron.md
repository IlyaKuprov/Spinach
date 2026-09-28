# kernel/overloads/@ttclass/kron.m

- Signature: `c=kron(a,b)`

## Purpose

Constructs the Kronecker product of two tensor-train matrix representations, core by core.

## Numerical / algorithmic content

Both operands are shrunk before their sizes and ranks are read. The operation requires the same number of cores, forms each output core from the Kronecker product of reshaped operand cores, multiplies the coefficients, and sets the output tolerance from the operand coefficients and tolerances.

## Parameters / inputs

- a,b -tensor train objects

## Outputs

- c -a tensor train object

## Ordering caveat

The result is not ordered like the flat matrix Kronecker product: it differs by a row and column permutation. Its element order is consistent with the output of `ttclass/vec`.
