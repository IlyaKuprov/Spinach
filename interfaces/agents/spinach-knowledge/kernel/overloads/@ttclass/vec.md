# kernel/overloads/@ttclass/vec.m

- Signature: `A=vec(A)`

## Purpose

Stretches arrays into vectors -useful for situations when the standard (:) syntax is not available. Syntax: A=vec(A)

## Physical / mathematical content

Reshapes numeric arrays or tensor-train cores as vectors.

## Numerical / algorithmic content

For a `ttclass` input, each core of each train is reshaped to a column-like core, preserving the order used by tensor-train Kronecker products. Other inputs are reshaped to a column vector with `reshape(A,[numel(A) 1])`.

## Parameters / inputs

- `A` — numeric array or `ttclass` array.

## Outputs

- `A` — numeric or `ttclass` array.

## Note

For tensor trains, stretching each core does not give the same element order as column-wise stretching the full matrix; it differs by an element permutation, while remaining consistent with tensor-train Kronecker-product output.
