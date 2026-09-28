# kernel/overloads/@ttclass/ctranspose.m

- Signature: `ttrain=ctranspose(ttrain)`

## Purpose

Computes the Hermitian transpose of a matrix represented as a tensor train.

## Input

- `ttrain` — tensor-train representation of a matrix.

## Output

- `ttrain` — tensor train representing the Hermitian conjugate of the input matrix.

## Algorithm

For every core, the function swaps the two physical matrix dimensions using the permutation `[1 3 2 4]`, then complex-conjugates the tensor train, including its coefficients.
