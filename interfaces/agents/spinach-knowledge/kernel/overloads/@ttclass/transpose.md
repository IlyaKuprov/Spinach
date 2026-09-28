# kernel/overloads/@ttclass/transpose.m

- Signature: `ttrain=transpose(ttrain)`

## Purpose

Transpose the matrix represented by a tensor train without complex conjugation.

## Parameters / inputs

- `ttrain` — tensor-train representation of a matrix.

## Outputs

- `ttrain` — tensor train representing the transpose of the input matrix.

## Implementation

For each core in each train, the function permutes the core dimensions with `[1 3 2 4]`, swapping the row and column physical dimensions while leaving the bond dimensions unchanged. This is a non-conjugating transpose.

## Source

D. Savostyanov and I. Kuprov, [`ttclass/transpose.m`](https://spindynamics.org/wiki/index.php?title=ttclass/transpose.m).
