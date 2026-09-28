# kernel/overloads/@ttclass/diag.m

- Signature: `tt=diag(tt)`

## Purpose

Applies the `diag` operation to a tensor-train matrix or vector.

## Input

- `tt` — tensor-train representation of a matrix or vector.

## Output

- `tt` — for a square-matrix input, its diagonal as a tensor-train vector; for a vector input, the corresponding diagonal matrix as a tensor train.

## Behavior

A vector is recognized when all row mode sizes or all column mode sizes are one. The function constructs a diagonal-matrix core for each vector core. For a matrix input, it requires each mode to be square and replaces each core's matrix-mode data with its diagonal. Inputs that are neither a vector nor a square matrix raise an error.
