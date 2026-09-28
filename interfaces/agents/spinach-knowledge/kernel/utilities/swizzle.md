# kernel/utilities/swizzle.m

- Signature: `tuples=swizzle(index_arrays)`

## Purpose

Flattens nested index lists into a matrix of tuples in random order, useful for distributing nested-loop iterations across parallel workers.

## Parameters / inputs

- `index_arrays` — a cell array of row vectors containing positive integers.

## Outputs

- `tuples` — a matrix with one tuple per row, in random order. Each column corresponds to an input index array.

## Algorithm

The function validates `index_arrays`, constructs all combinations of its entries using Kronecker products, then randomly permutes the resulting rows.

## Source

<https://spindynamics.org/wiki/index.php?title=swizzle.m>