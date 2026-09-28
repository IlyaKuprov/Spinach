# kernel/utilities/binpack.m

- Signature: `bins=binpack(box_sizes,bin_size)`

## Purpose

Greedily group box indices into bins without reordering the input. This is a simple, non-optimal one-dimensional packing algorithm.

## Parameters / inputs

- `box_sizes` - row vector of positive integers.
- `bin_size` - positive integer capacity.

## Outputs

- `bins` - cell array of index vectors, one per bin. Boxes larger than `bin_size` are returned in singleton bins, so those bins can exceed the specified capacity.

## Algorithm

The function first places each oversized box in its own bin. It then repeatedly places the longest remaining input prefix whose cumulative size does not exceed `bin_size` into the next bin. Consequently, it preserves input order and does not search for a globally optimal assignment.
