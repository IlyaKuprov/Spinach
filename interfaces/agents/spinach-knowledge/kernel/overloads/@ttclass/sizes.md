# kernel/overloads/@ttclass/sizes.m

- Signature: `modesizes=sizes(tt)`

## Purpose

Return the physical row and column dimensions for each core of a tensor train.

## Parameters / inputs

- `tt` — tensor train object.

## Outputs

- `modesizes` — an `ncores`-by-2 array; each row contains the second and third dimensions of the corresponding core in the first train.

## Implementation

The function obtains the number of cores from `tt.cores`, allocates the output array, and fills each row from `size(tt.cores{k,1},2)` and `size(tt.cores{k,1},3)`.

## Source

D. Savostyanov and I. Kuprov, [`ttclass/sizes.m`](https://spindynamics.org/wiki/index.php?title=ttclass/sizes.m).
