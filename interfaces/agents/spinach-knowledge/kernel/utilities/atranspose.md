# kernel/utilities/atranspose.m

## Purpose

`atranspose.m` performs an anti-diagonal transpose of a numeric array, returning the array reflected across its anti-diagonal.

## Behaviour

The function calls `grumble(M)` to enforce consistency: if `M` is not numeric, it errors with `'M must be a numeric array.'`. Otherwise, it computes the result as `transpose(rot90(M,2))` — rotating the array by 180 degrees and then applying the standard transpose — and returns the transformed array.

## Inputs and outputs

**Syntax:** `M=atranspose(M)`

- **Input:** `M` — a transposable (numeric) array.
- **Output:** `M` — the anti-diagonal transposed array.

## References

- Source: [kernel/utilities/atranspose.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/atranspose.m)
- Wiki: [atranspose.m — Spinach](https://spindynamics.org/wiki/index.php?title=atranspose.m)
