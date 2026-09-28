# kernel/utilities/atranspose.m

- Signature: `M=atranspose(M)`

## Purpose

Reflect a numeric matrix across its anti-diagonal. The implementation rotates it by 180 degrees with `rot90(M,2)` and then applies MATLAB's non-conjugating transpose.

## Parameters / inputs

- `M` - numeric array accepted by MATLAB's `rot90` and `transpose` operations.

## Output

- `M` - the anti-diagonally transposed array.
