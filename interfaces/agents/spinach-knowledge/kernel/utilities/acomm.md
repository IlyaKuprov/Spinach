# kernel/utilities/acomm.m

- Signature: `C=acomm(A,B)`

## Purpose

Computes the matrix anticommutator `C = A*B + B*A`.

## Parameters / inputs

- `A`, `B` - numeric square matrices of the same dimensions.

## Output

- `C` - square matrix, computed as `A*B + B*A`.

## Behavior

The function validates that both inputs are numeric square matrices with equal dimensions, then computes the anticommutator.
