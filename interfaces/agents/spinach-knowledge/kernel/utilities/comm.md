# kernel/utilities/comm.m

- Signature: `C=comm(A,B)`

## Purpose

Computes the commutator of two square matrices, `C=A*B-B*A`.

## Physical / mathematical content

The commutator measures the difference between the two possible matrix-product orders, `A*B` and `B*A`.

## Parameters / inputs

- `A`, `B` — numeric square matrices

## Outputs

- `C` — square matrix equal to `A*B-B*A`

## Implementation structure

The function checks that both inputs are numeric square matrices, then evaluates `A*B-B*A`.
