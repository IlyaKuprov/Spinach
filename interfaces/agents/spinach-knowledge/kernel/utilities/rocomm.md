# kernel/utilities/rocomm.m

- Signature: `C=rocomm(A)`

## Purpose

Computes the right-nested commutator of the matrices in cell array `A`: `[[...[A{1},A{2}],A{3}],...]`, where `[P,Q]=P*Q-Q*P`.

## Physical / mathematical content

This is a matrix-algebra utility. It applies successive commutators and does not impose a spin-specific physical interpretation on the input matrices.

## Numerical / algorithmic content

Starting with `C=A{1}`, the function updates `C=C*A{n}-A{n}*C` for each remaining cell entry, so the accumulated result is commuted on the right with the next matrix.

## Parameters / inputs

- `A` - cell array of numeric square matrices.

## Outputs

- `C` - right-ordered nested commutator. For a one-element cell array, the result is that sole matrix.

## Implementation structure

The function checks that `A` is a cell array whose entries are numeric square matrices, then performs the nested product-and-subtraction recurrence in cell order.

## Reference

[Spin Dynamics Wiki: rocomm.m](https://spindynamics.org/wiki/index.php?title=rocomm.m)
