# kernel/overloads/@ttclass/revert.m

- Signature: `tt=revert(tt)`

## Purpose

Reverse the order of the tensor-train cores and reverse each core's bond-index direction by swapping its outer bond dimensions.

## Parameters / inputs

- `tt` — tensor train operator.

## Outputs

- `tt` — tensor train operator with reversed core order and bond direction.

## Implementation

The function reads the number of cores and trains from `tt.cores`. For each train, it permutes every core with dimensions `[4,2,3,1]`, then reverses the core sequence. This performs the bit-revert permutation described by the source comments.

## Source

D. Savostyanov, [`ttclass/revert.m`](https://spindynamics.org/wiki/index.php?title=ttclass/revert.m).
