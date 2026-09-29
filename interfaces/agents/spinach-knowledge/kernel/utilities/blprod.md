# kernel/utilities/blprod.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/blprod.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/blprod.m)

## Purpose

Extends Blicharski's tensor invariants into scalar products of different spin interaction tensors using polarisation identities. The function computes first-rank and second-rank cross-correlation amplitudes between two real 3x3 interaction tensors.

## Behaviour

- Validates both inputs through an internal consistency check (`grumble`), which errors with `'A must be a real 3x3 matrix.'` or `'B must be a real 3x3 matrix.'` if an input is not numeric, not real, not a matrix, or not of size 3x3.
- Computes Blicharski invariants of the difference `A-B` and the sum `A+B` via `blinv`, returning `[LsqAmB,DsqAmB]` and `[LsqApB,DsqApB]` respectively.
- Applies the polarisation identity for the first rank: `X1_AB=(LsqApB-LsqAmB)/4`.
- Applies the polarisation identity for the second rank: `X2_AB=(DsqApB-DsqAmB)/4`.
- The function is not sensitive to the isotropic components of the `A` and `B` tensors.

## Inputs and outputs

**Inputs**

- `A` — a real 3x3 matrix.
- `B` — a real 3x3 matrix.

**Outputs**

- `X1_AB` — cross-correlation amplitude, first rank.
- `X2_AB` — cross-correlation amplitude, second rank.

## References

- Spinach Wiki: [blprod.m](https://spindynamics.org/wiki/index.php?title=blprod.m)
