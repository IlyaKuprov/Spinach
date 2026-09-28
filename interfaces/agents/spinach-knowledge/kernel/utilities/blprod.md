# kernel/utilities/blprod.m

- Signature: `[X1_AB,X2_AB]=blprod(A,B)`

## Purpose

Extends Blicharski's tensor invariants to scalar products of different spin interaction tensors using polarisation identities.

## Physical / mathematical content

The function returns first- and second-rank cross-correlation amplitudes for the supplied interaction tensors. The polarization identities obtain these amplitudes from the corresponding invariants of `A+B` and `A-B`; isotropic components of `A` and `B` do not affect the result.

## Numerical / algorithmic content

The inputs are checked as real numeric 3x3 matrices. The function evaluates `blinv` for `A-B` and `A+B`, then divides the difference between the corresponding first-rank invariants by four to obtain `X1_AB`, and does the same with the second-rank invariants for `X2_AB`.

## Parameters / inputs

- `A` — a real 3x3 matrix
- `B` — a real 3x3 matrix

## Outputs

- `X1_AB` — cross-correlation amplitude, first rank
- `X2_AB` — cross-correlation amplitude, second rank
- This function is not sensitive to the isotropic components of `A` and `B` tensors.
