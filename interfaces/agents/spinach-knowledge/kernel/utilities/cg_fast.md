# kernel/utilities/cg_fast.m

- Signature: `cg=cg_fast(L,M,L1,M1,L2,M2)`

## Purpose

Computes the Clebsch-Gordan coefficient of `Y(L,M)` in the product expansion of `Y(L1,M1)` and `Y(L2,M2)`. Equivalently, it gives the expansion coefficient of the angular-momentum or spin state `|L,M>` in the product basis `|L1,M1>|L2,M2>`.

## Physical / mathematical content

Only combinations allowed by the selection rules for spherical harmonics or spin states are admissible. The function returns zero for inadmissible indices.

## Numerical / algorithmic content

Staged zero tests reject inadmissible combinations before a log-factorial summation evaluates the coefficient. The double-precision result is accurate to about `1e-3` up to about `L=20`; for higher ranks, `clebsch_gordan.m` provides a slower machine-precision implementation.

## Parameters / inputs

- `L, M, L1, M1, L2, M2` — integer or half-integer indices of the angular-momentum or spin states

## Outputs

- `cg` — double-precision floating-point Clebsch-Gordan coefficient
