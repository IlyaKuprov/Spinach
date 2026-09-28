# kernel/utilities/clebsch_gordan.m

- Signature: `cg=clebsch_gordan(L,M,L1,M1,L2,M2)`

## Purpose

Computes the Clebsch-Gordan coefficient of `Y(L,M)` in the product expansion of `Y(L1,M1)` and `Y(L2,M2)`. Equivalently, it gives the expansion coefficient of the angular-momentum or spin state `|L,M>` in the product basis `|L1,M1>|L2,M2>`.

## Physical / mathematical content

Only combinations allowed by the selection rules for spherical harmonics or spin states are admissible; the function returns zero for inadmissible indices.

## Numerical / algorithmic content

The implementation checks admissibility before evaluating the coefficient in double precision. Its header reports machine-precision answers up to about `L=1e4`; `cg_fast.m` is a faster alternative for low ranks.

## Parameters / inputs

- `L, M, L1, M1, L2, M2` — integer or half-integer indices of the angular-momentum or spin states

## Outputs

- `cg` — double-precision floating-point Clebsch-Gordan coefficient
