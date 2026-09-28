# kernel/conventions/transforms/mat2sphten.m

- Signature: `[rank0,rank1,rank2]=mat2sphten(M)`

## Purpose

Converts a real 3x3 interaction matrix into coefficients of irreducible spherical-tensor operators: one rank-0, three rank-1, and five rank-2 components. The conventions are matched to Equation (22) of the paper by Len Mueller.

## Physical / mathematical content

The returned coefficients multiply the corresponding operators from `irr_sph_ten.m`. Their order is rank 0: (0,0); rank 1: (1,1), (1,0), (1,-1); rank 2: (2,2), (2,1), (2,0), (2,-1), (2,-2).

## Numerical / algorithmic content

The rank-0 coefficient is `trace(M)/3`; the rank-1 and rank-2 coefficients are the explicit linear combinations of matrix elements in the source. An empty input is replaced with a zero 3x3 matrix before validation.

## Syntax

```matlab
[rank0,rank1,rank2]=mat2sphten(M)
```

## Parameters / inputs

- `M` — a real numeric 3x3 interaction matrix. Empty input is accepted and treated as `zeros(3)`.

## Outputs

- `rank0` — scalar coefficient of T(0,0).
- `rank1` — 3x1 column vector, ordered as T(1,1), T(1,0), T(1,-1).
- `rank2` — 5x1 column vector, ordered as T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).

## Implementation structure

After replacing empty input, the function checks that M is real, numeric, and 3x3, then computes the nine coefficients.
