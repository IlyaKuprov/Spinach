# kernel/operators/irr_sph_ten.m

- Signature: `T=irr_sph_ten(mult,k)`

## Purpose

Constructs single-spin irreducible spherical tensor operators `T(k,m)`.

## Physical / mathematical content

The operators obey the commutation relation `[Lz,T(k,m)]=m*T(k,m)`. For a specified rank `k`, there are `2*k+1` components, returned in decreasing order of projection `m`. With only `mult`, the routine returns every rank from zero through `mult-1`, with decreasing projection within each rank.

The source notes that operator normalisation is not appropriate for its spin-dynamics convention: use identical commutation relations, rather than equal matrix norms, to make the formalism independent of the spin quantum number.

## Numerical / algorithmic content

For rank zero, the component is the identity matrix. For higher ranks, the highest-projection component is formed from the raising operator as `(-1)^k * 2^(-k/2) * L.p^k`; the remaining components are generated sequentially by the lowering-operator commutator using Racah's rule. The one-argument form concatenates these rank-specific cell arrays in increasing rank order.

## Parameters / inputs

- mult - spin multiplicity; a positive integer.
- k - irreducible spherical tensor rank (optional); an integer from zero to `mult-1`.

## Outputs

- T - with `mult` and `k`, a cell array of rank-`k` tensors in decreasing projection order; with `mult` alone, a cell array containing all ranks in increasing rank order and decreasing projection within each rank.

## Implementation structure

1. Validate `mult` and, when supplied, `k`.
2. For a one-argument call, recursively collect the tensors for all ranks from zero to `mult-1`.
3. For a two-argument call, return the identity for rank zero or generate the rank-`k` components from the highest projection by sequential lowering.
