# tests/kernel/test_operator_basis_suite.m

- Signature: `result=test_operator_basis_suite()`

## Purpose

Regression test for operator-basis construction and expansion helpers. Returns `result` with explanatory test messages.

## Checks

- For a spin of multiplicity 3, rank-2 irreducible spherical tensors satisfy `[Lz,T(k,m)]=m*T(k,m)` for projections 2 through -2; the complete single-spin tensor basis has `mult^2` operators. The rank-one, zero-projection Stevens operator equals `Lz`.
- In truncated bosonic bases, `weyl` satisfies `c*a=n`, `[n,c]=c`, and `[n,a]=-a`. `boson_mono` begins with the identity and creation operator, while distinct `boson_ortho` operators have zero Hilbert-Schmidt overlap.
- `sin_tran(3)` produces nine matrices, each with one nonzero element at its serpentine-indexed position. For multiplicity 4, `centrans` places its `z` populations and raising element only on the central two levels.
- Reconstructions from `oper2ist`, `ct2ist`, `enlev2ist`, and `bos2ist` match their source matrices or operators; reconstructions from `oper2bm` and `enlev2bm` match the requested projector.
- For two spin-half particles, `twospinist` rank 2, projection 0 matches its Cartesian-product expression. `mprealloc` returns sizes `[4 4]` in `zeeman-hilb` and `[16 16]` in `zeeman-liouv`.