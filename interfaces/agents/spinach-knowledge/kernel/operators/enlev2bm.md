# kernel/operators/enlev2bm.m

- Signature: `[states,coeffs]=enlev2bm(nlevels,lvl_num)`
- DIRECT source: [kernel/operators/enlev2bm.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/enlev2bm.m)
- Wiki: [enlev2bm.m](https://spindynamics.org/wiki/index.php?title=enlev2bm.m)

## Definition

The function represents one energy-level projector of a truncated bosonic mode in Spinach's bosonic-monomial basis. It creates an `nlevels`-by-`nlevels` diagonal matrix `P`, sets only `P(lvl_num,lvl_num)=1`, and passes `P` to [oper2bm](oper2bm.md). The level number is one-based: level 1 is the empty-mode state, and increasing indices count upward from it.

This is a basis expansion of a projector, not a propagator. The function does not exponentiate an operator or apply time evolution.

## Basis and outputs

`states` contains the Spinach BM-basis indices returned by `oper2bm(P)`; use [lin2kq](../indexing/lin2kq.md) for the K,Q bosonic-monomial indexing. `coeffs` is returned directly from the same conversion call. `enlev2bm` adds no scale factor or other coefficient normalisation of its own; the values are those calculated by `oper2bm` for this diagonal projector.

## Inputs and checks

`nlevels` is the mode's number of levels, intended as a positive integer. The source checks that it is numeric, scalar, real, and at least 1; it does not explicitly test integrality or finiteness. `lvl_num` must be numeric, scalar, real, and between 1 and `nlevels`; the source uses it directly as a MATLAB matrix index, so it must also be a valid integer index.
