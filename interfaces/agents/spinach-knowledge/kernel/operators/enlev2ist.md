# kernel/operators/enlev2ist.m

- Signature: `[states,coeffs]=enlev2ist(mult,lvl_num,particle)`
- DIRECT source: [kernel/operators/enlev2ist.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/operators/enlev2ist.m)
- Wiki: [enlev2ist.m](https://spindynamics.org/wiki/index.php?title=enlev2ist.m)

## Definition

This function expands the projector onto one selected level of a spin or boson in Spinach's irreducible spherical-tensor (IST) basis. It first forms a `mult`-by-`mult` diagonal matrix `P` with one unit diagonal entry, then returns the results of [oper2ist](oper2ist.md). The meaning of `lvl_num` depends on `particle`:

- For `particle='S'`, spin levels are numbered from the bottom up, and the matrix index set to one is `mult-lvl_num+1`. Thus the level numbering is reversed relative to the diagonal index.
- For `particle='B'`, bosonic levels are numbered from the top down, and the set diagonal index is `lvl_num`.

These are indexing and basis-conversion operations; the function does not exponentiate the projector or produce a propagator.

## Basis, coefficients, and order

`states` are the IST-basis indices supplied by `oper2ist(P)`; [lin2lm](../indexing/lin2lm.md) maps the linear indices to spherical-tensor L,M labels. `coeffs` contains the corresponding coefficients returned by `oper2ist`. This function applies no additional multiplier or normalisation: the coefficients are exactly the output of that conversion for `P`.

## Inputs and checks

`mult` is intended to be a positive integer and `lvl_num` a level in the range 1 through `mult`. The source checks each for numeric, scalar, real values and checks the stated range for `lvl_num`, but does not explicitly test integrality or finiteness. `particle` must be the character value `'S'` or `'B'`; other values raise an error.
