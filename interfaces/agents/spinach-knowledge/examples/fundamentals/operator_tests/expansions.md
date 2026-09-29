# examples/fundamentals/operator_tests/expansions.m

- Signature: `expansions()`

## Purpose

Checks the left- and right-product-action tables for irreducible spherical tensors (ISTs) and orthogonalised bosonic monomials against direct matrix multiplication.

## Mathematical content

For each ordered pair of basis elements, the script uses the table coefficients to reconstruct the product action and compares that matrix with the corresponding direct product. It therefore tests consistency of `ist_product_table` and `bos_product_table` with the matrix bases constructed in the same invocation. It is an algebraic consistency check, not a physical-system calculation.

## Callable context and model

Call the zero-input MATLAB function `expansions()` from a Spinach checkout with the project functions on the MATLAB path. It returns no values and reports success separately for the IST and bosonic-monomial groups; a failed comparison raises an error. The function selects `spin_mult=1+randi(9)` (a multiplicity from 2 through 10) and `bos_nlevels=2+randi(8)` (a level count from 3 through 10), then creates `irr_sph_ten(spin_mult)`, `boson_ortho(bos_nlevels)`, and their respective product-action tables.

## Checks encoded in the source

For every pair of IST basis indices, the left table is reconstructed from `ist_PTL` and the right table from `ist_PTR`; each table uses basis matrices divided by their Frobenius norms. The corresponding direct products are `T{n}*T{m}/norm(T{m},'fro')` and `T{m}*T{n}/norm(T{m},'fro')`. For each side, the Frobenius-norm difference is divided by `norm(T{m},'fro')` and must not exceed `sqrt(eps)`. The bosonic group performs the analogous comparisons using `bos_PTL`, `bos_PTR`, and `boson_ortho`'s matrices, with the same tolerance rule based on `norm(B{m},'fro')`.

## Assumptions and limits

The dimension choices are random and no seed is set in the function, so the selected dimensions can vary between calls. For the chosen dimensions, the loops cover every ordered basis-index pair; the checks compare the generated tables with direct finite-matrix products, not with an independent implementation or a physical observable. The source contains tolerances and conditional pass messages, not fixed residual values, and this description does not claim a run passed.

## Source

[`examples/fundamentals/operator_tests/expansions.m`](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/expansions.m)
