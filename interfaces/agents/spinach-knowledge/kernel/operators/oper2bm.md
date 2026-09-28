# kernel/operators/oper2bm.m

- Signature: `[states,coeffs]=oper2bm(A)`

## Purpose

Expands a square operator into bosonic monomials.

## Physical / mathematical content

The routine represents `A` in the bosonic-monomial basis for a mode whose level-space dimension is `size(A,1)`. Returned `states` use the Spinach BM basis indexing convention; use `lin2kq` to convert those indices to K,Q bosonic-monomial indices.

## Numerical / algorithmic content

The basis is generated with `boson_mono(size(A,1))`. The routine computes the basis-overlap matrix using `hdot`, obtains the coefficients by solving the overlap system, then drops terms whose coefficient magnitude is at most `10*eps('double')`.

## Parameters / inputs

- A - numeric square matrix to expand.

## Outputs

- states - BM basis indices for the terms retained in the expansion.
- coeffs - coefficients of the corresponding bosonic monomials in the linear combination.

## Implementation structure

1. Check that `A` is a numeric square matrix.
2. Generate the BM basis and its overlap matrix, then solve for the expansion coefficients.
3. Return the matching basis indices and coefficients after removing negligible terms.
