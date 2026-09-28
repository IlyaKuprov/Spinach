# kernel/operators/oper2ist.m

- Signature: `[states,coeffs]=oper2ist(A)`

## Purpose

Expands a square operator in the irreducible spherical tensor basis for a single spin.

## Physical / mathematical content

The tensor basis is generated for a spin multiplicity of `size(A,1)`. Returned `states` follow the Spinach IST basis indexing convention; use `lin2lm` to convert them to L,M spherical-tensor indices.

## Numerical / algorithmic content

For each tensor `X` from `irr_sph_ten(size(A,1))`, the coefficient is computed as `hdot(X,A)/hdot(X,X)` and converted to a full value. Terms with coefficient magnitude at most `10*eps('double')` are removed.

## Parameters / inputs

- A - numeric square matrix to expand.

## Outputs

- states - IST basis indices corresponding to retained terms.
- coeffs - coefficients of the corresponding irreducible spherical tensors in the expansion.

## Implementation structure

1. Check that `A` is a numeric square matrix.
2. Generate the single-spin irreducible spherical tensor basis with `irr_sph_ten`.
3. Compute each normalized inner-product coefficient, then return only terms above the numerical cutoff with their basis indices.
