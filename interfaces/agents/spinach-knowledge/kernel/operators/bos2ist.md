# kernel/operators/bos2ist.m

- Signature: `[states,coeffs]=bos2ist(prod_spec,nlevels)`

## Purpose

Convert a bosonic operator product specification into its contributing irreducible spherical-tensor (IST) states and coefficients.

## Physical / mathematical content

`prod_spec` specifies an ordered product of operators for a truncated bosonic mode: `C` is creation, `A` is annihilation, and `N` is the number operator. The allowed symbols are `C`, `A`, and `N`; an empty specification represents the identity operator.

## Numerical / algorithmic content

The function starts with an `nlevels`-by-`nlevels` sparse identity matrix, obtains the Weyl operators from `weyl(nlevels)`, and multiplies in the operators named by `prod_spec` in order. It converts the resulting matrix with `oper2ist`, which supplies `states` and `coeffs`.

## Parameters / inputs

- `prod_spec` - character specification of the ordered product, using `C` (creation), `A` (annihilation), and `N` (number); the empty character vector yields the identity.
- `nlevels` - number of energy levels in the truncated bosonic mode. The explicit check requires a numeric, real scalar at least 1; it does not test integrality, although the error message describes a positive real integer.

## Outputs

- `states` - contributing states in Spinach IST basis indexing; use `lin2lm` to convert to spherical-tensor `L,M` indices.
- `coeffs` - coefficients of the corresponding ISTs in the linear combination.

## Implementation structure

The routine checks that `prod_spec` is a character value and that each character is one of `C`, `A`, or `N`. It initializes the sparse identity, applies the selected Weyl operators sequentially, then calls `oper2ist` to produce the output state indices and coefficients.
