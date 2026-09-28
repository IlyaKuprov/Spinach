# kernel/operators/enlev2bm.m

- Signature: `[states,coeffs]=enlev2bm(nlevels,lvl_num)`

## Purpose

Expands the projector onto one energy level of a truncated bosonic mode as a linear combination of bosonic monomials.

## Physical / mathematical content

The selected level is numbered upward from the empty-mode state. The routine constructs the diagonal projector with its only nonzero element at (lvl_num,lvl_num), then expands that operator in the bosonic-monomial basis.

## Numerical / algorithmic content

A zero matrix of size nlevels by nlevels is created, the selected diagonal element is set to one, and oper2bm converts the projector to bosonic-monomial states and coefficients. Use lin2kq to convert the returned basis indices to K,Q indices when needed.

## Parameters / inputs

- nlevels - number of energy levels in the mode; a positive integer.
- lvl_num - energy-level index, counting upward from the empty-mode state and bounded by nlevels.

## Outputs

- states - bosonic-monomial basis indices contributing to the operator; use lin2kq to convert them to K,Q indices.
- coeffs - coefficients of those bosonic monomials in the linear combination.

## Implementation structure

1. Check that nlevels and lvl_num are numeric, scalar, real values and that the level is within the specified range.
2. Build the diagonal level projector.
3. Expand it with oper2bm and return the resulting states and coefficients.
