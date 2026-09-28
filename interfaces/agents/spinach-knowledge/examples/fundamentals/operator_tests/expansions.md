# examples/fundamentals/operator_tests/expansions.m

- Signature: `expansions()`

## Purpose

Checks the left and right product-action tables for irreducible spherical tensors and orthogonalised bosonic monomials.

## Physical / mathematical content

The test chooses random spin multiplicity `1+randi(9)` and bosonic level count `2+randi(8)`. It obtains IST and orthogonalised bosonic-monomial bases and their respective product-action tables, then compares table-derived actions with direct operator multiplication for every pair of basis elements.

## Numerical / algorithmic content

For both bases, left and right actions are checked separately. The difference is scaled by the Frobenius norm of the acted-on basis element, and the test fails if either relative discrepancy exceeds `sqrt(eps)`.

## Implementation structure

Nested loops over basis-index pairs reconstruct left and right products from the corresponding table coefficients. After all pairs pass, the script reports separate success messages for the IST and bosonic-monomial tables.
