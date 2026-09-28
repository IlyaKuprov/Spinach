# examples/fundamentals/operator_tests/commutation_7.m

- Signature: `commutation_7()`

## Purpose

Tests operator expansions between Zeeman, bosonic, IST, single-transition, and bosonic-monomial bases.

## Physical / mathematical content

The script reconstructs each of five Zeeman-level projectors from IST coefficients, each of six bosonic-level projectors from bosonic-monomial (BM) coefficients, and three bosonic products (`CA`, `ACCA`, and `CCAAA`) from IST coefficients. It also expands a random complex 7-by-7 matrix in the single-transition basis and a random complex bosonic matrix in the BM basis.

## Numerical / algorithmic content

The fixed reconstruction tolerance is `1e-10` for the projector and product expansions and the single-transition basis. The random bosonic matrix is checked using relative Frobenius error with threshold `1e-8`. Failed reconstructions raise errors.

## Implementation structure

Projector coefficients come from `enlev2ist` and `enlev2bm`; product coefficients come from `bos2ist`; and the random bosonic matrix uses `oper2bm`. Each expansion is reconstructed in its corresponding operator basis and compared with the original matrix.
