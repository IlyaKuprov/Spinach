# examples/fundamentals/convention_tests/tensors.m

- Signature: `tensors()`

## Purpose

Checks conversion of Stevens operator coefficients to irreducible spherical tensor (IST) coefficients by comparing the resulting operator matrices.

## Physical / mathematical content

For a spin-15/2 system, the test covers tensor ranks 1–6. At each rank it builds the same operator once from Stevens operators and once from ISTs, then compares the two matrices.

## Numerical / algorithmic content

Random real coefficient vectors are generated for each rank. The test passes when the 1-norm of the matrix difference is below `1e-6`; otherwise it displays both matrices and raises an error.

## Implementation structure

The code forms the Stevens combination, converts each coefficient vector with `stev2sph`, assembles the IST combination from `irr_sph_ten` operators, and checks the norm of their difference. The spin multiplicity is 16, corresponding to spin 15/2.
