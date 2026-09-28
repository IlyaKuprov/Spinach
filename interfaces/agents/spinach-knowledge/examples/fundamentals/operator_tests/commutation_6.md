# examples/fundamentals/operator_tests/commutation_6.m

- Signature: `commutation_6()`

## Purpose

Checks Pauli and central-transition operator identities, including central-transition IST expansions.

## Physical / mathematical content

For Pauli operators at multiplicities 2, 3, 4, and 7, the test verifies the cyclic SU(2) commutators and the definitions of the raising and lowering operators. For central-transition operators at multiplicities 4, 6, and 8, it checks the transverse commutator, both `z`-ladder commutators, `[CT+,CT-]=2CTz`, and `CT+CT-=CTz+P/2`, where `P` projects onto the two central levels. It also expands each central-transition Cartesian or ladder operator (`x`, `y`, `z`, `+`, `-`) into irreducible spherical tensor (IST) terms and reconstructs the matrix.

## Numerical / algorithmic content

Commutator, product, ladder-definition, and IST reconstruction discrepancies are compared with `tol=1e-10`. A violation raises the associated error; each test family prints a pass message after its checks.

## Implementation structure

The script first loops over Pauli multiplicities, then over central-transition multiplicities to check operator identities and build the central-level projector. A final nested loop over central-transition multiplicities and operator types obtains coefficients from `ct2ist`, reconstructs with `irr_sph_ten`, and checks the Frobenius-norm discrepancy.
