# examples/fundamentals/operator_tests/commutation_6.m

- Signature: `commutation_6()`
- Source: [examples/fundamentals/operator_tests/commutation_6.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/fundamentals/operator_tests/commutation_6.m)

## Purpose

Checks Pauli-operator SU(2) relations, central-transition (CT) commutators and products, and CT irreducible spherical tensor (IST) expansions.

## Operators and identities

The no-argument function uses the source-set tolerance `tol=1e-10`. For `pauli(mult)` at multiplicities `2, 3, 4, 7`, it checks the cyclic relations `[Sx,Sy]=i Sz`, `[Sy,Sz]=i Sx`, and `[Sz,Sx]=i Sy`, plus `S.p=S.x+i*S.y` and `S.m=S.x-i*S.y`.

For `centrans(mult,...)` at multiplicities `4, 6, 8`, it checks `[CTx,CTy]=i CTz`, `[CTz,CT+]=CT+`, `[CTz,CT-]=-CT-`, and `[CT+,CT-]=2 CTz`. It also checks `CT+*CT-=CTz+P/2`, where the source constructs `P` as a diagonal support projector with ones at indices `mult/2` and `mult/2+1`. A final nested loop over those three multiplicities and CT types `x,y,z,+,-` obtains states and coefficients from `ct2ist`, reconstructs each matrix from `irr_sph_ten(mult,mult)`, and compares it with the original CT operator.

## Calling and numerical checks

Call `commutation_6()` in MATLAB with the Spinach functions used by the source available. Each residual is a Frobenius norm; the operator identities and reconstructed IST expansions are compared with `tol`. The source prints `Pauli operator commutation test PASSED.`, `CT commutation test PASSED.`, and `CT IST expansion test PASSED.` after their respective loops if their checks remain within tolerance; violations raise the corresponding error. The function reports no residual values.

Coverage is limited to the listed multiplicities, identities, and five CT types. It does not sweep arbitrary dimensions or establish behaviour for other representations or tensor-conversion choices.
