# tests/kernel/test_operator_elementary_suite.m

- Signature: `result=test_operator_elementary_suite()`

## Purpose

Regression test for elementary operator generators in `kernel/operators`. It checks algebraic identities, basis indexing, orthogonality, and small explicit operator matrices.

## Checks

- For `pauli(3)`, verifies the spin-one commutator `[Sx,Sy]=i Sz`, the raising-operator definition `S+=Sx+i Sy`, and the identity operator.
- For `weyl(4)`, verifies `C A=N`, `[N,C]=C`, and `[N,A]=-A`.
- For `boson_mono(3)`, verifies that serpentine indices 1–3 give the identity, creation, and annihilation operators. For `boson_ortho(3)`, verifies zero off-diagonal Hilbert–Schmidt overlaps without requiring normalisation.
- For `sin_tran(4)`, verifies that indices 10 and 7 give the `(1,4)` and `(4,1)` single-transition matrices, respectively.
- For `centrans(4,...)`, compares the central-transition `z`, `+`, and `-` operators with explicit matrices acting on the middle two states.
- For `irr_sph_ten(3,2)`, checks `[Lz,T(2,m)]=m T(2,m)` for `m=2,1,0,-1,-2`.
- For `stevens(3,...)`, checks that the rank-zero operator is the identity and that the rank-two components with `q=+1` and `q=-1` are Hermitian.

## Output

`result` is a regression-test result containing explanatory messages for the checks.