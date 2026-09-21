# tests/kernel/test_operator_basis_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_operator_basis_suite.m`
- Signature: `result=test_operator_basis_suite()`
- Total lines: 141

## Purpose

Tests operator-basis construction and expansion helpers. Syntax: result=test_operator_basis_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `ist_reconstruct()`, `bm_reconstruct()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks single-spin tensor bases, bosonic bases, central
- transitions, Stevens operators, single-transition matrices, expansion
- helpers, and sparse preallocation dimensions.

## Implementation structure

- Tests operator-basis construction and expansion helpers. Syntax:
- result=test_operator_basis_suite()
- result -regression test result with explanatory messages
- The test checks single-spin tensor bases, bosonic bases, central
- transitions, Stevens operators, single-transition matrices, expansion
- helpers, and sparse preallocation dimensions.
- Announce the test target
- State the operator-basis target of the test
- Irreducible spherical tensors obey [Lz,T(k,m)]=m*T(k,m)
- Stevens rank-one zero-projection operator is Lz
- Weyl boson operators obey their defining number-operator commutators
- Bosonic monomials and orthogonalised monomials must have documented structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `pauli()`, `irr_sph_ten()`, `test_close()`, `int2str()`, `projections()`, `comm()`, `stevens()`, `weyl()`, `boson_mono()`, `speye()`, `boson_ortho()`, `hdot()`, `sin_tran()`, `lin2kq()`, `centrans()`.
