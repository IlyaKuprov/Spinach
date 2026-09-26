# tests/kernel/test_operator_basis_suite.m

- Signature: `result=test_operator_basis_suite()`

## Purpose

Tests operator-basis construction and expansion helpers. Syntax: result=test_operator_basis_suite()

## Physical / mathematical content

## Numerical / algorithmic content

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
