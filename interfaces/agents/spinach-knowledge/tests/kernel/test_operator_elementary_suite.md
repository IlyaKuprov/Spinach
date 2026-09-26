# tests/kernel/test_operator_elementary_suite.m

- Signature: `result=test_operator_elementary_suite()`

## Purpose

Tests elementary operator generators in kernel/operators. Syntax: result=test_operator_elementary_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks analytic commutation relations, indexing conventions,
- and small explicit matrices for low-level operator constructors.

## Implementation structure

- Tests elementary operator generators in kernel/operators. Syntax:
- result=test_operator_elementary_suite()
- result -regression test result with explanatory messages
- The test checks analytic commutation relations, indexing conventions,
- and small explicit matrices for low-level operator constructors.
- Announce the test target
- State the operator target of the test
- Check spin-one angular momentum commutation and ladder definitions
- Check finite-truncation Weyl algebra away from the unavoidable edge state
- Check bosonic monomial serpentine indexing
- Check Gram-Schmidt orthogonality without imposing normalisation
- Check single-transition basis indexing from the documented 4x4 map
