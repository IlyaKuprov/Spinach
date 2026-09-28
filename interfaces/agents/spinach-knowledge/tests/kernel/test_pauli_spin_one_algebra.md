# tests/kernel/test_pauli_spin_one_algebra.m

- Signature: `result=test_pauli_spin_one_algebra()`

## Purpose

Tests spin-one angular momentum matrices.

## Physical / mathematical content

The spin-one representation has `Sz` projections `+1`, `0`, and `-1`, ladder matrix elements of `sqrt(2)`, and `S^2=s(s+1)=2`. The operators must obey the `su(2)` commutation relation `[Sx,Sy]=iSz`.

## Numerical / algorithmic content

Generates the spin-one operators with `pauli(3)` and compares them with explicit textbook matrices for `Sz`, `S+`, `Sx`, and `Sy`. It also checks the commutator and the Casimir operator `S^2` against `2*S.u`. Each comparison uses absolute and relative tolerances of `1e-15`.

## Outputs

- `result` - regression test result with explanatory messages.

## Implementation structure

- Announce the test target and create the regression test result.
- Generate Spinach spin-one operators and construct the textbook matrices.
- Check matrix elements, Cartesian operators, and the commutator.
- Check the Casimir operator.