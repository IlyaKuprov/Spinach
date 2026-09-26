# tests/kernel/test_pauli_spin_one_algebra.m

- Signature: `result=test_pauli_spin_one_algebra()`

## Purpose

Tests spin-one angular momentum matrices. Syntax: result=test_pauli_spin_one_algebra()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks the spin-one representation: Sz projections are +1, 0,
- and -1, ladder matrix elements are sqrt(2), and S^2=s(s+1)=2.

## Implementation structure

- Tests spin-one angular momentum matrices. Syntax:
- result=test_pauli_spin_one_algebra()
- result -regression test result with explanatory messages
- The test checks the spin-one representation: Sz projections are +1, 0,
- and -1, ladder matrix elements are sqrt(2), and S^2=s(s+1)=2.
- Announce the test target
- State the physical target of the test
- Generate Spinach spin-one operators
- Write the textbook spin-one matrices explicitly
- Check matrix elements and commutators
- Check the Casimir operator
