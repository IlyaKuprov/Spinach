# tests/kernel/test_orientation_average_suite.m

- Signature: `result=test_orientation_average_suite()`

## Purpose

Tests orientation() and average() on small exact cases. Syntax: result=test_orientation_average_suite()

## Physical / mathematical content

- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks the zero-Euler-angle orientation contraction against an
- explicit diagonal Wigner sum and verifies that first-order average
- Hamiltonian theory leaves an unmodulated Hamiltonian unchanged.

## Implementation structure

- Tests orientation() and average() on small exact cases. Syntax:
- result=test_orientation_average_suite()
- result -regression test result with explanatory messages
- The test checks the zero-Euler-angle orientation contraction against an
- explicit diagonal Wigner sum and verifies that first-order average
- Hamiltonian theory leaves an unmodulated Hamiltonian unchanged.
- Announce the test target
- State the rotational-kernel target of the test
- Build a synthetic rank-one rotational basis
- Contract the zero-orientation Hamiltonian
- Check the zero-angle Wigner identity path
- Build a quiet spin system for average() diagnostics
