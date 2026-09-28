# tests/kernel/test_orientation_average_suite.m

- Signature: `result=test_orientation_average_suite()`

## Purpose

Regression-tests two exact limiting cases for `orientation()` and `average()`.

## Physical / mathematical content

At zero Euler angles, the Wigner rotation is the identity, so the orientation contraction retains the diagonal rotational components. With zero positive- and negative-frequency Hamiltonian components, first-order averaging returns the unmodulated component `H0`.

## Numerical / algorithmic content

The test builds a synthetic rank-one rotational basis of sparse 2×2 matrices and compares `orientation(Q,[0 0 0])` with the sum of its three diagonal components. It then constructs a quiet one-proton spin system, sets `Hp` and `Hm` to zero, and compares `average(spin_system,Hp,H0,Hm,2*pi*1000,'ah_first_order')` with `H0`. Both comparisons use absolute and relative tolerances of `1e-14`.

## Outputs

- `result` — regression-test result containing the outcomes and explanatory messages for the zero-Euler-orientation and unmodulated-average checks.