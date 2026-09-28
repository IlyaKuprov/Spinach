# tests/kernel/test_wigner_angular_suite.m

- Signature: `result=test_wigner_angular_suite()`

## Purpose

Tests Clebsch-Gordan coefficients, Wigner symbols and matrices, and spherical harmonics against elementary values.

## Physical / mathematical content

- The spin-half coupling checks give the `|1,0>` triplet and `|0,0>` singlet coefficients for `|alpha beta>` as `1/sqrt(2)` in the tested Spinach convention; forbidden projections return zero.
- The suite checks the `(1 1 0; 0 0 0)` Wigner 3j value `-1/sqrt(3)`, two Wigner 6j values, the rank-one Wigner matrix at `beta=pi/2`, identity at zero Euler angles, and Wigner-matrix unitarity.
- It checks `Y_0^0=1/sqrt(4*pi)` and `Y_1^0=sqrt(3/(4*pi))*cos(theta)`.

## Numerical / algorithmic content

- The test evaluates each helper and compares its result with the specified scalar or matrix reference using `test_close`; unitarity is checked against the 5-by-5 identity matrix.

## Outputs

- `result` - regression test result with explanatory messages.

## Implementation structure

- Creates the regression result for `kernel/wigner_angular_suite`, then checks Clebsch-Gordan coefficients, Wigner 3j and 6j symbols, Wigner D matrices, and spherical harmonics.
