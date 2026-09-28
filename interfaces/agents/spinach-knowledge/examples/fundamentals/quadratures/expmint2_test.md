# examples/fundamentals/quadratures/expmint2_test.m

- Signature: `expmint2_test()`

## Purpose

Verify expmint2 against a reference formed by nested numerical integration of matrix-valued integrands.

## Physical / mathematical content

This is a numerical test, not a spin-system simulation: it uses random matrices, with A, C, and E made Hermitian, and compares a double-exponential integral evaluated by Spinach with the MATLAB integral-based construction.

## Numerical / algorithmic content

The random matrix dimension is 11–15 and the upper integration limit is randomly chosen between 1 and 11. The Frobenius-norm difference is tested against the source's 10*n*eps('double') tolerance ratio, producing a pass or fail message.

## Implementation structure

- Generate the matrices and bootstrap Spinach.
- Call expmint2(spin_system,A,B,C,D,E,ul).
- Build the nested reference using MATLAB's array-valued integral calls and compare the two results.
