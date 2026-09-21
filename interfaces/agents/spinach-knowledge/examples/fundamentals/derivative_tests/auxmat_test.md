# examples/fundamentals/derivative_tests/auxmat_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/auxmat_test.m`
- Signature: `auxmat_test()`
- Total lines: 61

## Purpose

Testing IK's favourite equation numerically -auxiliary matrix expression against a high-accuracy finite diffe- rence approximation.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Testing IK's favourite equation numerically -auxiliary
- matrix expression against a high-accuracy finite diffe-
- rence approximation.
- Random matrices, arbitrary function
- Some weird ass function and its derivatives
- Real coefficients and fin. diff. increment
- First parameter
- Second parameter
- Third parameter

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `randi()`, `dL_da()`, `dL_db()`, `dL_dc()`.
