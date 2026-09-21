# examples/fundamentals/quadratures/expmint2_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/quadratures/expmint2_test.m`
- Signature: `expmint2_test()`
- Total lines: 48

## Purpose

Verification of expmint2 against numerical integration.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Verification of expmint2 against numerical
- integration.
- Random dimension and upper limit
- Generate random matrices
- Bootstrap Spinach
- Call Spinach function
- Inner integrand and its integral
- Outer integrand and its integral
- Call Matlab integrator
- Test the error norm

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `randi()`, `bootstrap()`, `expmint2()`, `integral()`, `int_inner()`, `int_outer()`, `eps()`.
