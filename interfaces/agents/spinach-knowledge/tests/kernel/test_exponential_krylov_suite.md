# tests/kernel/test_exponential_krylov_suite.m

- Signature: `result=test_exponential_krylov_suite()`

## Purpose

Tests exponential, Chebyshev, and Krylov numerical utilities. Syntax: result=test_exponential_krylov_suite()

## Physical / mathematical content

- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Outputs

- result -regression test result with explanatory messages
- The test checks Arnoldi basis identities, Chebyshev coefficients,
- exponential drop boundary values, and Van Loan exponential-integral
- helpers against small closed-form references.

## Implementation structure

- Tests exponential, Chebyshev, and Krylov numerical utilities. Syntax:
- result=test_exponential_krylov_suite()
- result -regression test result with explanatory messages
- The test checks Arnoldi basis identities, Chebyshev coefficients,
- exponential drop boundary values, and Van Loan exponential-integral
- helpers against small closed-form references.
- Announce the test target
- State the utility target of the test
- Check the Arnoldi orthogonality and projected recurrence identities
- Check exact Arnoldi breakdown from an eigenvector initial condition
- Check Chebyshev coefficients for a quadratic polynomial in T0,T1,T2
- Check the exponential drop against the boundary-value closed form
