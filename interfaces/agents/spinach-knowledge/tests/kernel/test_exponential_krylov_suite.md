# tests/kernel/test_exponential_krylov_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_exponential_krylov_suite.m`
- Signature: `result=test_exponential_krylov_suite()`
- Total lines: 125

## Purpose

Tests exponential, Chebyshev, and Krylov numerical utilities. Syntax: result=test_exponential_krylov_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.
- The file also defines local helper function(s): `local_spin_system()`, `local_expmint_ref()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `arnoldi()`, `test_close()`, `test_true()`, `isequal()`, `cheb_coeff()`, `expdrop()`, `all()`, `diff()`, `local_spin_system()`, `expmint()`, `local_expmint_ref()`, `expmint2()`, `right_freq()`, `left_freq()`, `int_ref()`.
