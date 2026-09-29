# tests/kernel/test_exponential_krylov_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_exponential_krylov_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_exponential_krylov_suite.m)

## Purpose

Regression test suite for the exponential, Chebyshev, and Krylov numerical utilities in Spinach. The suite verifies Arnoldi basis identities, Chebyshev coefficients, exponential drop boundary values, and Van Loan exponential-integral helpers against small closed-form references.

## Behaviour

- Announces the test target with `fprintf` and initialises a regression test result via `new_test_result` under the identifier `kernel/exponential_krylov_suite`, describing the requirement that Krylov bases and exponential-integral helpers reproduce closed-form matrix identities.
- **Arnoldi identities:** builds `krylov_mat=diag([1 2 4])` with the operator `@(x)krylov_mat*x` and calls `arnoldi(krylov_op,[1;2;3],2)`. Checks that `V'*V` equals the identity (tolerances `1e-13`) and that `krylov_mat*V(:,1:2)` equals `V*H` (tolerances `1e-13`), i.e. `A*V(:,1:n)` equals `V*H` for the extended Hessenberg matrix.
- **Arnoldi exact breakdown:** uses `break_mat=[5 0;0 7]` with initial vector `[1;0]` and `4` requested steps. Checks that `size(V_break,2)==1` with `size(H_break)` equal to `[1 1]` (an invariant one-dimensional Krylov subspace returned without zero padding) and that `H_break` equals `5` (tolerances `1e-15`), the eigenvalue on the invariant subspace.
- **Chebyshev coefficients:** defines `cheb_poly=@(x)2-3*x+4*(2*x.^2-1)` and calls `cheb_coeff(cheb_poly,-1,1,8)`. Checks the result against the reference `[2 -3 4 0 0 0 0 0]` (tolerances `1e-13`), since a degree-two Chebyshev polynomial has only its first three coefficients.
- **Exponential drop:** calls `expdrop(5,2,0.4,5,3)` and compares against the boundary-value closed form `5-fall_scale+fall_scale*exp(-3*fall_time)` with `fall_time=linspace(0,0.4,5)` and `fall_scale=(5-2)/(1-exp(-3*0.4))` (tolerances `1e-14`). Also checks that `all(diff(fall_obs)<0)`, i.e. a strictly monotonic fall for a positive drop rate from a larger value to a smaller value.
- **Exponential integral:** constructs a minimal spin system via `local_spin_system` and calls `expmint(spin_system,left_mat,mid_mat,right_mat,int_time)` with `left_mat=diag([1 2])`, `mid_mat=[1 2;3 4]/10`, `right_mat=diag([1/2 -3/2])`, and `int_time=0.35`. Compares against `local_expmint_ref` (tolerances `1e-12`), where each diagonal-matrix element integrates to a scalar complex exponential quotient.
- **Zero-time shortcut:** calls `expmint` with integration time `0` and checks the result equals the zero matrix of the same size as `mid_mat` (tolerances `1e-15`).
- **Nested exponential integral:** calls `expmint2(spin_system,0,2,0,3,0,nested_time)` with `nested_time=0.25` and compares against `2*3*nested_time^2/2` (tolerances `1e-14`); with zero generators the nested integral reduces to `B*D*T^2/2`.

## Inputs and outputs

- **Syntax:** `result=test_exponential_krylov_suite()`
- **Outputs:**
  - `result` — regression test result with explanatory messages, accumulated through `test_close` and `test_true` checks.
- The function takes no inputs.

## References

- Van Loan, C. F. — exponential-integral helper functions referenced in the suite header.
- Arnoldi iteration — basis for the Krylov subspace orthonormality and Hessenberg recurrence checks.
- Chebyshev polynomial expansions — basis for the `cheb_coeff` coefficient checks.
