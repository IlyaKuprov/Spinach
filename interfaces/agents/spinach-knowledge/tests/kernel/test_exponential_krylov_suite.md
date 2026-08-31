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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Exponential and Krylov numerical utilities\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/exponential_krylov_suite', 'Exponential and Krylov numerical utilities', 'Krylov bases and exponential-integral helpers must reproduce clo…`.
- Lines 25-26: Check the Arnoldi orthogonality and projected recurrence identities; implemented by `krylov_mat=diag([1 2 4])`.
- Lines 36-37: Check exact Arnoldi breakdown from an eigenvector initial condition; implemented by `break_mat=[5 0;0 7]`.
- Lines 46-47: Check Chebyshev coefficients for a quadratic polynomial in T0,T1,T2; implemented by `cheb_poly=@(x)2-3*x+4*(2*x.^2-1)`.
- Lines 54-55: Check the exponential drop against the boundary-value closed form; implemented by `fall_obs=expdrop(5,2,0.4,5,3)`.
- Lines 65-66: Check the single exponential integral against a diagonal closed form; implemented by `spin_system=local_spin_system()`.
- Lines 77-78: Check the zero-time shortcut in the exponential integral helper; implemented by `int_zero=expmint(spin_system,left_mat,mid_mat,right_mat,0)`.
- Lines 83-84: Check the nested exponential integral in the scalar nilpotent case; implemented by `nested_time=0.25`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/exponential_krylov_suite', 'Exponential and Krylov numerical utilities', 'Krylov bases and exponential-integral helpers must reproduce clo…`.
- Lines 26: computes `krylov_mat` using `krylov_mat=diag([1 2 4])`.
- Lines 27: computes `krylov_op` using `krylov_op=@(x)krylov_mat*x`.
- Lines 28: computes `[V,H]` using `[V,H]=arnoldi(krylov_op,[1;2;3],2)`.
- Lines 37: computes `break_mat` using `break_mat=[5 0;0 7]`.
- Lines 38: computes `[V_break,H_break]` using `[V_break,H_break]=arnoldi(@(x)break_mat*x,[1;0],4)`.
- Lines 47: computes `cheb_poly` using `cheb_poly=@(x)2-3*x+4*(2*x.^2-1)`.
- Lines 48: computes `cheb_obs` using `cheb_obs=cheb_coeff(cheb_poly,-1,1,8)`.
- Lines 49: computes `cheb_ref` using `cheb_ref=[2 -3 4 0 0 0 0 0]`.
- Lines 55: computes `fall_obs` using `fall_obs=expdrop(5,2,0.4,5,3)`.
- Lines 56: computes `fall_time` using `fall_time=linspace(0,0.4,5)`.
- Lines 57: computes `fall_scale` using `fall_scale=(5-2)/(1-exp(-3*0.4))`.
- Lines 58: computes `fall_ref` using `fall_ref=5-fall_scale+fall_scale*exp(-3*fall_time)`.
- Lines 66: computes `spin_system` using `spin_system=local_spin_system()`.
- Lines 67: computes `left_mat` using `left_mat=diag([1 2])`.
- Lines 68: computes `mid_mat` using `mid_mat=[1 2;3 4]/10`.
- Lines 69: computes `right_mat` using `right_mat=diag([1/2 -3/2])`.
- Lines 70: computes `int_time` using `int_time=0.35`.

### Local helper functions

- Line 94: `local_spin_system()` — `function spin_system=local_spin_system()`. Closed-form integral for diagonal A and C matrices
  - Representative operation: `spin_system.sys.enable={}`.
  - Representative operation: `spin_system.sys.disable={}`.
- Line 106: `local_expmint_ref()` — `function int_ref=local_expmint_ref(left_mat,mid_mat,right_mat,int_time)`.
  - Representative operation: `int_ref=zeros(size(mid_mat))`.
  - Representative operation: `left_freq=diag(left_mat)`.

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
