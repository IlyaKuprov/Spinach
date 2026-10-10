# tests/kernel/test_dynamic_remaining_regularisation_suite.m

## Purpose

Regression test for the remaining regularisation and inverse-problem utilities in Spinach. The suite checks L-curve corner detection, positivity-constrained Tikhonov inversion, and L1 sparsity targeting on compact analytical inverse problems.

## Behaviour

The function announces the test target with `fprintf('TESTING: Remaining regularisation utilities\n')` and initialises a test result via `new_test_result` for `kernel/dynamic_remaining_regularisation_suite`, describing the target as "Remaining regularisation utilities" with the requirement that regularisation helpers must recover compact analytical inverse-problem references.

Three checks are performed:

1. **L-curve analysis on a synthetic corner near lambda equals one.** The test builds `lam=logspace(-3,3,9)`, `err=sqrt(1+lam.^2)` and `reg=sqrt(1+lam.^-2)`, then calls `lcurve(lam,err,reg,'log')`. After `close all`, `test_true` asserts that the returned optimum `lam_opt` satisfies `lam_opt>0.1 && lam_opt<10`, with the message that the maximum-curvature point of the symmetric synthetic L-curve should lie near lambda=1.

2. **Positivity-constrained Tikhonov inversion against a scalar analytic solution.** With `K=1`, `D=1`, `KtK=1`, `DtD=1`, `H=4`, `y=2` and `lambda=1`, the test calls `[x_tikh,err_tikh,reg_tikh]=tikhonov(K,D,KtK,DtD,H,y,lambda)` and uses `test_close` (tolerances `1e-8`, `1e-8`) to verify:
   - `x_tikh` equals 1, since minimising `(x-2)^2+x^2` with `x>=0` gives `x=1`;
   - `err_tikh` equals 1, the residual error at `x=1` being `(1-2)^2`;
   - `reg_tikh` equals 1, the regularisation signal at `x=1` being `x^2`.

3. **L1 sparsity targeting on an identity sensing matrix.** With `A=eye(3)` and `y=[2;0;0]`, the test calls `[x_l1,err_l1,reg_l1]=tikhol1n(A,y,1)` and asserts via `test_true` that `nnz(abs(x_l1)>1e-8)==1 && x_l1(1)>0` (identity data with one non-zero target should keep only the populated component) and that `err_l1` and `reg_l1` are finite with `err_l1>=0` and `reg_l1>0` (finite non-negative error and positive L1 regularisation metrics).

## Inputs and outputs

**Syntax:**

```matlab
result = test_dynamic_remaining_regularisation_suite()
```

The function takes no inputs.

**Outputs:**

- `result` — regression test result structure with explanatory messages, accumulated through `test_true` and `test_close` assertions.

## References

- Source: [tests/kernel/test_dynamic_remaining_regularisation_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_remaining_regularisation_suite.m)

Functions exercised by the suite: `new_test_result`, `lcurve`, `tikhonov`, `tikhol1n`, `test_true`, `test_close`.
