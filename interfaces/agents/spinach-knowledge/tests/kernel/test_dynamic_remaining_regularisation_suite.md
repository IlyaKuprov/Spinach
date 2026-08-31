# tests/kernel/test_dynamic_remaining_regularisation_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_remaining_regularisation_suite.m`
- Signature: `result=test_dynamic_remaining_regularisation_suite()`
- Total lines: 62

## Purpose

Tests remaining regularisation and inverse-problem utilities. Syntax: result=test_dynamic_remaining_regularisation_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The code contains an inverse-problem or ill-conditioning aspect and therefore introduces explicit regularisation, model selection, or stabilisation logic.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Remaining regularisation utilities\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/dynamic_remaining_regularisation_suite', 'Remaining regularisation utilities', 'Regularisation helpers must recover compact analytical inv…`.
- Lines 25-26: Check L-curve analysis on a synthetic corner near lambda equals one; implemented by `lam=logspace(-3,3,9)`.
- Lines 34-35: Check positivity-constrained Tikhonov inversion against a scalar analytic solution; implemented by `K=1`.
- Lines 50-51: Check L1 sparsity targeting on an identity sensing matrix; implemented by `A=eye(3)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_remaining_regularisation_suite', 'Remaining regularisation utilities', 'Regularisation helpers must recover compact analytical inv…`.
- Lines 26: computes `lam` using `lam=logspace(-3,3,9)`.
- Lines 27: computes `err` using `err=sqrt(1+lam.^2)`.
- Lines 28: computes `reg` using `reg=sqrt(1+lam.^-2)`.
- Lines 29: computes `lam_opt` using `lam_opt=lcurve(lam,err,reg,'log')`.
- Lines 35: computes `K` using `K=1`.
- Lines 36: computes `D` using `D=1`.
- Lines 37: computes `KtK` using `KtK=1`.
- Lines 38: computes `DtD` using `DtD=1`.
- Lines 39: computes `H` using `H=4`.
- Lines 40: computes `y` using `y=2`.
- Lines 41: computes `lambda` using `lambda=1`.
- Lines 42: computes `[x_tikh,err_tikh,reg_tikh]` using `[x_tikh,err_tikh,reg_tikh]=tikhonov(K,D,KtK,DtD,H,y,lambda)`.
- Lines 46: computes `'the residual error at x` using `'the residual error at x=1 is (1-2)^2')`.
- Lines 48: computes `'the regularisation signal at x` using `'the regularisation signal at x=1 is x^2')`.
- Lines 51: computes `A` using `A=eye(3)`.
- Lines 53: computes `[x_l1,err_l1,reg_l1]` using `[x_l1,err_l1,reg_l1]=tikhol1n(A,y,1)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks L-curve corner detection, positivity-constrained
- Tikhonov inversion, and L1 sparsity targeting on compact analytical
- inverse problems.

## Implementation structure

- Tests remaining regularisation and inverse-problem utilities. Syntax:
- result=test_dynamic_remaining_regularisation_suite()
- result -regression test result with explanatory messages
- The test checks L-curve corner detection, positivity-constrained
- Tikhonov inversion, and L1 sparsity targeting on compact analytical
- inverse problems.
- Announce the test target
- State the utility target of the test
- Check L-curve analysis on a synthetic corner near lambda equals one
- Check positivity-constrained Tikhonov inversion against a scalar analytic solution
- Check L1 sparsity targeting on an identity sensing matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `lcurve()`, `test_true()`, `tikhonov()`, `test_close()`, `minimising()`, `tikhol1n()`, `nnz()`, `x_l1()`.
