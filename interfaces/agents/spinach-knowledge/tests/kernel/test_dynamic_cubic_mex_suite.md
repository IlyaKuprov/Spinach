# tests/kernel/test_dynamic_cubic_mex_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_cubic_mex_suite.m`
- Signature: `result=test_dynamic_cubic_mex_suite()`
- Total lines: 103

## Purpose

Tests the cubic-polynomial MEX helper used by eigenfields(). Syntax: result=test_dynamic_cubic_mex_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Cubic polynomial root MEX helper\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/dynamic_cubic_mex_suite', 'Cubic polynomial root MEX helper', 'eigenfields cubic root helper must match analytical and Matlab references.')`.
- Lines 25-26: Set the production root tolerance; implemented by `root_tol=sqrt(eps)`.
- Lines 28-32: Check three roots, including endpoints; implemented by `result=test_close(result,'cubic_roots three endpoint roots', cubic_roots([1 -1.5 0.5 0],root_tol),[0 0.5 1], 1e-12,1e-12, 'x*(x-1/2)*(x-1) must return all three unit-int…`.
- Lines 34-38: Check a triple root; implemented by `result=test_close(result,'cubic_roots triple root', cubic_roots([1 -1.5 0.75 -0.125],root_tol),0.5, 1e-9,1e-12, 'a cubic with a triple root should return one merged root…`.
- Lines 40-44: Check a double root plus an endpoint root; implemented by `result=test_close(result,'cubic_roots double root', cubic_roots([1 -1.4 0.49 0],root_tol),[0 0.7], 1e-9,1e-12, 'x*(x-0.7)^2 should return the endpoint and merged double…`.
- Lines 46-50: Check quadratic and linear degeneracies; implemented by `result=test_close(result,'cubic_roots quadratic degeneracy', cubic_roots([0 1 -1 0],root_tol),[0 1], 1e-12,1e-12, 'leading zero cubic coefficient should reduce to the qu…`.
- Lines 56-58: Check constant and zero polynomials; implemented by `result=test_true(result,'cubic_roots constant polynomial',isempty(cubic_roots([0 0 0 1],root_tol)), 'a non-zero constant polynomial has no roots')`.
- Lines 62-66: Check extreme coefficient scaling; implemented by `result=test_close(result,'cubic_roots large coefficient scale', cubic_roots([1e300 -1.5e300 0.5e300 0],root_tol),[0 0.5 1], 1e-12,1e-12, 'normalisation should prevent ov…`.
- Lines 72-73: Check derivative-root use case; implemented by `turn_ref=sort((3+[-1 1]*sqrt(3))/6)`.
- Lines 79-80: Check a random set of well separated roots; implemented by `rng_state=rng`.

### Control flow inferred from the code

- Line 84: `for` loop over `n=1:200`.
- Line 86: conditional branch on `min(diff(test_roots))<1e-3`.
- Line 90: conditional branch on `numel(obs_roots)~=3`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/dynamic_cubic_mex_suite', 'Cubic polynomial root MEX helper', 'eigenfields cubic root helper must match analytical and Matlab references.')`.
- Lines 26: computes `root_tol` using `root_tol=sqrt(eps)`.
- Lines 73: computes `turn_ref` using `turn_ref=sort((3+[-1 1]*sqrt(3))/6)`.
- Lines 80: computes `rng_state` using `rng_state=rng`.
- Lines 81: computes `rng_cleanup` using `rng_cleanup=onCleanup(@()rng(rng_state))`.
- Lines 83: computes `rand_ok` using `rand_ok=true; rand_err=0`.
- Lines 85: computes `test_roots` using `test_roots=sort(0.05+0.90*rand(1,3))`.
- Lines 89: computes `obs_roots` using `obs_roots=cubic_roots(poly(test_roots),root_tol)`.
- Lines 94: computes `rand_err` using `rand_err=max(rand_err,max(abs(obs_roots-test_roots)))`.

## Outputs

- result -regression test result with explanatory messages
- The test checks cubic roots, degenerate lower-order polynomials,
- repeated roots, endpoint roots, extreme coefficient scaling, and
- derivative-root use cases against explicit references.

## Implementation structure

- Tests the cubic-polynomial MEX helper used by eigenfields(). Syntax:
- result=test_dynamic_cubic_mex_suite()
- result -regression test result with explanatory messages
- The test checks cubic roots, degenerate lower-order polynomials,
- repeated roots, endpoint roots, extreme coefficient scaling, and
- derivative-root use cases against explicit references.
- Announce the test target
- State the utility target of the test
- Set the production root tolerance
- Check three roots, including endpoints
- Check a triple root
- Check a double root plus an endpoint root

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `cubic_roots()`, `test_true()`, `onCleanup()`, `rng()`, `diff()`, `poly()`, `clear()`.
