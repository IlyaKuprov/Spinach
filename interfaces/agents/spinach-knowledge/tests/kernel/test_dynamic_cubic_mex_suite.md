# tests/kernel/test_dynamic_cubic_mex_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_cubic_mex_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_cubic_mex_suite.m)

## Purpose

Regression test suite for the cubic-polynomial root MEX helper used by `eigenfields()`. The suite verifies that the helper's root extraction matches explicit analytical and MATLAB references across standard, degenerate, repeated, endpoint, extreme-scaling, and derivative-root cases.

## Behaviour

- Announces the test target with `fprintf('TESTING: Cubic polynomial root MEX helper\n')`.
- Initialises a test result object via `new_test_result('kernel/dynamic_cubic_mex_suite', ...)`, describing the requirement that the eigenfields cubic root helper match analytical and MATLAB references.
- Sets the production root tolerance to `root_tol = sqrt(eps)`.
- Checks three endpoint roots for `x*(x-1/2)*(x-1)` (coefficients `[1 -1.5 0.5 0]`) against `[0 0.5 1]` with tolerances `1e-12` (absolute and relative).
- Checks a triple root case `[1 -1.5 0.75 -0.125]` against `0.5` with tolerances `1e-9` and `1e-12`.
- Checks a double root plus endpoint case `x*(x-0.7)^2` (coefficients `[1 -1.4 0.49 0]`) against `[0 0.7]` with tolerances `1e-9` and `1e-12`.
- Checks quadratic degeneracy `[0 1 -1 0]` against `[0 1]` and linear degeneracy `[0 0 2 -1]` against `0.5`, both with tolerances `1e-12`.
- Verifies that constant polynomial `[0 0 0 1]` and zero polynomial `[0 0 0 0]` return empty root sets.
- Checks extreme coefficient scaling: large scale `[1e300 -1.5e300 0.5e300 0]` and small scale `[1e-300 -1.5e-300 0.5e-300 0]`, both expected to return `[0 0.5 1]` with tolerances `1e-12`, confirming that normalisation prevents overflow and underflow.
- Checks the derivative-root use case `[0 3 -3 0.5]` against `turn_ref = sort((3+[-1 1]*sqrt(3))/6)` with tolerances `1e-12`, expecting both turning points.
- Runs a randomised check over 200 trials: saves and restores the RNG state via `onCleanup`, seeds with `rng(1)`, generates three sorted roots in `[0.05, 0.95]`, skips sets with minimum root separation below `1e-3`, builds coefficients with `poly`, and requires exactly three returned roots with maximum absolute error below `1e-10`.
- Uses `test_close` for numerical comparisons and `test_true` for boolean conditions, accumulating results and explanatory messages into the returned test result object.

## Inputs and outputs

**Syntax:**

```matlab
result = test_dynamic_cubic_mex_suite()
```

**Outputs:**

- `result` — regression test result object with explanatory messages.

**Inputs:**

- None.

## References

- `cubic_roots` — cubic-polynomial root MEX helper under test.
- `eigenfields` — production module that uses the tested helper.
- `new_test_result`, `test_close`, `test_true` — test framework utilities.
