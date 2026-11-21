# tests/kernel/test_mecce_rates.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_rates.m)

## Purpose

Amplitude damping, dephasing, column convention, and T1 limitation.

## Behaviour

Checks the following named assertions in the source:

- `C7 numerical check`.
- `C7 contract check`.
- `Zero-rate limit`.
- `Evolving exterior refusal`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_mecce_rates()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_rates.m)

## Native contract

The registered test calls clexpn with free_evolution and native compile/subsystem/generator services. Natural density answers are time-by-column-vectorised observed matrices. It preserves amplitude damping, longitudinal dephasing, weak-rate propagation and error coverage, reciprocal operator/rate scaling, negative-rate refusal, exact zero-rate equivalence, active exterior-jump refusal, the population-product counterexample, nonnormal loss leakage, and independent three-level master-equation references. Zero-rate projected jumps disappear before manifold checks. The weak-rate stock-propagation assertion uses the native expert services with its original cleanup thresholds.
