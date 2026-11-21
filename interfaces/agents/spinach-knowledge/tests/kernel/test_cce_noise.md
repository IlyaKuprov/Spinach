# tests/kernel/test_cce_noise.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_noise.m)

## Purpose

Full order, disconnected cumulants, ordering, and commuting limit.

## Behaviour

Checks the following named assertions in the source:

- `C7 numerical check`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_noise()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_noise.m)

Checks sparse unconnected noise for independent polarised sites, time-dependent ordered/symmetrised moments, and connected covariance of an explicitly constructed finite sampled mixture under both averaging estimators.

Public optional-field contracts are explicit in reused fixtures: irrelevant method, event, geometry, and legacy numerical options are removed rather than passed silently.

Direct noise spin_configs with nontracial bath states are checked against dense ordered and symmetrised connected correlations; the initial moments come from the supplied bath state.

## Native contract

The registered test calls clexpn with noise_correlation, retaining both ensemble averaging orders, supplied nontracial bath states, ordered and symmetrised dense references, disconnected cross terms, finite-mixture covariance, allocation preflight, and unavailable internal connected-covariance uncertainty. Direct registered evaluator calls reject foreign experiment requests. The exact scaled operator 1e-20*[0 1;0 0] remains a non-Hermitian rejection case. Numeric correlation arrays and separate diagnostics replace wrapped results without changing physical assertions.
