# test_clexpn_ensemble

Syntax: `result=test_clexpn_ensemble()`

Registered native contract tests for noncommuting ensemble estimators, full-order density, physical partitions, and mask/error propagation.

## Additional contract assertions

Checks exact raw-sample standard-error arithmetic, NaN for truncated internal averaging and single draws, unchanged signals under a task ceiling, legacy-equivalent interlacing density errors, and a noncommuting bath requiring adaptive group merging before held-out validation. The adaptive complete-bath result is compared bitwise with an explicitly supplied complete group.
