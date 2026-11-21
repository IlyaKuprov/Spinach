# tests/kernel/test_gcce_exact.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_exact.m)

## Purpose

Four bath spins, complex pure/mixed states, and physical Hahn pulses, tested through native clexpn calls and independent dense Kronecker/expm references.

## Behaviour

Checks the following named assertions in the source:

- `raw full-observed-spin system exact`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_gcce_exact()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_gcce_exact.m)

Native preparation and observation are tested with independent receiver scales from 1e-20 to 1e20, true zero-reference rejection, and both row and column detector cells. The same 64 row/column comparisons cover free, directly exponentiated Hahn, finite-drive, master-equation, and raw interlaced evaluations; physical inputs and assertion tolerances are fixed.

The y-phase Hahn rotation is specified by explicit transition columns alpha=[1;0], beta=[0;1i], with the physical initial density supplied separately. Direct raw interlacing uses the registered callback and a context-generated ensemble, not a legacy evaluator. The context call used solely to obtain that ensemble has order zero; the subsequent raw calculation uses the requested physical cluster.

Native-only regressions reject unknown or Hilbert-only preferences on master-equation requests, and reject direct-segment Hahn propagation for projected or Liouville representations.

Finite-control invocations use finite_drive with native physical descriptors, not a raw sequence field. Independent references retain their physical inputs and assertions; unused raw-event declarations are not part of the reference construction.
