# tests/kernel/test_cce_checkpoint.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_checkpoint.m)

## Purpose

Partial ledger resume, identity rejection, and write integrity.

## Behaviour

Checks the following named assertions in the source:

- `partial raw restart`.
- `partial reduced restart`.
- `completed bitmap`.
- `time identity`.
- `model identity`.
- `seed identity`.
- `cluster identity`.
- `raw corruption`.
- `completed task immutability`.
- `failed replacement integrity`.
- `write failure integrity`.
- `actual save failure integrity`.
- `execution-neutral restart`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_checkpoint()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status, and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_checkpoint.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

## Native experiment coverage

The registered test enters through clexpn with free_evolution, finite_drive, and noise_correlation. It preserves interrupted raw/reduced replay, adaptive-partition pilot ownership, tensor/tolerance/seed/time identities, checksum and completed-record integrity, actual filesystem save failure, execution-neutral restart, and fresh-versus-cold physical-model checks. Projected/generalised Hilbert and master-equation requests are explicit native approximation and representation choices. Direct callbacks and interlacing refuse checkpoint ownership; an independent convergence request cannot share a ledger.

The real save-failure test denies file creation in its owned temporary directory: POSIX owner write permission on Unix, and a temporary NTFS creation-deny rule on Windows. Windows directory read-only attributes alone do not refuse file creation. The rule is removed before loading and comparing the known-good checkpoint; no production checkpoint behaviour or numerical tolerance is changed.
