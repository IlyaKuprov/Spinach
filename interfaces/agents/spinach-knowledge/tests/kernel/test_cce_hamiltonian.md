# tests/kernel/test_cce_hamiltonian.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_hamiltonian.m)

## Purpose

Check arbitrary spin, ordering, rotations, retention and fallback.

## Behaviour

Checks the following named assertions in the source:

- `dense propagator parity`.
- `transpose mutation rejected`.
- `nonsymmetric half cells`.
- `physical spin permutation`.
- `multiple central spins`.
- `central-only subset`.
- `exterior expectation`.
- `physical distance cutoff`.
- `cutoff omission count`.
- `memory ceiling rejected`.
- `unsupported assumption rejected`.
- `rank-three fallback`.
- `fallback recorded`.
- `unsupported direct rejection`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_hamiltonian()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_hamiltonian.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
