# tests/kernel/test_mecce_exact.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_exact.m)

## Purpose

Zero jumps, dense Liouville full bath and correlated support.

## Behaviour

Checks the following named assertions in the source:

- `C7 numerical check`.
- `ME-gCCE Hahn zero jumps`.
- `ME-CCE Hahn zero jumps`.
- `Full-bath projected ME reference`.
- `Dissipative Hahn reference`.
- `C7 contract check`.
- `PyCCE reconstructed ME-CCE2`.
- `PyCCE reconstructed ME-gCCE2`.
- `ME tutorial validity`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_mecce_exact()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_mecce_exact.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)

Supplied-spin_config regression checks compare gCCE and ME-gCCE density, signal and observables (raw and normalized) with an independent diagonal Hamiltonian solution. They also verify projected bath phases and swapped user projections, and require informative rejection of unsupported spin_config state/projection/reference aliases.

## Native contract

Native migration is blocked by the selective instantaneous bath-operator pulse tested at fraction 0.4. `free_evolution` accepts finite events only; the native call with this ideal event is rejected with `free_evolution:sequence`. Named Hahn sequences do not represent this arbitrary off-midpoint selective pulse. The test is retained unchanged; no physical assertion is removed.
