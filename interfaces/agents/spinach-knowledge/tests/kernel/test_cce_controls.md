# tests/kernel/test_cce_controls.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_controls.m)

## Purpose

Check finite drives, waveform refinement and named pulse schedules.

## Behaviour

Checks the following named assertions in the source:

- `Finite detuned drive matches independent dense expm`.
- `Public finite-control driver`.
- `Finite pulse approaches ideal rotation`.
- `Piecewise midpoint refinement is second order`.
- `XY8 timestamps`.
- `XY8 phases`.
- `UDD timestamps`.
- `Exterior ideal pulse flips the sampled field`.
- `Exterior pulse changes subsequent Hamiltonian`.
- `Finite exterior drive is refused rather than frozen`.
- `Conditional finite central drive refused`.
- `Selective register operator ordering`.
- `Physical spin-one pi pulse`.
- `Selective two-level pi pulse`.
- `Physical and selective pulses differ`.

The test uses `new_test_result`, `test_true` and/or `test_close`; numerical tolerances and independent reference calculations are defined in the linked source.

## Inputs and outputs

```matlab
result=test_cce_controls()
```

**Inputs**

- None; run with Spinach and its test libraries on the MATLAB path.

**Outputs**

- `result` — standard test-harness result with checks, pass/fail status and diagnostics. Printed numerical metrics are not additional result fields.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_cce_controls.m)
- [Test manifest](https://github.com/IlyaKuprov/Spinach/blob/main/tests/lib/test_manifest.m)
