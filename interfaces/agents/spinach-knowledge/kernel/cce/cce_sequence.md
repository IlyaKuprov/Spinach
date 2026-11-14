# kernel/cce/cce_sequence.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_sequence.m)

## Purpose

Compile whole-sequence sample times and ideal conditional pulse events.

## Behaviour

- Optional public fields are checked by `cce_validate`; unknown, reserved, and method-inapplicable requests are rejected rather than silently ignored.

- Positive absolute timestamps are accepted: each sample includes only events at or before its duration, while the zero sample always has no pulses.

- Checks finite, nonnegative, sorted row-vector sample times in seconds and explicit event ordering.
- Ideal events use `type=ideal`, physical `target` IDs, axis, angle in rad, and `placement=fraction|absolute`; `at` is a fraction of each total sequence time or seconds, respectively.
- Conditional CCE accepts only complete-register ideal pi swaps. Named sequences, finite pulses and DEER requests are routed through `cce_controls`.

- Explicitly permitted empty DEER selections are no-op events. After selector validation, empty-target pulses are removed before conditional/ME/noise schedule recursion; generalized control compilation already treats them as no-ops.

- Explicit pulse-generator Hermiticity is checked relative to its norm, with no absolute unit floor; reciprocal operator/angle rescaling cannot conceal a non-Hermitian generator.

## Inputs and outputs

```matlab
schedule=cce_sequence(model,parameters)
```

**Inputs**

- `model` — prepared register and bath IDs.
- `parameters` — `times`, `sequence` (`[]` for FID), method and central-state declarations.

**Outputs**

- `schedule` — sample times, compiled events and conditional/generalized control metadata.

## References

- [Source file on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/cce/cce_sequence.m)
- [cce_controls.m](cce_controls.md)
- [cce_cluster.m](cce_cluster.md)

An empty event list remains valid when an enabled DEER selector has no pulse to select; the conditional finite-event check does not dereference an empty numeric array. Zero-valued unused padding fields in mixed event arrays are harmless; nonzero unused control components are rejected.
