# tests/kernel/test_dynamic_spin_edit_suite.m

Source: [tests/kernel/test_dynamic_spin_edit_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_spin_edit_suite.m)

## Purpose

Regression test suite for the deterministic spin-system editing support utilities in the Spinach kernel: `kill_spin`, `dilute`, `dictum`, and `merge_inp`. The suite verifies that local spin-system editing helpers update dependent metadata without touching external state.

## What the suite checks

- **Removing a spin:** `kill_spin` preserves consistent isotope/label identity, coordinate and chemical-partition dimensions, and renumbers spin-dependent indices such as SRSK sources. Obsolete basis and interaction metadata cannot survive a changed spin list. A numeric spin index and the equivalent logical deletion mask must yield the same reduced system.
- **Isotopic dilution and assumptions:** `dilute` constructs separate subsystems with a single selected low-abundance `13C` at each labelled site. `dictum` changes the intended Zeeman or coupling strength through numeric spin/pair selectors or an isotope selector without changing unrelated assignments.
- **Merging input structures:** `merge_inp` concatenates extensive spin fields, block-merges pair matrices, offsets spin and chemical-subsystem indices, and retains common settings such as the shared temperature. It must reject conflicting common values, partially specified nested groups, unknown fields and missing isotope lists rather than silently producing inconsistent combined inputs.

## Inputs and outputs

```matlab
result=test_dynamic_spin_edit_suite()
```

- **Output**: `result` — regression result reporting which spin-editing invariants pass and any failed assertions.
- The function takes no inputs.

## References

- Tested functions: `kill_spin`, `dilute`, `dictum`, `merge_inp`.
- Test infrastructure: `new_test_result`, `test_true`.
- Support utilities referenced by the fixture: `spin`, `md5_hash`.
