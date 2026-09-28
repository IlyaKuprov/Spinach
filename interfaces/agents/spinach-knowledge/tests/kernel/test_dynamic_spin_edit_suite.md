# tests/kernel/test_dynamic_spin_edit_suite.m

- Signature: `result=test_dynamic_spin_edit_suite()`

## Purpose

Regression tests for deterministic spin-system editing utilities.

## Tests

- `kill_spin` checks numeric and logical spin removal, including updated spin identities, labels, coordinates, particle types, isotope hash, coupling/proximity dimensions, chemical-part indices, and relaxation-source indices. It also checks removal across multiple parts, including an empty part, and deletion of stale basis, connectivity, symmetry, and assumption data.
- Checks dilute single-spin subsystem construction and numeric/isotope `dictum` overrides.
- Checks `merge_inp` on system fields, coordinates, common values, square blocks, rate arrays, index offsets, susceptibility centres, column-oriented parts, and partless chemistry; it also checks rejection of mismatched common values, partial nested groups, unknown fields, and missing isotope lists.

## Outputs

- result -regression test result with explanatory messages
- The test checks spin removal, dilute-isotope subsystem generation,
- assumption overrides, and merging of Spinach input structures.
