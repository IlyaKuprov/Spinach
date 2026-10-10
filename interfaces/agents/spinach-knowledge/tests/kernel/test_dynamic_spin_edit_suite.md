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

The basis-bearing removal fixture is built by `test_spin_system`; removal must
rebuild its 16-state two-spin basis while clearing stale symmetry and assumptions.

Longitudinal and zero-quantum isotope filters are tested after partial isotope removal, loss of the last local match while another substance retains it, removal of all matching isotopes globally, and creation of a spin-free substance. Rebuilt descriptors, offsets, and hashes must equal an explicit rebuild with the surviving numeric and isotope filters.

IK-0, IK-1, and IK-2 basis-bearing fixtures test local depth clipping through both `kill_spin` and `dilute`. The affected three-spin substance shrinks to two spins, while an unaffected three-spin substance retains its smaller depth of two. Both `prox_level` and the `space_level` alias are exercised for IK-1/IK-2; retained depth cells, descriptors, offsets, and hashes must match explicitly capped basis settings.

IK-DNP and IK-SBS vector-depth fixtures separately remove an electron/mode or a nucleus/spin while retaining both particle classes. The resulting depth vectors, descriptors, and hashes must match a rebuild with explicitly bounded settings.

Last-class fixtures remove all electrons or nuclei from IK-DNP, and all modes
or spins from IK-SBS. Each surviving substance must match an explicit IK-0
rebuild at its surviving class depth, including descriptors, offsets, and hash;
a separate substance must retain its descriptor and settings.

## Reaction-record fixtures

Synthetic editing objects carry an explicit empty reaction list. Merge tests retain empty chemistry, spin offsets, column-oriented part lists, and all non-chemical assertions. An explicit selector-loss record tests simultaneous substance and selector-spin offsets; with declared parts, correlation times concatenate per substance rather than using the retired partless radical-pair fields.
