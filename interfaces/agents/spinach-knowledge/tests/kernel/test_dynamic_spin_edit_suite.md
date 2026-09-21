# tests/kernel/test_dynamic_spin_edit_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_dynamic_spin_edit_suite.m`
- Signature: `result=test_dynamic_spin_edit_suite()`
- Total lines: 270

## Purpose

Tests deterministic spin-system editing support utilities. Syntax: result=test_dynamic_spin_edit_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file also defines local helper function(s): `local_merge_fails()`, `local_edit_spin_system()`, `local_assumed_spin_system()`, `local_merge_parts()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Outputs

- result -regression test result with explanatory messages
- The test checks spin removal, dilute-isotope subsystem generation,
- assumption overrides, and merging of Spinach input structures.

## Implementation structure

- Tests deterministic spin-system editing support utilities. Syntax:
- result=test_dynamic_spin_edit_suite()
- result -regression test result with explanatory messages
- The test checks spin removal, dilute-isotope subsystem generation,
- assumption overrides, and merging of Spinach input structures.
- Announce the test target
- State the spin editing target of the test
- Build a small but structurally complete spin-system descriptor
- Check numeric spin removal updates identities, labels, coordinates, and parts
- Check chemical part renumbering across multiple subsystems
- Check destruction of stale basis, connectivity, symmetry, and assumption data
- Check logical spin removal follows the same path

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `local_edit_spin_system()`, `kill_spin()`, `test_true()`, `isequal()`, `strcmp()`, `md5_hash()`, `logical()`, `speye()`, `true()`, `isfield()`, `dilute()`, `cellfun()`, `nnz()`, `all()`, `local_assumed_spin_system()`.
